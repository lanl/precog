// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
#include <iostream>
#include <memory>
#include <cstdio>
#include <unordered_set>
#include <vector>
#include <cstdint>
#include <cmath>
#include <cstring>
#ifdef _OPENMP
  #include <omp.h>
#endif
#include "nanoflann.hpp"

// [[Rcpp::export]]
void eoe_omp_set_threads(int n) {
#ifdef _OPENMP
  omp_set_num_threads(n);
#endif
}

// [[Rcpp::export]]
int eoe_omp_max_threads() {
#ifdef _OPENMP
  return omp_get_max_threads();
#else
  return -1;
#endif
}

// [[Rcpp::export]]
int eoe_omp_threads_in_parallel() {
#ifdef _OPENMP
  int nthreads = 1;
  #pragma omp parallel
  {
    #pragma omp single
    nthreads = omp_get_num_threads();
  }
  return nthreads;
#else
  return -1;
#endif
}

// ========================== Helper functions =================================
// scale to 0-1 ----------------------------------------------------------------
//[[Rcpp::export]]
arma::mat scale_cpp(arma::mat M) {
  arma::uword n = M.n_rows;
  arma::rowvec l = arma::min(M,0);
  arma::rowvec h = arma::max(M,0);
  arma::vec Ones = arma::ones(n);
  M = (M-Ones*l) / (Ones*(h-l));
  return M;
}
// find unique rows ------------------------------------------------------------
template <typename T>
inline bool rows_equal(const T& lhs, const T& rhs, double tol=1e-12) {
  return arma::approx_equal(lhs, rhs, "absdiff", tol);
}
arma::mat unique_rows(const arma::mat& x) {
  unsigned int count = 1, i = 1, j = 1, nr = x.n_rows, nc = x.n_cols;
  arma::mat result(nr, nc);
  result.row(0) = x.row(0);
  
  for ( ; i < nr; i++) {
    bool matched = false;
    if (rows_equal(x.row(i), result.row(0))) continue;
    
    for (j = i + 1; j < nr; j++) {
      if (rows_equal(x.row(i), x.row(j))) {
        matched = true;
        break;
      }
    }
    
    if (!matched) result.row(count++) = x.row(i);
  }
  
  return result.rows(0, count-1);
}

// Hash combine helper
static inline std::size_t hash_combine(std::size_t seed, std::size_t h) {
  return seed ^ (h + 0x9e3779b97f4a7c15ULL + (seed<<6) + (seed>>2));
}

struct RowKey {
  std::vector<long long> v;  // quantized coordinates
};

struct RowKeyHash {
  std::size_t operator()(RowKey const& k) const noexcept {
    std::size_t seed = 0;
    for (auto &x : k.v) {
      seed = hash_combine(seed, std::hash<long long>{}(x));
    }
    return seed;
  }
};

struct RowKeyEq {
  bool operator()(RowKey const& a, RowKey const& b) const noexcept {
    return a.v == b.v;
  }
};

arma::mat unique_rows_tol_hash(const arma::mat& x, double eps) {
  arma::uword nr = x.n_rows, nc = x.n_cols;
  arma::mat out(nr, nc);
  arma::uword count = 0;

  std::unordered_set<RowKey, RowKeyHash, RowKeyEq> seen;
  seen.reserve(static_cast<std::size_t>(nr * 1.3));

  for (arma::uword i = 0; i < nr; ++i) {
    RowKey key;
    key.v.resize(nc);
    for (arma::uword j = 0; j < nc; ++j) {
      // Quantize to grid of width eps
      key.v[j] = static_cast<long long>(std::llround(x(i, j) / eps));
    }
    if (seen.emplace(std::move(key)).second) {
      out.row(count++) = x.row(i);
    }
  }
  return out.rows(0, count - 1);
}



// find a vector in a matrix ---------------------------------------------------
arma::uword compare_v_m (arma::rowvec v, arma::mat M){
  // find the first row that match v
  arma::uword idx = 0;
  for(arma::uword i=0; i<M.n_rows; i++){
    if(arma::all(v==M.row(i))){
      idx = i;
      break;
    }
  }
  return (idx);
}
// knn -------------------------------------------------------------------------
class DF
{
private:
  std::shared_ptr<arma::mat> df_;
  
public:
  void import_data(arma::mat& df)
  {
    df_ = std::make_shared<arma::mat>(df.t());
  }
  
  unsigned int kdtree_get_point_count() const
  {
    return df_->n_cols ;
  }
  
  double kdtree_get_pt(const unsigned int idx, const unsigned int dim) const 
  {
    return (*df_)(dim, idx);
  }
  
  const double* get_row(const unsigned int idx) const
  {
    return &(*df_)(0, idx);
  }
  
  template <class BBOX>
  bool kdtree_get_bbox(BBOX&) const 
  { 
    return false; 
  }
};
typedef nanoflann::KDTreeSingleIndexDynamicAdaptor<nanoflann::L2_Adaptor<double, DF>, DF, -1, unsigned int> kdTree;
class KDTree
{
private:
  const unsigned int dim_; // dimension of data
  // const unsigned int N_; //data size
  const unsigned int n_; //query size
  DF data_;
  DF query_;
  const unsigned int k_;
  
public:
  KDTree(arma::mat& data, arma::mat& query, int k) : 
  dim_(data.n_cols), n_(query.n_rows), k_(k)
  {
    if(static_cast<unsigned int>(query.n_cols) != dim_)
      Rcpp::Rcerr << "\nDimensions do not match.\n";
    else
    {
      data_.import_data(data);
      query_.import_data(query);
    }
  }
  
  Rcpp::List knn_cpp(bool knnx_flag=true)
  {
    kdTree tree(dim_, data_, nanoflann::KDTreeSingleIndexAdaptorParams(2*k_));
    
    arma::umat nn_index(n_,k_);
    arma::mat nn_dist(n_,k_);
    
    nanoflann::KNNResultSet<double> resultSet(k_);
    
    for(unsigned int i = 0; i<n_; i++)
    { 
      std::size_t *index = new std::size_t[k_];
      double *distance = new double[k_];
      resultSet.init(index, distance);
      tree.findNeighbors(resultSet, query_.get_row(i), nanoflann::SearchParams());
      for (unsigned int kth = 1; kth<=k_;kth++){
        nn_index(i,kth-1) = index[kth-1];
        nn_dist(i,kth-1) = std::sqrt(distance[kth-1]);
      } 
      delete[] index;
      delete[] distance;
    }
    Rcpp::List L;
    if (knnx_flag){
      L = Rcpp::List::create(Rcpp::Named("nn_index")=nn_index,Rcpp::Named("nn_dist")=nn_dist);
    }else{
      L = Rcpp::List::create(Rcpp::Named("nn_index")=nn_index.cols(1,k_-1),Rcpp::Named("nn_dist")=nn_dist.cols(1,k_-1));
    }
    return L;
  }
};

//[[Rcpp::export]]
Rcpp::List knnx_my(arma::mat& data, arma::mat& points, int k)
{
  KDTree kdt(data, points,k);
  return kdt.knn_cpp(true);
}
//[[Rcpp::export]]
Rcpp::List knn_my(arma::mat& data, int k)
{
  KDTree kdt(data, data, k+1);
  return kdt.knn_cpp(false);
}

// generate random orthonormal matrix
arma::mat randortho_cpp(arma::uword n){
  if (n==1){
    return arma::ones(1,1);
  }
  arma::mat Z = arma::randn(n, n) / std::sqrt(2);
  arma::mat Q;
  arma::mat R;
  // QR decompostion
  arma::qr(Q, R, Z);
  arma::vec d = R.diag();
  arma::vec ph = d/arma::abs(d);
  return (Q*arma::diagmat(ph));
}
//  ============================== Main functions ==============================
// generate uniform random point in ball----------------------------------------
//[[Rcpp::export]]
arma::mat runif_in_sphere_cpp(arma::uword n, arma::uword p){
  arma::mat out(n, p);
  arma::rowvec Ones = arma::ones<arma::rowvec>(p);
  arma::mat Z = arma::randn(n, p);
  arma::vec U = arma::randu(n);
  arma::vec Z_norm(n);
  for (arma::uword i=0; i<n; i++){
    Z_norm(i) = arma::norm(Z.row(i));// the norm of the Z vector
  }
  Z = Z / (Z_norm * Ones);
  U = arma::pow(U, 1./p);
  out = (U * Ones) % Z;
  return (out);
}

// generate uniform points in a (manifold) ball --------------------------------
//[[Rcpp::export]]
arma::mat ball_gen(arma::rowvec cen, double rad, arma::uword n, arma::uword p=0,
                   Nullable<NumericMatrix> sub_=R_NilValue,
                   bool rand=true,
                   Nullable<NumericMatrix> twinsample_=R_NilValue){
  if(p==0){
    p = cen.n_cols;
  }
  arma::uword pp = cen.n_cols; //output dimension, p should be no larger than pp
  arma::mat y(n, p);
  arma::mat yy(n, pp);
  arma::mat M(n, pp); // returned matrix
  arma::mat sub(pp, pp); // vectors defining the subspace
  arma::vec Ones = arma::ones(n);
  
  if (rand){
    y = runif_in_sphere_cpp(n, p);
  }else{
    if (twinsample_.isNotNull()){
      NumericMatrix twinsample_tempt(twinsample_);       // casting to arma matrix
      y = as<arma::mat>(wrap(twinsample_tempt));
      n = y.n_rows;
      y = y * randortho_cpp(p);
    }else{
      // std::cout<<"No twinning sample provided!";
      Rcpp::Rcout << "No twinning sample provided!" << std::endl;
    }
  }
  // check if we need to generate on tangent space or not
  if(p == pp){
    M = Ones*cen + y*rad;
  }else{
    NumericMatrix sub_tempt(sub_);       // casting to arma matrix
    arma::mat subb = as<arma::mat>(wrap(sub_tempt));
    yy = arma::join_rows(y, arma::zeros(n) * arma::ones<arma::rowvec>(1));//y:nxp--->yy:nxpp
    sub = arma::join_rows(subb, arma::zeros(pp)*arma::ones<arma::rowvec>(pp-p)); //sub_:ppxp-->sub:ppxpp
    M = Ones*cen + yy*sub.t()*rad;
  }
  return (M);
}


// Approx points for output space ----------------------------------------------
//[[Rcpp::export]]
arma::mat approx_gen(arma::mat X, Rcpp::List d, arma::uword p, bool rand=true,
                     Nullable<NumericMatrix> twinsample_=R_NilValue){
  arma::uword n = X.n_rows;
  arma::uword q = X.n_cols; // dimension in the output space
  
  arma::umat d_nn_index = d["nn_index"];
  arma::mat d_nn_dist = d["nn_dist"];
  
  arma::uword no_nb = d_nn_index.n_cols;// number of neighbor; d is the knn between output points
  arma::uword no_sub = std::min(std::min(p,q), no_nb); //number of points in one subspace simplex
  arma::uword n_ball = 2*(std::min(p,q)+1) + 1 + no_nb; //number of points in one ball
  if (twinsample_.isNotNull()){
    NumericMatrix twinsample_tempt(twinsample_);
    n_ball = twinsample_tempt.nrow();
  }
  arma::uword N = no_nb*n + (no_sub+2)*n + n*n_ball; // total number of approx points
  
  arma::mat pca(1+no_sub, p); // pca coefficients
  
  arma::mat point_approx(N, q);
  arma::uword count = 0;
  
  // mid points
  for (arma::uword i=0; i<no_nb; i++){
    point_approx.rows(count, count+n-1) = (X + X.rows(d_nn_index.col(i))) / 2;
    count += n;
  }
  // simplex and ball
  double delta = 0.5;
  arma::vec rad = d_nn_dist.col(0);
  arma::mat tempt(1+no_sub, q); // points defining the subspace
  arma::uvec idx_axial(no_sub);
  
  for (arma::uword i=0; i<n; i++){
    tempt.row(0) = X.row(i);
    tempt.rows(1, no_sub) = X.rows(d_nn_index.row(i).cols(0, no_sub-1));
    // centroid
    point_approx.row(count) = arma::mean(tempt, 0);
    count += 1;
    
    // axial points
    idx_axial = arma::linspace<arma::uvec>(1, no_sub, no_sub);
    point_approx.row(count) = (1+delta) * arma::mean(tempt.rows(idx_axial),0) - delta*tempt.row(0);
    count += 1;
    for (arma::uword j=1; j<(no_sub); j++){
      idx_axial.rows(0, j-1) = arma::linspace<arma::uvec>(0, j-1, j);
      idx_axial.rows(j, no_sub-1) = arma::linspace<arma::uvec>(j+1, no_sub, no_sub-j);
      point_approx.row(count) = (1+delta) * arma::mean(tempt.rows(idx_axial),0) - delta*tempt.row(j);
      count += 1;
    }
    
    idx_axial = arma::linspace<arma::uvec>(0, no_sub-1, no_sub);
    point_approx.row(count) = (1+delta) * arma::mean(tempt.rows(idx_axial), 0) - delta*tempt.row(no_sub);
    count += 1;
    // ball
    if(p<q){
      pca = arma::princomp(tempt);
      Rcpp::NumericMatrix pca_rcpp = wrap(pca.cols(0,p-1));
      point_approx.rows(count,count + n_ball-1) = ball_gen(X.row(i),rad(i), n_ball, 
                        p, pca_rcpp, rand, twinsample_);
    }else{
      point_approx.rows(count,count+n_ball-1) = ball_gen(X.row(i), rad(i), n_ball, 
                        q, R_NilValue, rand,twinsample_);
    }
    count += n_ball;
  }
  
  return unique_rows_tol_hash(point_approx, 1e-7);
}





// // Approx points for output space ----------------------------------------------
// //[[Rcpp::export]]
// arma::mat approx_gen(arma::mat X, Rcpp::List d, arma::uword p, bool rand=true,
//                      Nullable<NumericMatrix> twinsample_=R_NilValue, int n_threads=1){
//   arma::uword n = X.n_rows;
//   arma::uword q = X.n_cols; // dimension in the output space
  
//   arma::umat d_nn_index = d["nn_index"];
//   arma::mat d_nn_dist = d["nn_dist"];
  
//   arma::uword no_nb = d_nn_index.n_cols;// number of neighbor; d is the knn between output points
//   arma::uword no_sub = std::min(std::min(p,q), no_nb); //number of points in one subspace simplex
//   arma::uword n_ball = 2*(std::min(p,q)+1) + 1 + no_nb; //number of points in one ball
//   if (twinsample_.isNotNull()){
//     NumericMatrix twinsample_tempt(twinsample_);
//     n_ball = twinsample_tempt.nrow();
//   }
//   arma::uword N = no_nb*n + (no_sub+2)*n + n*n_ball; // total number of approx points
//   arma::mat pca(1+no_sub, p); // pca coefficients
//   arma::mat point_approx(N, q);
//   arma::uword count = 0;
  
//   // mid points
//   for (arma::uword i=0; i<no_nb; i++){
//     point_approx.rows(count, count+n-1) = (X + X.rows(d_nn_index.col(i))) / 2;
//     count += n;
//   }
//   // simplex and ball
//   double delta = 0.5;
//   arma::vec rad = d_nn_dist.col(0);
//   arma::mat tempt(1+no_sub, q); // points defining the subspace
//   arma::uvec idx_axial(no_sub);
  
//   arma::uword base  = no_nb * n;
//   arma::uword block = (no_sub + 2) + n_ball;

//   // I need everything in openMP loop to work with arma (no Rcpp)
//   arma::mat twinsample_mat;
//   const arma::mat* twinsample_ptr = nullptr;

//   if (!rand && twinsample_.isNotNull()) {
//     Rcpp::NumericMatrix ts_nm(twinsample_);
//     twinsample_mat = Rcpp::as<arma::mat>(ts_nm);
//     twinsample_ptr = &twinsample_mat;
//   }

//   #ifdef _OPENMP
//   omp_set_num_threads(n_threads);
//   #endif
  
//   Rcpp::Rcout << "n_threads passed = " << n_threads << "\n";
//   Rcpp::Rcout.flush();

//   int threads_used = 1;
//   #ifdef _OPENMP
//   #pragma omp parallel
//   {
//     #pragma omp single
//     threads_used = omp_get_num_threads();
//   }
//   #endif

//   Rcpp::Rcout << "threads_used = " << threads_used << "\n";
//   Rcpp::Rcout.flush();


//   #ifdef _OPENMP
//   #pragma omp parallel for schedule(static)
//   #endif
//   for (arma::uword i = 0; i < n; i++) {

//     arma::uword off = base + i * block;
//     arma::uword c = off;

//     arma::mat tempt(1 + no_sub, q);     // thread-local
//     arma::uvec idx_axial(no_sub);       // thread-local

//     tempt.row(0) = X.row(i);
//     tempt.rows(1, no_sub) = X.rows(d_nn_index.row(i).cols(0, no_sub - 1));

//     // centroid
//     point_approx.row(c) = arma::mean(tempt, 0);
//     c += 1;

//     // axial points
//     idx_axial = arma::linspace<arma::uvec>(1, no_sub, no_sub);
//     point_approx.row(c) = (1 + delta) * arma::mean(tempt.rows(idx_axial), 0) - delta * tempt.row(0);
//     c += 1;

//     for (arma::uword j = 1; j < no_sub; j++) {
//       idx_axial.rows(0, j - 1) = arma::linspace<arma::uvec>(0, j - 1, j);
//       idx_axial.rows(j, no_sub - 1) = arma::linspace<arma::uvec>(j + 1, no_sub, no_sub - j);
//       point_approx.row(c) = (1 + delta) * arma::mean(tempt.rows(idx_axial), 0) - delta * tempt.row(j);
//       c += 1;
//     }

//     idx_axial = arma::linspace<arma::uvec>(0, no_sub - 1, no_sub);
//     point_approx.row(c) = (1 + delta) * arma::mean(tempt.rows(idx_axial), 0) - delta * tempt.row(no_sub);
//     c += 1;

//     arma::rowvec center = X.row(i);
//     double r = rad(i);

//     if (p < q) {
//       arma::mat pca   = arma::princomp(tempt);
//       arma::mat basis = pca.cols(0, p - 1);
//       point_approx.rows(c, c + n_ball - 1) =
//         ball_gen_arma(center, r, n_ball, p, &basis, rand, twinsample_ptr);
//     } else {
//       point_approx.rows(c, c + n_ball - 1) =
//         ball_gen_arma(center, r, n_ball, q, nullptr, rand, twinsample_ptr);
//     }
//     c += n_ball;  // <-- outside the if/else


//   }

//     Rcpp::Rcout << "we made it out of the loop! \n";
//     Rcpp::Rcout.flush();
  
//   return unique_rows_tol_hash(point_approx, 1e-6);
// }


// Fill distance ---------------------------------------------------------------
//[[Rcpp::export]]
arma::mat filldist_cpp(arma::mat M, arma::uword p, bool rand=true, 
                       Nullable<NumericMatrix> twinsample_=R_NilValue){
  
  arma::uword m = M.n_rows;
  arma::uword q = M.n_cols;
  arma::vec filldist(m);
  //  unique output points
  arma::mat M_u = unique_rows(M);
  arma::uword no_u = M_u.n_rows;
  arma::vec filldist_u(no_u);
  
  if (no_u == 1){
    filldist = arma::ones(m) / m;
    return (filldist);
  }
  // NN will be used in approx_gen for mid points
  arma::uword no_nb = std::min(2*std::min(q, p), no_u-1);
  Rcpp::List dout = knn_my(M_u, no_nb) ;
  arma::umat dout_nn_index = dout["nn_index"];
  arma::mat dout_nn_dist = dout["nn_dist"];
  
  // approx_point
  arma::mat point_approx = approx_gen(M_u, dout, p, rand, twinsample_);
  
  // find the local filldist for each design output
  Rcpp::List d = knnx_my(M_u, point_approx, 1);   //needs special care as k=1, the output matrix format
  arma::umat d_nn_index = d["nn_index"];
  arma::mat d_nn_dist = d["nn_dist"];
  
  for(arma::uword i=0; i<no_u; i++){
    if (!arma::any(arma::any(d_nn_index==i))){
      filldist_u(i) = arma::mean(dout_nn_dist.row(i)) / 2;
    }else{
      filldist_u(i) = arma::max(d_nn_dist(arma::find(d_nn_index==i))); //needs special care; nn_dist need to be a col vec
    }
  }
  
  arma::rowvec tempt(q);
  for (arma::uword i=0; i<m; i++){
    tempt = M.row(i);
    filldist(i) = filldist_u(compare_v_m(tempt, M_u));
    if (filldist(i)==0){
      filldist(i) = 1e-30;
    }
  }
  return (filldist);
}

// Nearest neighbor with uq for EI ----------------------------------------------------------------
arma::rowvec nn_EI(arma::mat train, arma::mat test, arma::vec y){
  arma::uword m = train.n_rows;
  // arma::uword p = train.n_cols;
  arma::uword n = test.n_rows;
  // arma::uword K = std::min(2*p,m-1);
  arma::uword K = 1;
  arma::vec pred(n);
  arma::vec s(n);
  // float alpha = 1;
  
  // variance
  Rcpp::List din = knn_my(train, K); 
  arma::umat din_nn_index = din["nn_index"];
  arma::mat din_nn_dist = din["nn_dist"];
  
  arma::vec sigma2_p = arma::square(y - y(din_nn_index.col(0))) / din_nn_dist.col(0);
  arma::vec sigma2 = arma::ones(m) * arma::mean(sigma2_p);
  
  // knn between train and test
  Rcpp::List dx = knnx_my(train, test, 1); 
  arma::umat dx_nn_index = dx["nn_index"];
  arma::mat dx_nn_dist = dx["nn_dist"];
  
  pred = y(dx_nn_index);
  
  s = arma::sqrt(sigma2(dx_nn_index) % dx_nn_dist);
  double fmax = arma::max(y);
  arma::vec u = (pred-fmax) / s.clamp(1e-20, 1e20);
  arma::vec ei = s % (u % arma::normcdf(u) + arma::normpdf(u));
  
  arma::uword idx = ei.index_max();
  arma::rowvec u_new = test.row(idx);
  return (u_new);
}
// Criterion for IOSFD
arma::rowvec cri_iosfd(arma::mat D, arma::mat cand, 
                       Rcpp::List d_D, arma::vec filldist_out, 
                       arma::uword q, double lambda=0.5){
  arma::uword m = D.n_rows;
  arma::uword p = D.n_cols;
  
  // Compute nearest neighbors for CAND in D (replace with your NN library function)
  Rcpp::List d_D_cand = knnx_my(D, cand, 1);
  arma::umat d_D_cand_nn_index = d_D_cand["nn_index"];
  arma::mat d_D_cand_nn_dist = d_D_cand["nn_dist"];
  
  // Compute input fill distance
  arma::vec filldist_in(m);
  for (arma::uword i=0; i<m; i++) {
    if (!arma::any(arma::any(d_D_cand_nn_index==i))) {
      // Compute average distance of k-nearest neighbors
      filldist_in(i) = arma::mean(d_D_cand_nn_dist.row(i)) / 2.0;
    } else {
      filldist_in(i) = arma::max(d_D_cand_nn_dist(arma::find(d_D_cand_nn_index==i)));
    }
  }
  
  int intrinsic_dim = std::min(p, q);
  // double out_vol = arma::sum(arma::pow(filldist_out, intrinsic_dim)) * 
  //   std::pow(M_PI, intrinsic_dim/2.0) / std::tgamma(intrinsic_dim/2.0 + 1);
  // double alpha = arma::as_scalar(arma::pow(arma::mean(arma::pow(filldist_out, intrinsic_dim)), 1.0/intrinsic_dim)) /
  //   arma::as_scalar(arma::pow(arma::mean(arma::pow(filldist_in, p)), 1.0/p));
  
  double alpha =
    std::pow(arma::as_scalar(arma::mean(arma::pow(filldist_out, intrinsic_dim))), 1.0/intrinsic_dim) /
      std::pow(arma::as_scalar(arma::mean(arma::pow(filldist_in, p))), 1.0/p);
  
  
  arma::vec balance_cri = lambda * filldist_in + (1-lambda) * filldist_out / alpha;
  
  // Find the index with the maximum balance criterion
  arma::uword idx_perturb = balance_cri.index_max();
  arma::uvec sel2 = arma::find(d_D_cand_nn_index == idx_perturb);
  
  // Find the furthest point in the Voronoi cell of that point
  arma::rowvec new_point;
  if (!sel2.is_empty()) {
    // Select points in cand matching the criteria
    arma::mat u = cand.rows(sel2);
    arma::uword ind2 = d_D_cand_nn_dist.rows(sel2).index_max();
    new_point = u.row(ind2);
  } else {
    arma::mat query = D.row(idx_perturb);
    Rcpp::List tempt = knnx_my(cand, query, 1);
    arma::umat tempt_nn_index = tempt["nn_index"];
    new_point = cand.row(tempt_nn_index(0, 0));
  }
  
  return new_point;
}
// perturbation ----------------------------------------------------------------
//[[Rcpp::export]]
arma::rowvec perturb_cpp(arma::mat D, arma::vec filldist, int q, 
                         Nullable<NumericMatrix> CAND_=R_NilValue, 
                         std::string cri="EI", double balance_ratio=0.5,
                         bool rand=true, 
                         Nullable<NumericMatrix> twinsample_=R_NilValue){
  arma::uword m = D.n_rows;
  arma::uword p = D.n_cols;
  arma::uword K = std::min(2*p, m-1); //number of neighbors between design D
  arma::uword ind = filldist.index_max();
  bool din_flag = false; 
  
  arma::mat CAND;
  Rcpp::List din;
  // Candidate generation
  if (CAND_.isNotNull()){
    NumericMatrix CAND_tempt(CAND_);  // casting to arma matrix
    CAND = as<arma::mat>(wrap(CAND_tempt));
  }else{
    // find the nearest neighbor around the point with largest filldist
    din = knn_my(D, K); // knn among design points. need modifying
    arma::umat din_nn_index = din["nn_index"];
    arma::mat din_nn_dist = din["nn_dist"];
    din_flag = true;
    arma::uvec sel(K+1);
    sel(0) = ind;
    sel.rows(1, K) = din_nn_index.row(ind).t();
    arma::mat cen = D.rows(sel);
    arma::vec rad(K+1);
    rad(0) = din_nn_dist.row(ind)(0);
    rad.rows(1, K) = din_nn_dist.row(ind).t();
    
    if (cri=="greedy"){
      // greedy
      // uniform points around the largest filldist input
      double r = rad(K);
      // hypercube points
      arma::rowvec l(p);
      arma::rowvec h(p);
      l = clamp((cen.row(0)-r), 0, 1);
      h = clamp((cen.row(0)+r), 0, 1);
      arma::mat cube = arma::randu(10*p*(K+1), p);
      cube.each_row() %= (h-l);
      cube.each_row() += l;
      // ball points
      arma::mat Ball_point(p*10*(K+1), p); 
      for (int i=0; i<(K+1); i++){
        Ball_point.rows(p*10*i, p*10*(i+1)-1) = ball_gen(cen.row(i), rad(i), 
                        p*10, 0, R_NilValue, rand, twinsample_);
      }
      CAND = join_cols(cube,Ball_point);
    }else{
      // hypercube points
      arma::mat cube = arma::randu(10*m, p);
      // ball points
      arma::mat Ball_point(p*10*(K+1), p);
      for (int i=0; i<(K+1); i++){
        Ball_point.rows(p*10*i, p*10*(i+1)-1) = ball_gen(cen.row(i), rad(i), 
                        p*10, 0, R_NilValue, rand, twinsample_);
      }
      // mid points
      arma::mat mid_point(K*m, p);
      for (int i=0; i<K; i++){
        mid_point.rows(m*i, m*(i+1)-1) = (D+D.rows(din_nn_index.col(i)))/2;
      }
      CAND = join_cols(cube, Ball_point, mid_point);
    }
    CAND = clamp(CAND, 0, 1);
    CAND = unique_rows(CAND);
  }
  
  // Choose the next design point
  arma::rowvec u_new;
  if(cri=="EI"){
    u_new = nn_EI(D, CAND, filldist);
  }else if(cri=="iosfd"){
    if(!din_flag){
      din = knn_my(D, K);
    }
    u_new = cri_iosfd(D, CAND, din, filldist, q, balance_ratio);
  }else{
    Rcpp::List d = knnx_my(D, CAND, 1); 
    arma::umat d_nn_index = d["nn_index"];
    arma::mat d_nn_dist = d["nn_dist"];
    arma::uvec sel2 = arma::find(d_nn_index==ind);
    if (sel2.empty()){
      arma::mat query = D.row(ind);
      Rcpp::List tempt = knnx_my(CAND, query, 1);
      arma::umat tempt_nn_index = tempt["nn_index"];
      u_new = CAND.row(tempt_nn_index(0, 0));
    }else{
      arma::mat u = CAND.rows(sel2);
      arma::uword ind2 = d_nn_dist.rows(sel2).index_max();
      u_new = u.row(ind2);
    }
    
  }
  return (u_new);
}
