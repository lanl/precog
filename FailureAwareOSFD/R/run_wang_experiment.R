## Replicate Wang et al. (2024) first experiment: inverse-radius (Figure 7)
## - p = 2, q = 2
## - epsilon = 0.1
## - initial design size n0 = 10 (random LHD)
## - 20 replications
## - report output-space fill distance vs n
##
## Paper refs: inverse-radius definition (eq 4) and Figure 7 setup. 
## See: :contentReference[oaicite:5]{index=5} and :contentReference[oaicite:6]{index=6}

suppressPackageStartupMessages({
  library(lhs)
  library(RANN)
  library(ggplot2)
  library(patchwork)
  library(parallel)
  library(doParallel)
})

## ---- 1) Test function (inverse-radius) ----
f_ir <- function(x, eps = 0.1) {
  # x: numeric vector length 2 in [0,1]^2
  stopifnot(length(x) == 2)
  x1 <- x[1]; x2 <- x[2]
  y1 <- 1 / sqrt(x1^2 + x2^2 + 2*eps)
  y2 <- atan2(x2, x1)  # stable arctan(x2/x1) on quadrant; in [0, pi/2] here
  c(y1, y2)
}

## Vectorized mapping for matrices (n x 2 -> n x 2)
f_ir_mat <- function(X, eps = 0.1) {
  X <- as.matrix(X)
  x1 <- X[,1]; x2 <- X[,2]
  cbind(
    1 / sqrt(x1^2 + x2^2 + 2*eps),
    atan2(x2, x1)
  )
}

## ---- 2) Output-space fill distance estimator ----
## phi_mM(M) = max_{y in Y} min_{yi in M} ||y - yi||.
## Approximate by large MC sample A ~ f(U([0,1]^2)).
fill_distance_mc <- function(M, eps = 0.1, Nmc = 200000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  M <- as.matrix(M)
  stopifnot(ncol(M) == 2)

  Xmc <- matrix(runif(Nmc * 2), ncol = 2)
  Amc <- f_ir_mat(Xmc, eps = eps)

  # nearest neighbor distance from each a in A to design outputs M
  nn <- RANN::nn2(data = M, query = Amc, k = 1)
  max(nn$nn.dists[,1])
}

## ---- 3) Baseline: random LHD in input, evaluated in output ----
isfd_random_lhd_curve <- function(n_grid, eps = 0.1, Nmc = 200000) {
  vapply(n_grid, function(n) {
    X <- lhs::randomLHS(n, 2)
    M <- f_ir_mat(X, eps = eps)
    fill_distance_mc(M, eps = eps, Nmc = Nmc)
  }, numeric(1))
}

## ---- 4) Your OSFD-EI wrapper (using your constrained_osfd_ei_with_retry) ----
## This assumes constrained_osfd_ei_with_retry() returns a design in order:
## res$D_success (n x p) and res$Y_success (n x q)
osfd_ei_curve_via_your_code <- function(
  n0 = 10,
  n_max = 150,
  eps = 0.1,
  Ncand = 200000,      # candidate pool size for EI search
  cand_batch = 2000,   # how many candidates to score per iteration inside your routine
  Nmc = 200000,        # MC size for fill-distance evaluation
  seed = 1,
  mc.cores = 1,
  verbose = TRUE
) {
  set.seed(seed)

  # f <- function(x) f_ir(x, eps = eps)

  # finite candidate pool, as your implementation prefers
  CAND <- lhs::randomLHS(Ncand, 2)

  # Sequential (one point at a time) to match the paper’s sequential design spirit

  res <- constrained_osfd_ei_with_retry(
    f = f_ir,
    p = 2,
    q = 2,
    CAND = CAND,
    n_success = n_max,
    n_ini_success = n0,
    batch_size = 50,          # add one per iteration
    cand_batch = cand_batch,
    beta = 2.0,              # if your EI/UCB logic uses it; keep fixed for sanity check
    p_floor = 0.0,
    tau_hard = NULL,         # no constraints in this benchmark
    repel_radius = 0.1,      # paper doesn’t impose repulsion in input; set 0 for fidelity
    update_feas_every = 1,
    mc.cores = mc.cores,
    verbose = verbose
  )

  D <- as.matrix(res$D_success)
  Y <- as.matrix(res$Y_success)

  if (nrow(D) < n_max) {
    warning(sprintf("OSFD returned only %d points (expected %d). Using what we have.",
                    nrow(D), n_max))
    n_max <- nrow(D)
  }

  n_grid <- seq.int(n0, n_max)

  # fill distance after each prefix of the sequential design
  fd <- vapply(n_grid, function(n) {
    fill_distance_mc(Y[1:n, , drop = FALSE], eps = eps, Nmc = Nmc)
  }, numeric(1))

  list(n = n_grid, fill = fd, D = D, Y = Y)
}

## ---- 5) Run the experiment (20 replications) ----
run_experiment <- function(
  R = 20,
  n0 = 10,
  n_max = 150,
  eps = 0.1,
  Ncand = 200000,
  cand_batch = 2000,
  Nmc = 200000,
  mc.cores = 1
) {
  n_grid <- seq.int(n0, n_max)

  # OSFD-EI replications
  osfd_mat <- matrix(NA_real_, nrow = R, ncol = length(n_grid))
  for (r in seq_len(R)) {
    cat(sprintf("OSFD-EI replication %d/%d\n", r, R))
    out <- osfd_ei_curve_via_your_code(
      n0 = n0, n_max = n_max, eps = eps,
      Ncand = Ncand, cand_batch = cand_batch,
      Nmc = Nmc, seed = 1000 + r,
      mc.cores = mc.cores, verbose = TRUE
    )
    # align in case your routine returned fewer points
    idx <- match(out$n, n_grid)
    osfd_mat[r, idx] <- out$fill
  }

  # ISFD (random LHD) replications
  isfd_mat <- matrix(NA_real_, nrow = R, ncol = length(n_grid))
  for (r in seq_len(R)) {
    cat(sprintf("ISFD (random LHD) replication %d/%d\n", r, R))
    set.seed(2000 + r)
    isfd_mat[r, ] <- isfd_random_lhd_curve(n_grid, eps = eps, Nmc = Nmc)
  }

  list(n = n_grid, osfd = osfd_mat, isfd = isfd_mat)
}

summarize_curves <- function(n_grid, mat, method_name) {
  data.frame(
    n = n_grid,
    method = method_name,
    mean = apply(mat, 2, mean, na.rm = TRUE),
    q05  = apply(mat, 2, quantile, probs = 0.05, na.rm = TRUE),
    q95  = apply(mat, 2, quantile, probs = 0.95, na.rm = TRUE)
  )
}

plot_summary <- function(df) {
  ggplot(df, aes(x = n, y = mean)) +
    geom_ribbon(aes(ymin = q05, ymax = q95), alpha = 0.2) +
    geom_line(linewidth = 1) +
    facet_wrap(~method, scales = "free_y") +
    labs(
      x = "Run size (n)",
      y = "Estimated output-space fill distance",
      title = "Inverse-radius experiment: fill distance vs run size",
      subtitle = "Bands = 5th–95th quantiles over replications"
    ) +
    theme_bw()
}

## ----------------------------
filldist_mc <- function(Y_design, eps = 0.1, Nmc = 200000, seed = 999) {
  set.seed(seed)
  Xmc <- matrix(runif(Nmc * 2), ncol = 2)
  Ymc <- f_ir_mat(Xmc, eps = eps)

  # nearest neighbor distance from each Ymc point to the design outputs
  # (simple brute force would be slow; use RANN if available)
  if (requireNamespace("RANN", quietly = TRUE)) {
    nn <- RANN::nn2(data = Y_design, query = Ymc, k = 1)
    return(max(nn$nn.dists[,1]))
  } else {
    # fallback (slow): compute in chunks
    maxd <- 0
    chunk <- 5000
    for (i in seq(1, nrow(Ymc), by = chunk)) {
      j <- min(i + chunk - 1, nrow(Ymc))
      Q <- Ymc[i:j, , drop = FALSE]
      dmin <- apply(Q, 1, function(q) min(sqrt(rowSums((Y_design - q)^2))))
      maxd <- max(maxd, max(dmin))
    }
    return(maxd)
  }
}

output_boundary <- function(eps = 0.1, ngrid = 600) {
  th <- seq(0, pi/2, length.out = ngrid)
  rmax1 <- 1 / pmax(cos(th), 1e-12)
  rmax2 <- 1 / pmax(sin(th), 1e-12)
  rmax  <- pmin(rmax1, rmax2)

  y1_min <- 1 / sqrt(rmax^2 + 2*eps)             # lower boundary
  y1_max <- rep(1 / sqrt(0 + 2*eps), length(th)) # upper boundary

  data.frame(
    th = c(th, rev(th)),
    y1 = c(y1_min, rev(y1_max))
  )
}

# ## ---- 6) Actually run (edit sizes for speed vs fidelity) ----
# ## Notes:
# ## - Nmc controls accuracy of fill-distance evaluation (bigger = closer to paper, slower).
# ## - Ncand/cand_batch control how aggressively your EI searches.
# ## - Start smaller, then crank up if curves look noisy.
# usethis::use_description()
# usethis::use_namespace()
# pkgload::load_all()
# exp1 <- run_experiment(
#   R = 20,
#   n0 = 10,
#   n_max = 150,
#   eps = 0.1,
#   Ncand = 200000,
#   cand_batch = 2000,
#   Nmc = 200000,
#   mc.cores = 50
# )

# df_osfd <- summarize_curves(exp1$n, exp1$osfd, "OSFD-EI (your code)")
# df_isfd <- summarize_curves(exp1$n, exp1$isfd, "ISFD (random LHD)")
# df_all <- rbind(df_osfd, df_isfd)

# print(plot_summary(df_all))

# ## If you want a single overlay plot instead of facets:
# ggplot(df_all, aes(x = n, y = mean, linetype = method)) +
#   geom_ribbon(aes(ymin = q05, ymax = q95, fill = method), alpha = 0.15) +
#   geom_line(linewidth = 1) +
#   labs(
#     x = "Run size (n)",
#     y = "Estimated output-space fill distance",
#     title = "Inverse-radius experiment (Figure 7 style): OSFD-EI vs random LHD"
#   ) +
#   theme_bw() +
#   guides(fill = "none")



## ----------------------------
## 4) Run ISFD and OSFD-EI (your function) for Figure-1-like settings
## ----------------------------
usethis::use_description()
usethis::use_namespace()
pkgload::load_all()
eps <- 0.1
n   <- 50
n0  <- 5      # Figure 1 caption says init size 5 :contentReference[oaicite:5]{index=5}

set.seed(123)

# ISFD: input space-filling, size n
X_isfd <- lhs::maximinLHS(n, 2)
Y_isfd <- f_ir_mat(X_isfd, eps = eps)

# OSFD-EI using YOUR routine
# Important: use batch_size=1 here to keep it truly sequential like the paper’s OSFD setup. :contentReference[oaicite:6]{index=6}
# Candidate pool: big LHS menu in [0,1]^2 (scaled already)
Ncand <- 300000
CAND  <- lhs::randomLHS(Ncand, 2)

# Assumes you already have constrained_osfd_ei_with_retry() loaded in your session
res_osfd <- constrained_osfd_ei_with_retry(
  f = f_ir,
  p = 2,
  q = 2,
  CAND = CAND,
  n_success = n,           # target number of successful points
  n_ini_success = n0,   # init successes
  batch_size = 1,       # sequential
  cand_batch = 5000,    # score this many candidates per iteration (increase if needed)
  beta = 2.0,           # doesn’t matter much here since p_feas ~ 1
  p_floor = 0.0,
  tau_hard = NULL,
  repel_radius = 0.0,   # not needed for batch_size=1
  update_feas_every = 1,
  mc.cores = 1,
  verbose = TRUE
)

X_osfd <- as.matrix(res_osfd$D_success)
Y_osfd <- as.matrix(res_osfd$Y_success)

# Numeric sanity checks (smaller output fill distance is better)
cat("\nApprox output fill distance (Monte Carlo):\n")
cat(sprintf("  ISFD:    %.4f\n", filldist_mc(Y_isfd, eps = eps)))
cat(sprintf("  OSFD-EI: %.4f\n\n", filldist_mc(Y_osfd, eps = eps)))

## ----------------------------
## 5) Plot in a Figure-1-like layout (inputs top row, outputs bottom row)
## ----------------------------
boundary <- output_boundary(eps = eps)

df_in <- rbind(
  transform(data.frame(x1 = X_isfd[,1], x2 = X_isfd[,2]), method = "ISFD (input LHS)"),
  transform(data.frame(x1 = X_osfd[,1], x2 = X_osfd[,2]), method = "OSFD-EI (your code)")
)

df_out <- rbind(
  transform(data.frame(y1 = Y_isfd[,1], th = Y_isfd[,2]), method = "ISFD (input LHS)"),
  transform(data.frame(y1 = Y_osfd[,1], th = Y_osfd[,2]), method = "OSFD-EI (your code)")
)

p_in <- ggplot(df_in, aes(x1, x2)) +
  geom_point(size = 2, alpha = 0.9) +
  # coord_equal() +
  facet_wrap(~method, nrow = 1) +
  labs(title = "Input space designs", x = expression(x[1]), y = expression(x[2])) +
  theme_bw()

p_out <- ggplot(df_out, aes(y1, th)) +
  geom_point(size = 2, alpha = 0.9) +
  # coord_equal()+
  geom_path(data = boundary, aes(y1, th), inherit.aes = FALSE,
            linetype = "dashed", linewidth = 0.6) +
  facet_wrap(~method, nrow = 1) +
  labs(title = "Output space images", x = expression(y[1]), y = expression(theta)) +
  theme_bw()

print(p_in)
print(p_out)

p_in / p_out
