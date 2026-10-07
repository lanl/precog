# depends: 
create_embed_matrix               <- function(synthetic, h, k = 4){
  s_idx                           <- 1
  for (s in synthetic){
    s$ts_id                       <- s_idx
    synthetic[[s_idx]]            <- s
    s_idx                         <- s_idx + 1
  }
  embed_mat                       <- lapply(synthetic,function(x){ embed( pmax(1e-8,x$ts),k+h)})
  embed_mat                       <- do.call(rbind,embed_mat)
  embed_mat                       <- embed_mat[,ncol(embed_mat):1] #casey
  RowVar                          <- function(x, ...) {
    rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
  }
  rows_to_delete = RowVar(embed_mat[,1:k])
  embed_mat                       <- embed_mat[which(rows_to_delete > 0),]
  embed_mat_X                     <- embed_mat[,1:k]
  embed_mat_y                     <- embed_mat[,(k+1):(k+h)]
  ret_list                        <- list()
  ret_list[[1]]                   <- embed_mat_X
  ret_list[[2]]                   <- embed_mat_y
  return (ret_list)
}