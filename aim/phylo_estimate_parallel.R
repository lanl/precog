# parallelize parameter estimation from trees, using TreePar
# Julie Spencer
# Scottie Alexander
# 09/12/2024

# load libraries
library(parallel)
library(doParallel)
library(MASS)

### run tree estimates ###
# ============================================================================ #
# function to estimate maximum likelihood #
process_tree <- function(tree) {
  start_lk <- LikTypesSTT(c(3,3,3,3,3), phylo=tree,
                      fix=rbind(c(6,7,8),c(-5,0,0),c(1,1,1)),
                      sampfrac=s,survival=0,posR=0)

  if (is.finite(start_lk)) {
    p_H0 <- optim(c(3,3,3,3,3), LikTypesSTT, 
                  phylo=tree, fix=rbind(c(6,7,8), c(-5,0,0), c(1,1,1)),
                  sampfrac=c(1,1), survival=0, posR=0, control=list(maxit=10000))
    out <- p_H0$par
  } else {
    out <- rep(NaN, 5)
  }
  
  return(out)
}

# ============================================================================ #
###all subtrees###
all_subtrees = list(subtrees_H0, subtrees_H1, subtrees_H2, subtrees_H3)

odir <- file.path(getwd(), "results_subtrees")

for (k in 1:4) {
  
  trees <- all_subtrees[[k]]
  
  tmp2 <- mclapply(trees, process_tree, mc.cores = 8)
  
  df <- as.data.frame(do.call(rbind, tmp2))
  colnames(df) <- c("beta11", "beta12", "beta21", "beta22", "gamma")
  
  write.csv(df, file.path(odir, paste0("result_18Dec24_H", k-1, ".csv")), row.names=FALSE)
  
}


###all trees###
all_trees = list(trees_H0, trees_H1, trees_H2, trees_H3)

odir <- file.path(getwd(), "results_trees")

for (k in 1:1) {
  
  dur <- system.time({
  trees <- trees_H0[1:12]#all_trees[[k]]
  
  tmp2 <- mclapply(trees, process_tree, mc.cores = 6)

  df <- as.data.frame(do.call(rbind, tmp2))
  colnames(df) <- c("beta11", "beta12", "beta21", "beta22", "gamma")
  })
  
  #write.csv(df, file.path(odir, paste0("result_17Dec24_H", k-1, ".csv")), row.names=FALSE)
  
}

# ============================================================================ #
