library(parallel)
library(foreach)
library(doParallel)
library(GGally)
library(data.table)

# Load helper functions
source(here::here("R", "source_helpers.R"))
source_helpers()

setwd(here::here())
# usethis::use_description()
# usethis::use_namespace()
# pkgload::load_all(this.path::here())

p <- 8
q <- 10
nsuccesses = 175000
pkg_location = this.path::here()
mc.cores = 99

# Example: user-supplied f that sometimes fails
f <- function(x) {
  generate_smoa_synthetic(x, "sir_rollercoaster", "OSFD", pkg_location)
}

print("Clearing previously built synthetic data and visualizations.")
unlink(paste0(pkg_location, "/data/OSFD", '/', '*'))
unlink(paste0(pkg_location, "/data/ISFD", '/', '*'))
unlink(paste0(pkg_location, "/viz", '/', '*'))

# Finite candidate pool (strongly recommended for blacklist logic)
set.seed(13)
CAND <- lhs::randomLHS(200000, 8)

print("About to enter OSFD function.")
res <- constrained_osfd_ei_with_retry(
  f, p, q,
  CAND,
  nsuccesses,
  n_ini_success = floor(3 * nsuccesses / 4),
  batch_size = mc.cores,
  n_replicates = 100,
  cand_batch = 50000,     # score this many candidates per iteration
  beta = 2.0,
  p_floor = 0.05,
  tau_hard = NULL,
  repel_radius = 0.05,    # input-space diversity scale
  update_feas_every = 1,
  mc.cores = mc.cores,
  verbose = TRUE, 
  pkg_location = pkg_location
)
print("We have now exited the OSFD function.")

D <- res$D_success
Y <- res$Y_success


# ISFD: input space-filling, size n
X_isfd <- lhs::randomLHS(nsuccesses, p) 

f <- function(x) {
  generate_smoa_synthetic(x, "sir_rollercoaster", "ISFD", pkg_location)
}

print("About to sample the ISFD.")
cl <- parallel::makeCluster(mc.cores)
doParallel::registerDoParallel(cl)

# send the variable to workers
parallel::clusterExport(cl, varlist = c("pkg_location", "f", "q"), envir = environment())

parallel::clusterEvalQ(cl, {
  # suppressPackageStartupMessages(library(pkgload))
  source(here::here("R", "source_helpers.R")); source_helpers()
  # pkgload::load_all(pkg_location, quiet = TRUE)
  NULL
})
results <- foreach::foreach(
  i = seq_len(nrow(X_isfd)),
  .combine = 'rbind',
  .multicombine = TRUE,
  .inorder = TRUE
) %dopar% {
  # source(save_eval_location)
  # source(save_generate_location)
  # pkgload::load_all(pkg_location, quiet = TRUE)
  x <- X_isfd[i, , drop = FALSE]
  matrix(safe_eval_f(f, x, q)$y, nrow = 1)
}
parallel::stopCluster(cl) 

Y_isfd <- results[complete.cases(results),]

save(Y_isfd, file=here::here("data", "Y_isfd.RData"))
save(Y, file=here::here("data", "Y_osfd.RData"))
load(file=here::here("data", "Y_isfd.RData"))
load(file=here::here("data", "Y_osfd.RData"))

dfY <- as.data.frame(Y_isfd)
p1 = GGally::ggpairs(dfY)
p <- GGally::ggpairs(dfY) +
  ggtitle("Pairs plot of Y ISFD") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(
  filename = here::here("viz", "Y_isfd_ggpairs.png"),
  plot = p,
  width = 12, height = 12, units = "in", dpi = 300
)

dfY <- as.data.frame(Y)
p <- GGally::ggpairs(dfY) +
  ggtitle("Pairs plot of Y OSFD") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(
  filename = here::here("viz", "Y_osfd_ggpairs.png"),
  plot = p,
  width = 12, height = 12, units = "in", dpi = 300
)

###################### UMAP
## Load Disease
DISEASE_KEY = "Ebola virus disease"
FILE_NAME = "Ebola_ginkgo.RDS"
### demonstrating using Global_Covid
tmp = dplyr::as_tibble(readRDS(paste0(here::here("data", "processed_features_data"),"/", FILE_NAME)))
if(nrow(tmp) > 10000){
  tmp <- tmp[sample(1:nrow(tmp), 10000, replace = F),]
}
tmp = tmp[!is.na(tmp$entropy),]
train_disease_data = tmp %>%
  subset(h == 1) %>%
  dplyr::select(gr12_div_23, last_div_max, coefvar, 
                gam_with_div_without, avg_recent_div_avg_global,
                diff_zscore, entropy, relative_increases,
                prop_since_peak, seasonality)%>%
  dplyr::distinct()
# train_disease_data = screen_multivariate_outliers(train_disease_data, threshold = outlier_threshold)

## combine with synthetic
colnames(Y_isfd) = colnames(train_disease_data)
colnames(Y) = colnames(train_disease_data)
all_other_data = rbind(train_disease_data, Y_isfd, Y)

global_coviddf <- as.matrix(rbind(all_other_data, train_disease_data))

## fit UMAP (2-dimensions)
umap_df = create_umap_df(DISEASE_KEY, FILE_NAME, Y_isfd, Y)

final_plot = create_umap_plot(umap_df, "Ebola", data_color = "olivedrab")

ggsave(
  filename = here::here("viz", "umap.png"),
  plot = final_plot,
  width = 20, height = 14, units = "in", dpi = 300
) 
