################################
### New SIR Experiment Plots ###
################################
library(ggplot2)
library(ggplot2)
library(patchwork)
library(tidyr)
library(latex2exp)
setwd(here::here())
source(here::here('R', 'sir_experiment_setup.R'))

# alpha_bounds = c(0.001, 100)
# reproduction_number_bounds = c(1.001, 100)
# pit_bounds = c(0, 50)
# piv_bounds = c(0.005, 0.75)
# s0_bounds = c(0.95, 0.999)
n = 100000
alpha_lo <- alpha_bounds[1];  alpha_hi <- alpha_bounds[2]
rho_lo   <- reproduction_number_bounds[1];  rho_hi   <- reproduction_number_bounds[2]
s0_lo    <- s0_bounds[1];  s0_hi    <- s0_bounds[2]


load(file = here::here('data', paste0('fullOutput_OSFD_basic_n',format(n, scientific = FALSE),'.RData')))
X_OSFD = data.frame(res$D_success)
colnames(X_OSFD) = c('V1','V2','V3')
X_OSFD$alpha <- alpha_lo + X_OSFD$V1 * (alpha_hi - alpha_lo)
X_OSFD$s0    <- s0_lo    + X_OSFD$V3 * (s0_hi    - s0_lo)
X_OSFD$rho   <- rho_lo   + X_OSFD$V2 * (rho_hi   - rho_lo)
X_OSFD$beta <- X_OSFD$alpha / X_OSFD$rho
X_OSFD$logtrans = log(X_OSFD$alpha/X_OSFD$beta)

y_OSFD = res$Y_success

X_all = data.frame(res$X_all)
colnames(X_all) = c('V1','V2','V3')
X_all$alpha <- alpha_lo + X_all$V1 * (alpha_hi - alpha_lo)
X_all$s0    <- s0_lo    + X_all$V3 * (s0_hi    - s0_lo)
X_all$rho   <- rho_lo   + X_all$V2 * (rho_hi   - rho_lo)
X_all$beta <- X_all$alpha / X_all$rho
X_all$logtrans = log(X_all$alpha/X_all$beta)


load(file = here::here('data', paste0('fullOutput_OSFD_wSIRStuff_n',format(n, scientific = FALSE),'.RData')))
X_OSFDsir = data.frame(res$D_success)
colnames(X_OSFDsir) = c('V1','V2','V3')
X_OSFDsir$alpha <- alpha_lo + X_OSFDsir$V1 * (alpha_hi - alpha_lo)
X_OSFDsir$s0    <- s0_lo    + X_OSFDsir$V3 * (s0_hi    - s0_lo)
X_OSFDsir$rho   <- rho_lo   + X_OSFDsir$V2 * (rho_hi   - rho_lo)
X_OSFDsir$beta <- X_OSFDsir$alpha / X_OSFDsir$rho
X_OSFDsir$logtrans = log(X_OSFDsir$alpha/X_OSFDsir$beta)

y_OSFDsir = res$Y_success


load(file = here::here('data', paste0('inputs_and_outputs_sirBasicLHS_n',format(n, scientific = FALSE),'.RData')))
results_ISFD = results
results_ISFD$logtrans = log(results_ISFD$alpha/results_ISFD$beta)
results_ISFD$rho = results_ISFD$alpha/results_ISFD$beta


load(file = here::here('data', paste0('inputs_and_outputs_sirMappings_n',format(n, scientific = FALSE),'.RData')))
results_sirMappings = results
results_sirMappings$logtrans = log(results_sirMappings$alpha/results_sirMappings$beta)
results_sirMappings$rho = results_sirMappings$alpha/results_sirMappings$beta

#################################
### Plot OSFD Failure Surface ###
#################################


# TO_PLOT = X_all
# TO_PLOT$pred = predict(res$feas_model, newdata = X_all, type = 'response')

# ggplot(TO_PLOT)+
#   geom_point(aes(x=beta, y=logtrans, col = pred))+
#   theme_classic()+
#   scale_x_continuous(expand=c(0,0))+
#   scale_y_continuous(expand=c(0,0))+
#   xlab('beta')+
#   ylab('log(alpha/beta)')


# Combine both LHS and SIR inverse maps for validation set
TO_TEST = rbind(results_ISFD, results_sirMappings)
cat(sprintf("Validation set: %d total samples (%d from LHS + %d from SIR Maps)\n",
            nrow(TO_TEST), nrow(results_ISFD), nrow(results_sirMappings)))
NSAMPLES = seq(1000, 5000, 1000)  # Only use n where data files exist
logit = function(x){log(x/(1-x))}
expit = function(x){exp(x)/(1+exp(x))}
TO_TEST$y = logit(TO_TEST$piv)

### OSFD
# TO_TEST = rbind(data.frame(piv = y_OSFD[,1], pit = y_OSFD[,2],
#                            alpha = X_OSFD$alpha, s0 = X_OSFD$s0, rho = X_OSFD$rho, beta = X_OSFD$beta)[15000:nrow(y_OSFD),],
#                 results_ISFD[15000:nrow(results_ISFD),c('piv','pit','alpha','s0','rho','beta')])
ggplot(TO_TEST)+
  geom_point(aes(x=alpha, y=rho, col = piv))+
  theme_classic()+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  xlab('alpha')+
  ylab('rho')


plot(TO_TEST$s0, TO_TEST$piv)


# Set up parallel cluster
library(foreach)
library(doSNOW)
library(dplyr)
library(tidyr)

# Create grid of all n x replicate combinations
n_rep_grid <- expand.grid(
  n = NSAMPLES,
  rep = 1:20,
  stringsAsFactors = FALSE
)

cat(sprintf("Processing %d n x replicate combinations in parallel...\n", nrow(n_rep_grid)))

n_cores <- min(parallel::detectCores() - 2, 20)
cl <- parallel::makeCluster(n_cores)
doSNOW::registerDoSNOW(cl)

cat(sprintf("Running emulation experiment in parallel on %d cores\n", n_cores))

# Export necessary variables to workers
parallel::clusterExport(cl, varlist = c("n_rep_grid", "TO_TEST", "alpha_lo", "alpha_hi",
                                         "rho_lo", "rho_hi", "s0_lo", "s0_hi",
                                         "logit", "expit"),
                        envir = environment())

# Set up progress bar
pb <- txtProgressBar(max = nrow(n_rep_grid), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# Run parallel loop - one iteration per n-replicate combination
results_list <- foreach::foreach(
  i = 1:nrow(n_rep_grid),
  .packages = c("mgcv", "here"),
  .combine = rbind,
  .inorder = FALSE,
  .errorhandling = "pass",
  .options.snow = opts
) %dopar% {

  # Wrap in tryCatch for error handling
  tryCatch({
    n <- n_rep_grid$n[i]
    rep_num <- n_rep_grid$rep[i]

    mae_row <- data.frame(n = n, replicate = rep_num, isfd = NA, osfd = NA, wangorig = NA, error_message = NA_character_)

    # Check if all three methods have files for this n and replicate
    budget_file <- Sys.glob(here::here('data', sprintf('fullOutput_OSFD_basic_budget%d_n*_rep%d.RData', n, rep_num)))
    isfd_file <- Sys.glob(here::here('data', sprintf('inputs_and_outputs_sirBasicLHS_n%s_rep%d.RData', format(n, scientific = FALSE), rep_num)))
    wangorig_file <- Sys.glob(here::here('data', sprintf('fullOutput_OSFD_WangOrig_budget%d_n*_rep%d.RData', n, rep_num)))

    # Only process if all three methods have data
    if (length(budget_file) == 0 || length(isfd_file) == 0 || length(wangorig_file) == 0) {
      # Skip this replicate if any method is missing
      return(NULL)
    }

    ### OSFD (using budget-based files)
    res <- NULL
    suppressWarnings({
      load(file = budget_file[1])
    })
    X_OSFD <- data.frame(res$D_success)
    colnames(X_OSFD) <- c('V1','V2','V3')
    X_OSFD$alpha <- alpha_lo + X_OSFD$V1 * (alpha_hi - alpha_lo)
    X_OSFD$s0    <- s0_lo    + X_OSFD$V3 * (s0_hi    - s0_lo)
    X_OSFD$rho   <- rho_lo   + X_OSFD$V2 * (rho_hi   - rho_lo)
    X_OSFD$beta <- X_OSFD$alpha / X_OSFD$rho
    X_OSFD$logtrans <- log(X_OSFD$alpha/X_OSFD$beta)
    y_OSFD <- res$Y_success

    # Use all available data (not just first n)
    actual_n <- nrow(y_OSFD)
    TO_FIT <- data.frame(piv = y_OSFD[,1], pit = y_OSFD[,2],
                         alpha = X_OSFD$alpha, s0 = X_OSFD$s0, rho = X_OSFD$rho)
    TO_FIT$y <- logit(TO_FIT$piv)
    fit <- mgcv::bam(y~s(alpha, rho, bs = 'tp') + s0, data = TO_FIT)
    preds <- expit(predict(fit, newdata = TO_TEST))
    mae_row$osfd <- mean(abs(TO_TEST$piv - preds))

    ### ISFD (Basic LHS)
    results <- NULL
    load(file = isfd_file[1])
    results_ISFD <- results
    results_ISFD$logtrans <- log(results_ISFD$alpha/results_ISFD$beta)
    results_ISFD$rho <- results_ISFD$alpha/results_ISFD$beta
    TO_FIT <- results_ISFD
    TO_FIT$y <- logit(TO_FIT$piv)
    fit <- mgcv::bam(y~s(alpha, rho, bs = 'tp') + s0, data = TO_FIT)
    preds <- expit(predict(fit, newdata = TO_TEST))
    mae_row$isfd <- mean(abs(TO_TEST$piv - preds))

    ### Wang-like OSFD
    res <- NULL
    suppressWarnings({
      load(file = wangorig_file[1])
    })
    X_WangOrig <- data.frame(res$D_success)
    colnames(X_WangOrig) <- c('V1','V2','V3')
    X_WangOrig$alpha <- alpha_lo + X_WangOrig$V1 * (alpha_hi - alpha_lo)
    X_WangOrig$s0    <- s0_lo    + X_WangOrig$V3 * (s0_hi    - s0_lo)
    X_WangOrig$rho   <- rho_lo   + X_WangOrig$V2 * (rho_hi   - rho_lo)
    X_WangOrig$beta <- X_WangOrig$alpha / X_WangOrig$rho
    X_WangOrig$logtrans <- log(X_WangOrig$alpha/X_WangOrig$beta)
    y_WangOrig <- res$Y_success

    # Use all available data (not just first n)
    actual_n_wang <- nrow(y_WangOrig)
    TO_FIT_WANG <- data.frame(piv = y_WangOrig[,1], pit = y_WangOrig[,2],
                         alpha = X_WangOrig$alpha, s0 = X_WangOrig$s0, rho = X_WangOrig$rho)
    TO_FIT_WANG$y <- logit(TO_FIT_WANG$piv)
    fit_wang <- mgcv::bam(y~s(alpha, rho, bs = 'tp') + s0, data = TO_FIT_WANG)
    preds_wang <- expit(predict(fit_wang, newdata = TO_TEST))
    mae_row$wangorig <- mean(abs(TO_TEST$piv - preds_wang))

    # Return the row for this iteration
    mae_row

  }, error = function(e) {
    # On error, return error info
    return(data.frame(
      n = n,
      replicate = rep_num,
      isfd = NA,
      osfd = NA,
      wangorig = NA,
      error_message = as.character(e$message),
      stringsAsFactors = FALSE
    ))
  })
}

# Close progress bar and stop the cluster
close(pb)
parallel::stopCluster(cl)

# Clean up results - use dplyr::bind_rows for proper handling
MAE <- dplyr::bind_rows(results_list)

# Ensure proper column types
MAE$n <- as.numeric(MAE$n)
MAE$replicate <- as.numeric(MAE$replicate)
MAE$isfd <- as.numeric(MAE$isfd)
MAE$osfd <- as.numeric(MAE$osfd)
MAE$wangorig <- as.numeric(MAE$wangorig)
MAE$error_message <- as.character(MAE$error_message)

# Separate errors from valid data
errors <- MAE[!is.na(MAE$error_message), , drop = FALSE]
MAE <- MAE[is.na(MAE$error_message), , drop = FALSE]

cat(sprintf("Total rows: %d replicates across %d sample sizes\n", nrow(MAE), n_distinct(MAE$n)))

# Report any errors
if (nrow(errors) > 0) {
  cat(sprintf("\n*** WARNING: %d n-replicate combinations encountered errors ***\n", nrow(errors)))
  cat("\nError details:\n")
  print(errors[, c("n", "replicate", "error_message")])
  cat("\n")
} else {
  cat("No errors encountered during processing\n")
}

# Reshape data to long format for easier plotting
MAE_long <- MAE %>%
  select(-error_message) %>%
  pivot_longer(
    cols = c(isfd, osfd, wangorig),
    names_to = "method",
    values_to = "mae"
  ) %>%
  mutate(
    method = case_when(
      method == "osfd" ~ "Failure-aware OSFD",
      method == "isfd" ~ "LHS",
      method == "wangorig" ~ "Wang-like OSFD",
      TRUE ~ method
    ),
    method = factor(method, levels = c("LHS", "Failure-aware OSFD", "Wang-like OSFD"))
  )

# Calculate central 90% region (5th to 95th percentile) and median for each method-n combo
MAE_summary <- MAE_long %>%
  group_by(n, method) %>%
  summarize(
    median_mae = median(mae, na.rm = TRUE),
    lower_90 = quantile(mae, 0.05, na.rm = TRUE),
    upper_90 = quantile(mae, 0.95, na.rm = TRUE),
    .groups = "drop"
  )

cat(sprintf("Computed central 90%% regions for %d method-n combinations\n", nrow(MAE_summary)))

# Define consistent color mapping (matching other plots)
viridis_colors <- scales::viridis_pal(option = "D")(5)
method_colors <- c(
  "LHS" = viridis_colors[1],                     # Purple/dark blue
  "Failure-aware OSFD" = viridis_colors[2],      # Blue/cyan
  "Wang-like OSFD" = viridis_colors[4]           # Yellow/green
)

p1 <- ggplot() +
  # Add shaded ribbon for central 90% region
  geom_ribbon(
    data = MAE_summary,
    aes(x = n, ymin = lower_90, ymax = upper_90, fill = method, group = method),
    alpha = 0.25
  ) +
  # Add median line
  geom_line(
    data = MAE_summary,
    aes(x = n, y = median_mae, color = method, group = method),
    linewidth = 1.2
  ) +
  # Add points at median
  geom_point(
    data = MAE_summary,
    aes(x = n, y = median_mae, color = method, group = method),
    size = 3
  ) +
  scale_color_manual(values = method_colors, name = "Method") +
  scale_fill_manual(values = method_colors, name = "Method") +
  theme_classic(base_size = 16) +
  xlab('Number of Samples') +
  ylab('Prediction MAE') +
  labs(title = 'PIV Emulator Out of Sample MAE') +
  theme(
    legend.position = 'bottom',
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 22, face = 'bold'),
    legend.key.size = unit(1.5, 'lines'),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 20, face = 'bold'),
    plot.title = element_text(size = 22, face = 'bold', hjust = 0.5)
  )

# Calculate paired differences between Failure-aware OSFD and Wang-like OSFD
mae_differences <- MAE %>%
  select(n, replicate, osfd, wangorig) %>%
  # Remove rows where either method is missing
  filter(!is.na(osfd) & !is.na(wangorig)) %>%
  # Calculate signed difference (Failure-aware - Wang-like)
  mutate(
    mae_diff = osfd - wangorig
  )

cat(sprintf("\nComputed MAE differences for %d paired observations\n", nrow(mae_differences)))

# Create difference boxplots
p2 <- ggplot(mae_differences, aes(x = factor(n), y = mae_diff)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_boxplot(fill = "gray80", color = "black", outlier.shape = 1) +
  labs(
    x = "Number of Samples",
    y = "MAE Difference\n(Failure-aware - Wang-like)"
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 20, face = 'bold'),
    legend.position = "none"
  )

# Combine the two plots side by side using patchwork
p_combined <- p1 + p2 +
  plot_layout(ncol = 2, widths = c(1, 1), guides = "collect") &
  theme(legend.position = "bottom")

scale = 3.5
ggsave(
  filename = here::here("viz", "emu_mae.png"),
  plot = p_combined,
  width = 8*scale,      # doubled width for 1x2 layout
  height = 3*scale,
  dpi = 300
)





