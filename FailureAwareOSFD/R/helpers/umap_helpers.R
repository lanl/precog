# depends: 
#' @export
create_umap_df = function(DISEASE_KEY, FILE_NAME, Y_ISFD, Y_OSFD){
  suppressPackageStartupMessages(library(umap))
  suppressPackageStartupMessages(library(patchwork))
  suppressPackageStartupMessages(library(dplyr))
  ### demonstrating using Global_Covid
  tmp = as_tibble(readRDS(paste0("data/processed_features_data/", FILE_NAME)))
  if(nrow(tmp) > 10000){
    tmp <- tmp[sample(1:nrow(tmp), 10000, replace = F),]
  }
  tmp = tmp[!is.na(tmp$entropy),]
  train_disease_data = tmp %>%
    # subset(disease_name == DISEASE_KEY) %>%
    dplyr::select(gr12_div_23, last_div_max, coefvar, 
                  gam_with_div_without, avg_recent_div_avg_global,
                  diff_zscore, entropy, relative_increases,
                  prop_since_peak, seasonality)%>%
    dplyr::distinct()
  
  ## combine with synthetic
  all_other_data = rbind(Y_ISFD, Y_OSFD)
  
  global_coviddf <- as.matrix(rbind(all_other_data, train_disease_data))
  
  ## fit UMAP (2-dimensions)
  set.seed(1128)
  global_covid_umap <- umap(global_coviddf, n_components = 2)
  global_covid_umap_df <- data.frame(global_covid_umap$layout)
  names(global_covid_umap_df) <- paste0("X",1:ncol(global_covid_umap_df))
  # global_covid_umap_df = screen_multivariate_outliers(global_covid_umap_df, threshold = outlier_threshold)
  global_covid_umap_df$type <- DISEASE_KEY
  global_covid_umap_df[1:nrow(all_other_data),]$type <- "OSFD"
  global_covid_umap_df[1:nrow(Y_ISFD),]$type <- "ISFD"
  global_covid_umap_df$is_synthetic = TRUE
  global_covid_umap_df[(nrow(all_other_data) + 1):nrow(global_coviddf),]$is_synthetic <- FALSE
  return(global_covid_umap_df)
}

#' @export
create_umap_plot = function(umap_df, disease_name, data_color = "#fb8072"){
  rubella_umap_df = umap_df
  p1 <- ggplot() +
    geom_point(aes(x = X1, y = X2),
               data = subset(rubella_umap_df, type == "OSFD"),
               color = I("black")) +
    geom_point(aes(x = X1, y = X2),
               data = subset(rubella_umap_df, is_synthetic == FALSE),
               color = I(data_color)) +
    ggtitle("") + theme_bw() +
    xlim(range(rubella_umap_df$X1)) +
    ylim(range(rubella_umap_df$X2)) +
    xlim(-15, 15) +
    ylim(-15, 15) +
    theme(plot.title = element_text(hjust = 0.5)) + 
    ylab("OSFD Sampled Synthetic") + xlab("")+
    theme(
      axis.title   = element_text(size = 20),
      axis.text    = element_text(size = 14),
      plot.title   = element_text(size = 18),
      plot.subtitle = element_text(size = 14),
      legend.title = element_text(size = 14),
      legend.text  = element_text(size = 12)
    )
  
  p2 <- ggplot()  +
    geom_point(aes(x = X1, y = X2),
               data = subset(rubella_umap_df, is_synthetic == FALSE),
               color = I(data_color)) +
    geom_point(aes(x = X1, y = X2),
               data = subset(rubella_umap_df, type == "OSFD"),
               color = I("black")) +
    ggtitle("") + theme_bw() +
    xlim(range(rubella_umap_df$X1)) +
    ylim(range(rubella_umap_df$X2)) +
    xlim(-15, 15) +
    ylim(-15, 15)+
    theme(plot.title = element_text(hjust = 0.5))+ 
    ylab("") + xlab("")+
    theme(
      axis.title   = element_text(size = 20),
      axis.text    = element_text(size = 14),
      plot.title   = element_text(size = 18),
      plot.subtitle = element_text(size = 14),
      legend.title = element_text(size = 14),
      legend.text  = element_text(size = 12)
    )
  
  p3 <- ggplot() +
    geom_point(aes(x = X1, y = X2),
               data = subset(rubella_umap_df, type == "ISFD"),
               color = I("black")) +
    geom_point(aes(x = X1, y = X2),
               data = subset(rubella_umap_df, is_synthetic == FALSE),
               color = I(data_color)) +
    ggtitle("") + theme_bw() +
    xlim(range(rubella_umap_df$X1)) +
    ylim(range(rubella_umap_df$X2)) +
    xlim(-15, 15) +
    ylim(-15, 15)+
    theme(plot.title = element_text(hjust = 0.5))+ 
    xlab(paste(disease_name, "Data on Top")) + ylab("ISFD Sampled Synthetic")+
    theme(
      axis.title   = element_text(size = 20),
      axis.text    = element_text(size = 14),
      plot.title   = element_text(size = 18),
      plot.subtitle = element_text(size = 14),
      legend.title = element_text(size = 14),
      legend.text  = element_text(size = 12)
    )
  
  p4 <- ggplot()  +
    geom_point(aes(x = X1, y = X2),
               data = subset(rubella_umap_df, is_synthetic == FALSE),
               color = I(data_color)) +
    geom_point(aes(x = X1, y = X2),
               data = subset(rubella_umap_df, type == "ISFD"),
               color = I("black")) +
    ggtitle("") + theme_bw() +
    xlim(range(rubella_umap_df$X1)) +
    ylim(range(rubella_umap_df$X2)) +
    xlim(-15, 15) +
    ylim(-15, 15)+
    theme(plot.title = element_text(hjust = 0.5))+ 
    ylab("") + xlab("Synthetic Data on Top")+
    theme(
      axis.title   = element_text(size = 20),
      axis.text    = element_text(size = 14),
      plot.title   = element_text(size = 18),
      plot.subtitle = element_text(size = 14),
      legend.title = element_text(size = 14),
      legend.text  = element_text(size = 12)
    ) 
  
  # 2×2 grid
  pp = (p1 | p2) /
    (p3 | p4) +
    plot_annotation(
      title = paste(disease_name, "UMAP"),
      theme = theme(
        plot.title = element_text(
          hjust = 0.5,
          size = 20,
          face = "bold"
        )
      )
    )
  return(pp)
  
}