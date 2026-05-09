# Look through the outfiles to see which ones did and did not complete, color input parameter
# space based on this to see if there is a pattern.
## Author: AC Murph
## Date 08/05/2025
library(ggplot2)
library(gridExtra)
library(ggExtra)
library(parallel)
library(doParallel)
library(tidyverse)
library(this.path)
library(data.table)
library(phytools)
library(parallel)
library(doParallel)
setwd(paste0(this.path::here(), '/../'))
source("R/emma_helpers.R")

# Number of iterations performed in this experiment
num_iterations = length(list.files('input_files'))


#### Dave output file analysis stuff:
# Iterate through all outfiles expected
all_outfiles = list.files('outfiles/')
all_outfiles_split = strsplit(all_outfiles, split = '\\.')
has_all_files = c()
has_timeseries_not_trees = c()
num_of_files = c()
full_file_df = NULL
lhs_samples = as_tibble(read.csv('lhs_samples.csv'))
lhs_samples$sim_number = 1:nrow(lhs_samples)

# go through tabulate_df and find when prop_non_zero_cases < 0.01
tabulate_df = as_tibble(read.csv('data_post_processing/tabulate_df.csv'))
# only_csv = tabulate_df %>%
#   subset(csv == TRUE) %>%
#   select(prop_non_zero_cases, common_name)
# only_csv$has_cases = only_csv$prop_non_zero_cases > 0.01
# only_csv$sim_number = as.numeric(sub(".*_(\\d+)", "\\1", only_csv$common_name))

lhs_merged = lhs_samples #%>%
#   left_join(only_csv, by = "sim_number")
lhs_merged = lhs_merged[order(lhs_merged$sim_number),]

#### Emma output file analysis stuff:
# Iterate through everything represented in lhs_merged:
lhs_merged$has_both_outfiles = FALSE
lhs_merged$has_atleast_five_cases = FALSE
lhs_merged$has_atleast_five_seq = FALSE
source_path = this.path::here()

parfctn = function(runnum){
  runnum = as.character(runnum)
  outdir <- paste0(source_path, '/../outfiles')
  
  mafile <- file.path(outdir, paste0("out_", runnum,
                                     ".out.simmapAntigenic"))
  nxfile <- file.path(outdir, paste0("out_", runnum,
                                     ".simmapAntigenic.nex"))
  tpfile <- file.path(outdir, paste0("out_", runnum, ".tips"))
  csfile <- file.path(outdir, paste0("out_", runnum, ".timeseries"))
  
  if(!file.exists(csfile) | !file.exists(mafile)){
    return(list(
      has_both_outfiles = FALSE,
      has_atleast_five_cases = FALSE,
      has_atleast_five_seq = FALSE
    ))
  }
  
  tree_to_nexus(mafile, nxfile)
  ans <- tips_to_bins(tpfile, nxfile, csfile)
  write_csv(ans, paste0("data_post_processing/count_files/counts", runnum, ".csv"))
  
  # drop the first rows, where tips are not yet sampled
  ans1 <- ans
  while(anyNA(ans1[1,]) == T)
    ans1 <- ans1[-1,]
  
  nt <- nrow(ans1)
  
  # how many times have at least 5 cases?
  c5 <- ans1 %>% filter(Cases > 4) %>% nrow
  has_atleast_five_cases = c5 / nt > 0.9 # want this to be true
  
  # how many times have at least 5 sequences?
  s5 <- ans1 %>% filter(NumSeq > 4) %>% nrow
  has_atleast_five_seq = s5 / nt > 0.9 # want this to be true
  
  return(list(
    has_both_outfiles = TRUE,
    has_atleast_five_cases = has_atleast_five_cases,
    has_atleast_five_seq = has_atleast_five_seq
  ))
}

# java -cp dist/MutAntiGen.jar:/lib/colt.jar Mutantigen "parameters_load.yml" 1 
ncores <- 10
sockettype <- "FORK"
cl <- parallel::makeCluster(spec = ncores,type = sockettype, outfile = "") #, outfile=""
setDefaultCluster(cl)
registerDoParallel(cl)
sim_ts <- foreach(i=1:nrow(lhs_merged),
                  .verbose = T)%dopar%{
                  library(ggplot2)
                  library(gridExtra)
                  library(ggExtra)
                  library(parallel)
                  library(doParallel)
                  library(tidyverse)
                  library(this.path)
                  library(data.table)
                  library(phytools)
                  library(parallel)
                  library(doParallel)
                  # setwd(paste0(this.path::here(), '/../'))
                  source("R/emma_helpers.R")
                    parfctn(i)	    
                  }
stopCluster(cl)

for(sim_num in 1:length(sim_ts)){
  lhs_merged$has_both_outfiles[sim_num] = sim_ts[[sim_num]]$has_both_outfiles
  lhs_merged$has_atleast_five_cases[sim_num] = sim_ts[[sim_num]]$has_atleast_five_cases
  lhs_merged$has_atleast_five_seq[sim_num] = sim_ts[[sim_num]]$has_atleast_five_seq
}

write.csv(lhs_merged, file = 'data_post_processing/lhs_merged.csv')


# Now graph all of the data.
lhs_merged = as_tibble(read.csv(file = 'data_post_processing/lhs_merged.csv'))
lhs_merged$X = NULL
credible_marginals = lhs_merged %>%
  subset(!is.na(has_cases))%>%
  GGally::ggpairs(
    columns = 1:14,
    aes(color = has_cases, alpha = 0.4)
) + ggtitle("Variables when both files are present, but there are <0.01 proportion of cases")

ggsave(
  filename = "data_post_processing/low_proportion_cases.png",
  plot = credible_marginals,
  width = 14*1.5,
  height = 14*1.5,
  dpi = 300 * 2
)



lhs_filtered <- lhs_merged %>%
  filter(!is.na(has_both_outfiles)) %>%
  mutate(has_both_outfiles = factor(has_both_outfiles, levels = c(TRUE, FALSE)))

# Create the plot
credible_marginals <- GGally::ggpairs(
  lhs_filtered,
  columns = 1:14,
  aes(color = has_both_outfiles, alpha = 0.4)
) +
  ggtitle("Variables by when both outfiles are present")

ggsave(
  filename = "data_post_processing/when_outfiles_present.png",
  plot = credible_marginals,
  width = 14*1.5,
  height = 14*1.5,
  dpi = 300 * 2
)

lhs_filtered <- lhs_merged %>%
  filter(!is.na(has_both_outfiles)) %>%
  mutate(has_both_outfiles = factor(has_both_outfiles, levels = c(FALSE, TRUE)))

# Create the plot
credible_marginals <- GGally::ggpairs(
  lhs_filtered,
  columns = 1:14,
  aes(color = has_both_outfiles, alpha = 0.4)
) +
  ggtitle("(reversed) Variables by when both outfiles are present")

scale = 1.25
ggsave(
  filename = "data_post_processing/REVERSEDwhen_outfiles_present.png",
  plot = credible_marginals,
  width = 14 * scale,
  height = 14 * scale,
  dpi = 300 * scale
)




credible_marginals = lhs_merged %>%
  subset(!is.na(has_atleast_five_cases) & has_both_outfiles) %>%
  GGally::ggpairs(
        columns = 1:14,
        aes(color = has_atleast_five_cases, alpha = 0.4)
) + ggtitle("Variables when both files are present, but there are >=5 cases")

ggsave(
  filename = "data_post_processing/low_count_cases.png",
  plot = credible_marginals,
  width = 14*scale,
  height = 14*scale,
  dpi = 300 * scale
)

credible_marginals = lhs_merged %>%
  subset(!is.na(has_atleast_five_cases) & has_both_outfiles) %>%
  GGally::ggpairs(
      columns = 1:14,
      aes(color = has_atleast_five_seq, alpha = 0.4)
) + ggtitle("Variables when both files are present, but there are >=5 sequences")

ggsave(
  filename = "data_post_processing/low_count_sequences.png",
  plot = credible_marginals,
  width = 14*scale,
  height = 14*scale,
  dpi = 300 * scale
)


lhs_merged %>% subset(has_both_outfiles & has_atleast_five_cases & has_atleast_five_seq) %>% head(1)
lhs_merged %>% subset(has_both_outfiles & has_atleast_five_cases )%>% head(1)
lhs_merged %>% subset(!has_both_outfiles)%>% head(1)
lhs_merged %>% subset(has_both_outfiles & !has_atleast_five_cases)%>% head(1)

