# Look through output files of Mutantigen Parallel to see what finished and what 
# did not finish.
## Author: AC Murph
## Date: January 2025
library(ggplot2)
library(gridExtra)
library(ggExtra)
library(parallel)
library(doParallel)
library(this.path)
setwd(paste0(this.path::here(), '/../'))
# Number of iterations performed in this experiment
num_iterations = length(list.files('input_files'))

# Iterate through all outfiles expected
all_outfiles = list.files('outfiles/')
all_outfiles_split = strsplit(all_outfiles, split = '\\.')
has_all_files = c()
has_timeseries_not_trees = c()
num_of_files = c()
full_file_df = NULL

parfctn = function(curr_iter){
# for(curr_iter in 1:num_iterations){
  # print(paste("Investigating output file", curr_iter, "of", num_iterations))
  curr_pattern = paste0("out_", curr_iter)
  curr_pattern_2 = paste0("out", curr_iter)
  full_file_df_instance = NULL
  
  # Check to see if all 20 files are present.
  if( (paste0(curr_pattern, '.antigenicDistances')%in%all_outfiles)&
      (paste0(curr_pattern, '.antigenicSamples')%in%all_outfiles)&
      (paste0(curr_pattern, '.branches')%in%all_outfiles)&
      (paste0(curr_pattern, '.mrcaSeries')%in%all_outfiles)&
      (paste0(curr_pattern, '.mutationSeries')%in%all_outfiles)&
      (paste0(curr_pattern, '.out.simmapAntigenic')%in%all_outfiles)&
      (paste0(curr_pattern, '.simmapLoad')%in%all_outfiles)&
      (paste0(curr_pattern, '.summary')%in%all_outfiles)&
      (paste0(curr_pattern, '.timeseries')%in%all_outfiles)&
      (paste0(curr_pattern, '.tips')%in%all_outfiles)&
      (paste0(curr_pattern, '.trees')%in%all_outfiles)&
      (paste0(curr_pattern, '.trunkAntigenicShifts')%in%all_outfiles)&
      (paste0(curr_pattern, '.viralFitnessSeries')%in%all_outfiles)&
      (paste0(curr_pattern_2, '.fitnessSamples_t1000.0')%in%all_outfiles)&
      (paste0(curr_pattern_2, '.fitnessSamples_t2000.0')%in%all_outfiles)&
      (paste0(curr_pattern_2, '.fitnessSamples_t3000.0')%in%all_outfiles)&
      (paste0(curr_pattern_2, '.fitnessSamples_t4000.0')%in%all_outfiles)&
      (paste0(curr_pattern_2, '.fitnessSamples_t5000.0')%in%all_outfiles)&
      (paste0(curr_pattern_2, '.fitnessSamples_t6000.0')%in%all_outfiles)&
      (paste0(curr_pattern_2, '.fitnessSamples_t7000.0')%in%all_outfiles)
  ) {
    has_all_files_instance = 1
  } else {
    has_all_files_instance = 0
  }
  
  # Check to see if the timeseries file is there but the trees file is not.
  if( (paste0(curr_pattern, '.timeseries')%in%all_outfiles)&
      (!(paste0(curr_pattern, '.trees')%in%all_outfiles))
  ) {
    has_timeseries_not_trees_instance = 1
  } else {
    has_timeseries_not_trees_instance = 0
  }
  
  # Count the number of files currently for this simulation run.
  outfile_name = c()
  for(inner_iter in 1:length(all_outfiles)){
    if(all_outfiles_split[[inner_iter]][1] == curr_pattern | 
       all_outfiles_split[[inner_iter]][1] == curr_pattern_2){
      if(length(all_outfiles_split[[inner_iter]])==3) {
        if(all_outfiles_split[[inner_iter]][3]=="0"){
          tmp_file_name = all_outfiles_split[[inner_iter]][2]
        } else {
          tmp_file_name = all_outfiles_split[[inner_iter]][3]
        }
      }
      if(length(all_outfiles_split[[inner_iter]])==2) {
        tmp_file_name = all_outfiles_split[[inner_iter]][2]
      }
      
      outfile_name = c(outfile_name, tmp_file_name)
      full_file_df_instance = data.frame(iter_num = curr_iter,
                                      file_name = tmp_file_name)
    } 
  }
  
  return(list(num_of_files_instance = length(outfile_name),
              has_all_files_instance = has_all_files_instance,
              has_timeseries_not_trees_instance = has_timeseries_not_trees_instance,
              full_file_df_instance = full_file_df_instance
        ))
}

sockettype <- "PSOCK"
ncores <- 10
cl <- parallel::makeCluster(spec = ncores,type = sockettype) #, outfile=""
setDefaultCluster(cl)
registerDoParallel(cl)
sim_ts <- foreach(i=1:num_iterations,
                  .verbose = T)%dopar%{
                    parfctn(i)
                  }
stopCluster(cl)    

for(i in 1:length(sim_ts)){
  num_of_files = c(num_of_files, sim_ts[[i]]$num_of_files_instance)
  has_all_files = c(has_all_files, sim_ts[[i]]$has_all_files_instance)
  has_timeseries_not_trees = c(has_timeseries_not_trees, sim_ts[[i]]$has_timeseries_not_trees_instance)
  full_file_df = rbind(full_file_df,
                       sim_ts[[i]]$full_file_df_instance)
}


p1 = ggplot(full_file_df, aes(x = file_name)) + 
  geom_bar(aes(y = after_stat(count) / 100), fill = "steelblue") + theme_bw()+                              
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  ylab('proportion of sims \nthat had this file') + xlab('file name')

p2 = ggplot(data.frame(x = as.factor(num_of_files)), aes(x = x)) + 
  geom_bar(fill = "steelblue") + theme_bw()+                              
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  ylab('number of sims that \nproduced this many outfiles') + 
  xlab('number of outfiles produced by a sim')

png("mutantigen_outfiles_analysis.png")
grid.arrange(p1,p2, layout_matrix = matrix(1:2))
dev.off()

