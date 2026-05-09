# Mapping Input files to replicate number
## AC Murph
## Date: Feb 2026

library(tidyr)
library(this.path)

setwd(paste0(this.path::here(), '/../'))

# java -cp dist/MutAntiGen.jar:/lib/colt.jar Mutantigen "parameters_load.yml" 1 
ncores <- 99
sockettype <- "PSOCK"

number_of_reps = 20

good_input_files = as_tibble(read.csv('antipeak.csv'))

hey_murph_what_input_file_was_this_output_id_for = function(output_file_id){
  i = output_file_id
  num_of_yaml_row = (i-1) %% nrow(good_input_files) + 1
  num_of_yaml = good_input_files$runnum[num_of_yaml_row]
  print(paste("That output id corresponds to input file number:", num_of_yaml))

  name_of_yaml = paste0("input_files/parameters_load_updated_", num_of_yaml, ".yml")
  print(name_of_yaml)
  return(num_of_yaml)
}
# hey_murph_what_input_file_was_this_output_id_for(300)
# [1] "That output id corresponds to input file number: 218"
# [1] "input_files/parameters_load_updated_218.yml"
# [1] 218

hey_murph_give_me_all_replicates_that_correspond_to_this_input_file = function(input_file_id){
  outfiles_list = c()
  for( i in 1:(nrow(good_input_files)*number_of_reps) ){
    num_of_yaml_row = (i-1) %% nrow(good_input_files) + 1
    num_of_yaml = good_input_files$runnum[num_of_yaml_row]
    name_of_yaml = paste0("input_files/parameters_load_updated_", num_of_yaml, ".yml")

    if(num_of_yaml == input_file_id){
      file1 = paste0("outfiles/out_", i, ".timeseries")
      file2 = paste0("outfiles/out_", i, ".out.simmapAntigenic")
      if(file.exists(file1) & file.exists(file2)){
        outfiles_list = c(outfiles_list, i)
      }
    }
  }

  print(paste("The outfile ids that correspond to that input file id are:"))
  print(outfiles_list)
  return(outfiles_list)
}

# hey_murph_give_me_all_replicates_that_correspond_to_this_input_file(3787)
# [1] "The outfile ids that correspond to that input file id are:"
#  [1]    1  153  305  457  609  913 1065 1217 1369 1673


get_number_of_files_completed = function(){
  count = 0
  for( i in 1:(nrow(good_input_files)*number_of_reps) ){
    file1 = paste0("outfiles/out_", i, ".timeseries")
    file2 = paste0("outfiles/out_", i, ".out.simmapAntigenic")
    if(file.exists(file1) & file.exists(file2)){
      count = count + 1
    }
  }
  print(paste(count, "files completed out of", nrow(good_input_files)*number_of_reps))
  print(paste("Propostion completed:", count / (nrow(good_input_files)*number_of_reps) ))
}

# get_number_of_files_completed()
