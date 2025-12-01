# Script to print out initial evaluation of sMOA study.

setwd("~/GitLab/smoa")

full_data = NULL
for(file_name in list.files('data/results')){
  temp_data                  <- read.csv(paste("data/results/", file_name, sep = ""))
  split_string               <- unlist(strsplit(file_name, split = "_"))
  if(split_string[1]=='wHistories') next
  temp_data$k = as.numeric(split_string[2])
  temp_data$num_curves       <- as.numeric(split_string[5])
  temp_data$how_many_closest <- as.numeric(unlist(strsplit(split_string[7],
                                                           split=".c"))[1])
  if(ncol(temp_data) < 12) next
  print(temp_data)
  full_data = rbind(full_data, temp_data)
}
full_data$X                  <- NULL
print(full_data)
summary(full_data)

# Now let's print out the quick and dirty by model:
#model_names = c("BPagano-RtDriven","CEID-Walk","COVIDhub-baseline","COVIDhub-4_week_ensemble","COVIDhub-trained_ensemble","Covid19Sim-Simulator",
#                "CU-select","JHUAPL-Bucky","JHU_CSSE-DECOM","JHU_IDD-CovidSP","Karlen-pypm","LNQ-ens1",
#                "LANL-GrowthRate","Microsoft-DeepSTIA","MOBS-GLEAM_COVID","RobertWalraven-ESG","UVA-Ensemble")
#full_data = NULL
#for(file_name in list.files('data')){
#  split_string <- unlist(strsplit(file_name, split = "_"))
#  if(split_string[1]=="k"){
#    for(state_file_name in list.files(paste('data/', file_name, sep = ""))){
#      load(paste('data/', file_name, "/",state_file_name, sep = ""))
#      #browser()
#    }
#  }
#
#
#}
#full_data$X                  <- NULL
#print(full_data)
#summary(full_data)


