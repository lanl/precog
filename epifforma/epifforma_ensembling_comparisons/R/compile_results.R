# Compile the results of the parallelization for the covidHub comparison of sMOA.
## Author: Alexander C. Murph
## Date: July 2024

compile_results = function(dir_name, verbose=T){
	total_wins <- c()
	average_totals <- c()
	wins_by_model <- list()
	wins_by_model_wis <- list()
	my_model_wins_wis <- c()
	total_wins_wis <- c()

	scores                          <- read.csv("data/scores_tot.csv")

	model_names <- unique(scores$model)
	for(curr_model_name in model_names){
	  name = paste('justresults_', curr_model_name, sep = "")
	  wins_by_model[[name]] <- c()
	  wins_by_model_wis[[name]] <- c()
	}
	
	for(file_name in list.files(dir_name)){
		results_list <- NULL
		load(paste(dir_name, "/", file_name, sep=""))
		
		model_name = unlist(strsplit(file_name, split = ".R"))[1]
		for(curr_model_name in model_names){
		  name = paste('justresults_', curr_model_name, sep = "")
		  #if( (length(results_list[[name]]$moa_wis)>0)&(length(results_list[[name]]$covidhub_wis)>0) ){
		    wins_by_model_wis[[name]] <- c(wins_by_model_wis[[name]], results_list[[name]]$moa_wis < results_list[[name]]$covidhub_wis)
		  #}

		  #if( (length(results_list[[name]]$moa_abs)>0)&(length(results_list[[name]]$covidhub_abs)>0) ){
                    wins_by_model[[name]] <- c(wins_by_model[[name]], results_list[[name]]$moa_abs < results_list[[name]]$covidhub_abs)
                  #}
		}

		if(verbose){
  		print("-----------------------")
  		print(paste("Average wins wrt MAE for", model_name, "is", mean(results_list[['wins_vector_abs']]), "over a total of", length(results_list[['wins_vector_abs']]), "model comparisons"))
          	print(paste("The probability that sMOA would beat this many models as a random coin flip is", pbinom(q = mean(results_list[['wins_vector_abs']])*length(results_list[['wins_vector_abs']]), 
  														     prob = 0.5, size = length(results_list[['wins_vector_abs']]), lower.tail = F)))

		print(paste("Average wins wrt MAE for", model_name, "is", mean(results_list[['wins_vector_wis']]), "over a total of", length(results_list[['wins_vector_wis']]), "model comparisons"))
                print(paste("The probability that sMOA would beat this many models as a random coin flip is", pbinom(q = mean(results_list[['wins_vector_wis']])*length(results_list[['wins_vector_wis']]),
                                                                                                                     prob = 0.5, size = length(results_list[['wins_vector_wis']]), lower.tail = F)))
  		print("-----------------------")
		}

		total_wins = c(total_wins, results_list[['wins_vector_abs']])
		average_totals = c(average_totals, results_list[['proportion_wins']])
		total_wins_wis <- c(total_wins_wis, results_list[['wins_vector_wis']])
	}

	if(verbose) print("The MAE wins by model are:")
	for(curr_model_name in model_names){
	  name = paste('justresults_', curr_model_name, sep = "")
	  print(paste(name, ": ", mean(wins_by_model[[name]]), sep=""))
	}
	if(verbose) print("-----------------------")
	if(verbose) print("The WIS wins by model are:")
	
	for(curr_model_name in model_names){
	  name = paste('justresults_', curr_model_name, sep = "")
	  if(verbose) print(paste(name, ": ", mean(wins_by_model_wis[[name]]), sep=""))
	}
	if(verbose) print("-----------------------")
	
	if(verbose) {
  	print(paste("Total average wins is", mean(total_wins), "over a total of", length(total_wins), "model comparisons"))
  	print(paste("The probability that sMOA would beat this many models as a random coin flip is", pbinom(q = mean(total_wins)*length(total_wins), prob = 0.5, size = length(total_wins), lower.tail = F)))
  	print(paste("This same value for wis is", mean(total_wins_wis), "over a total of", length(total_wins_wis), "model comparisons"))
  	print(paste("The probability that sMOA would beat this many models as a random coin flip is", 
  	            pbinom(q = mean(total_wins_wis)*length(total_wins_wis), prob = 0.5, size = length(total_wins_wis), lower.tail = F)))
	}
	return(data.frame(average_wins = mean(total_wins), total_comparisons = length(total_wins),
	       total_wins_wis = mean(total_wins_wis)))
}


