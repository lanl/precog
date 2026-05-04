###########################
### Create Location LUT ###
###########################


library(ggplot2)

savepath <- "~/GitLab/frodo/raw_data/"
FILES = list.files(savepath)


LOCATIONS = NULL
for(i in 1:length(FILES)){
  data_list = readRDS(paste0(savepath, FILES[i]))
  LOCATIONS_SUB = NULL
  for(j in 1:length(data_list)){
    LOCATIONS_SUB = c(LOCATIONS_SUB, data_list[[j]]$ts_geography)
  }
  LOCATIONS = unique(c(LOCATIONS, LOCATIONS_SUB))
  gc()
}

LOCATIONS_CLEANED = gsub('_NA','',LOCATIONS)
LOCATIONS_CLEANED = gsub('_null','',LOCATIONS_CLEANED)
LOCATIONS_CLEANED = gsub('_',' ', LOCATIONS_CLEANED)

LOCATIONS_TO_MERGE = unique(LOCATIONS_CLEANED)

# 
# map_world <- map_data("world")
# 
# ### Organize Names 
# map_world$region[!is.na(map_world$region)] = gsub("'","", as.character(map_world$region[!is.na(map_world$region)]))
# map_world$subregion[!is.na(map_world$subregion)] = gsub("'","", as.character(map_world$subregion[!is.na(map_world$subregion)]))
# map_world$region[!is.na(map_world$region)] = gsub("Ostrov ","", as.character(map_world$region[!is.na(map_world$region)]))
# map_world$subregion[!is.na(map_world$subregion)] = gsub("Ostrov ","", as.character(map_world$subregion[!is.na(map_world$subregion)]))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# cbind(LOCATIONS[1001:length(LOCATIONS)])