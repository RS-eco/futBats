# Create prediction maps for the different taxa, GCMs, rcps and timesteps

# Empty R environment
rm(list=ls())

# Load libraries
packages <- c("tidyverse", "data.table")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages); rm(new.packages)

# Load required packages
l <- sapply(packages, require, character.only = TRUE); rm(packages, l)

# Set file directory
filedir <- "F:/" # ceremony - mammals

# Set taxa
taxa <- "Mammal"

# Model type
k <- 1; model_type <- c("GAM", "GBM")[k]

# Final species names
spNames <- sapply(list.files(paste0(filedir, "/", taxa, "_", 
                                    model_type, "_results_climate")), function(x){
                                      paste0(strsplit(x, split="_")[[1]][1:2], collapse="_")
                                    })

# Load bat file
bat_iucn <- readr::read_csv("data/IUCN_All_bat-species_Code_01.csv.xz") %>% select(-X1)
head(bat_iucn)
length(unique(bat_iucn$layer))

batNames <-  sub(" ", "_", unique(bat_iucn$layer))
sp_bat <- spNames[spNames %in% batNames]

# Final prediction files
predFiles <- list.files(paste0(filedir, "/", taxa, "_", model_type, "_results_climate"), 
                        full.names=T)
predFiles <- sapply(sp_bat, function(x) predFiles[grepl(predFiles, pattern=x)][[1]])
length(predFiles)

# Read csv files
library(snowfall)
sfInit(parallel=TRUE, cpus=ceiling(0.95*parallel::detectCores()))
predData <- sfLapply(predFiles, function(x){readr::read_csv(x)})
sfStop()
#predData <- dplyr::bind_rows(predData)
predData <- data.table::rbindlist(predData)

# Calculate summed probability and save to .csv.xz file
lapply(c("presence", "dispersal1", "dispersal2", "dispersal3", "dispersal4", "fulldisp"), function(disp){
  if(!file.exists(paste0("data/bat_prob_", model_type, "_", disp, ".csv.xz"))){
    if(disp != "fulldisp"){
      sumProb_nodisp <- predData %>% group_by(x, y) %>% filter(1 == !!as.name(disp)) %>% 
        select(x,y,GFDL.ESM2M_piControl_1845:MIROC5_rcp26_2250) %>% 
        mutate_if(is.character, as.numeric) %>% 
        summarise_all(sum, na.rm=T) %>% as.data.frame()
      readr::write_csv(sumProb_nodisp, paste0("data/bat_prob_", 
                                              model_type, "_", disp, ".csv.xz"))
    } else{
      sumProb_disp <- predData %>% group_by(x, y) %>% 
        select(x,y,GFDL.ESM2M_piControl_1845:MIROC5_rcp26_2250) %>% 
        mutate_if(is.character, as.numeric) %>% 
        summarise_all(sum, na.rm=T) %>% as.data.frame()
      readr::write_csv(sumProb_disp, paste0("data/bat_prob_", 
                                            model_type, "_fulldisp.csv.xz"))
    }
  }
  return(NULL)
}); q(save="no")
