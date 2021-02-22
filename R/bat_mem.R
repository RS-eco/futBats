#' ---
#' title: "R-Code for MEM models of Bats"
#' author: "RS-eco"
#' ---

#' setup, eval=F
#Load required packages

rm(list=ls()); invisible(gc())

#Automatically install required packages, which are not yet installed
packages <- c("sp", "raster", "tidyverse", "ggpmisc", "ModelMetrics",
              "gbm", "mgcv", "magrittr", "devtools")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages); rm(new.packages)

#Load packages
l <- sapply(packages, require, character.only = TRUE, quietly=TRUE); rm(packages, l)

# Set working directory
setwd("/home/matt/Documents/futBats/")

#' species_data, eval=F
#Load species data
bats_dist <- readr::read_csv("data/IUCN_All_bat-species_Code_01.csv") %>% select(-X1)

#Create presence column
bats_dist$presence <- 1
length(unique(bats_dist$layer))

#' species_richness, eval=F
#Calulate bat SR
sr_observed <- bats_dist %>% group_by(x, y) %>% 
  summarise(sum = sum(presence))

#' load_climate, eval=F
#Load climate predictor combs
climCombs <- c("bio4","bio5","bio12","bio15") 

## Load climate data and blocks
climVar <- read.csv(paste0("extdata/Blocking_SU_1995_", paste0(climCombs, collapse="_"), ".csv.xz"))[,c("x","y", climCombs, "block")]
colnames(climVar) <- c("x","y", climCombs, "block")

#' ## Macroecological models

#' mem_eco_mod, eval=F

#Set results Path
plotPath_MEM <- paste0(getwd(), "/data")
resultsPath_MEM <- paste0(getwd(), "/data")

# Create model, family, dataset combinations
model_fam <- c("GAM_poisson", "GBM_poisson")

#Load Model function
library(mgcv); library(gbm); library(PresenceAbsence); 
library(dplyr); library(snowfall)
source("R/GAM_func.R"); source("R/GBM_func.R")

#Run model for each algorithm
sfInit(parallel=TRUE, cpus=length(model_fam))
sfLibrary(mgcv); sfLibrary(dplyr); sfLibrary(gbm); sfLibrary(PresenceAbsence)
sfExport(list=c("sr_observed", "GAM_eco","GAM_split", 
                "resultsPath_MEM", "filedir", 
                "climCombs", "climVar", "GBM_eco", "GBM_split", 
                "plotPath_MEM")) 
sfLapply(model_fam, function(y){
  model_type <- strsplit(as.character(y), split="_")[[1]][1]
  family <- strsplit(as.character(y), split="_")[[1]][2]
  species.data <- sr_observed
  
  #Run model
  clim.var <- as.character(unlist(climCombs))
  if(!file.exists(paste(getwd(), resultsPath_MEM, "/", paste(clim.var, collapse="_"), "_model_output_", y, "_Eco_block.RData",sep=""))){
    ## Get species data
    spPseudoRep <- na.omit(species.data)
    colnames(spPseudoRep)[1:3] <-c("x", "y", "presence")
    spPseudoRep <- merge(spPseudoRep, climVar,by=c("x","y"), all.x=T)
    spPseudoRep <- dplyr::select(spPseudoRep, x,y,presence, 
                                 dplyr::one_of("block"), 
                                 dplyr::one_of(clim.var))
    spPseudoRep["cell.id"] <- seq(1,nrow(spPseudoRep))
    
    spPseudoRep <- na.omit(spPseudoRep)
    
    ## Model function GAM
    if(model_type == "GAM"){
      if(family == "nb"){family <- nb()}
      m <- GAM_eco(data.model=spPseudoRep, 
              outDir=resultsPath_MEM, 
              family=family, 
              clim.var=clim.var,
              fx=FALSE, k=-1, bs="tp",
              blocks=c(1:10))
    } else if(model_type == "GBM"){
      ## Model function GBM
      m <- GBM_eco(data.model=spPseudoRep, 
              outDir=resultsPath_MEM,
              plotPath=plotPath_MEM, eval=FALSE,
              learn.rate=c(0.05, 0.01),
              distribution=family,
              clim.var=clim.var,
              blocks=c(1:10))
    }
  }
})
sfStop()

#' ## Runs predictions for 1995 and 2080!

#' mem_eco_pred, eval=F

# Create model, family, dataset combinations
model_fam <- c("GAM_poisson", "GBM_poisson")

#Run model for each algorithm
library(snowfall); library(mgcv); library(dplyr); library(gbm)
sfInit(parallel=TRUE, cpus=nrow(df))
sfLibrary(mgcv); sfLibrary(dplyr); sfLibrary(gbm)
sfLapply(model_fam, function(y){
  if(!file.exists(paste0("data/sr_predicted_mem_", y, "_eco_2080.csv.xz"))){
    # Read models into R
    m <- get(load(paste0("data/bio4_bio5_bio12_bio15_model_output_", y,"_Eco_block.RData")))
    if(!file.exists(paste0("data/sr_predicted_mem_", y, "_eco_1995.csv.xz"))){
      # Run current predictions
      data <- readr::read_csv("extdata/bioclim_EWEMBI_1995_landonly.csv.xz")
      
      # Predict function depends on mgcv library
      pred <- lapply(1:length(m), function(x){
        predict <- round(as.numeric(
          predict(m[[x]]$mod, newdata=data, 
                  type="response",se.fit=FALSE)),4)
        predict <- as.data.frame(cbind(data[,c("x","y")], predict))
        colnames(predict) <- c("x", "y", x)
        return(predict)
      })
      pred <- Reduce(function(...) merge(..., by=c("x","y"), all.x=T), pred)
      
      # Calculate mean prediction
      pred$mean <- rowMeans(pred[,c(3:12)])
      
      # Merge predicted richness
      pred$model <- y
      pred$taxa <- "bats"
      readr::write_csv(pred, paste0("data/sr_predicted_mem_", y, "_eco_1995.csv.xz")); rm(pred)
    }
    
    # Select climate data
    climatenames<- c("extdata/bioclim_GFDL-ESM2M_rcp26_2080_landonly.csv.xz",
                     "extdata/bioclim_GFDL-ESM2M_rcp60_2080_landonly.csv.xz", 
                     "extdata/bioclim_GFDL-ESM2M_rcp85_2080_landonly.csv.xz", 
                     "extdata/bioclim_HadGEM2-ES_rcp26_2080_landonly.csv.xz", 
                     "extdata/bioclim_HadGEM2-ES_rcp60_2080_landonly.csv.xz", 
                     "extdata/bioclim_HadGEM2-ES_rcp85_2080_landonly.csv.xz", 
                     "extdata/bioclim_IPSL-CM5A-LR_rcp26_2080_landonly.csv.xz", 
                     "extdata/bioclim_IPSL-CM5A-LR_rcp60_2080_landonly.csv.xz", 
                     "extdata/bioclim_IPSL-CM5A-LR_rcp85_2080_landonly.csv.xz", 
                     "extdata/bioclim_MIROC5_rcp26_2080_landonly.csv.xz", 
                     "extdata/bioclim_MIROC5_rcp60_2080_landonly.csv.xz",
                     "extdata/bioclim_MIROC5_rcp85_2080_landonly.csv.xz")
    
    # Turn data into list
    climatedata <- lapply(climatenames, function(x){
      data <- readr::read_csv(x)
      dplyr::select(data, "x", "y", "bio4", "bio5", 
                    "bio12", "bio15", "bio18", "bio19")
    })
    
    # Predict function depends on mgcv library
    sr_predicted <- lapply(climatedata, function(data){
      pred <- lapply(1:length(m), function(x){
          predict <- round(as.numeric(
            predict(m[[x]]$mod, newdata=data, 
                    type="response",se.fit=FALSE)), 4)
          predict <- as.data.frame(cbind(data[,c("x","y")],predict))
          colnames(predict) <- c("x", "y", x)
          return(predict)
        })
      pred <- Reduce(function(...) merge(..., by=c("x","y"), all.x=T), pred)
      
      # Calculate mean prediction
      pred$mean <- rowMeans(pred[,c(3:12)])
      
      # Merge predicted richness
      pred$model <- model_fam
      pred$taxa <- "bats"
      return(pred)
    })
    
    # Change column names
    sr_predicted <- lapply(1:length(sr_predicted), function(x){
      data <- sr_predicted[[x]] %>% select(x,y,mean,model,taxa)
      data$gcm <- paste0(
        strsplit(climatenames[x], split="_")[[1]][2:4], collapse="_")
      data <- tidyr::spread(data, gcm, mean)
      return(data)
    })
    
    # Merge data
    sr_predicted <- Reduce(function(...) dplyr::left_join(..., by=c("x","y", "model", "taxa"), all.x=TRUE), sr_predicted); rm(m)
    readr::write_csv(sr_predicted, paste0("data/sr_predicted_mem_", y, "_eco_2080.csv.xz"))
  } else{
    NULL
  }
})
sfStop()
