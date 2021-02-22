#' ---
#' title: "Future global bat richness using MEMs"
#' author: "RS-eco"
#' ---

#' This is a first analysis based on our MEM models for the 1265 bat species.
#' The MEMs were run with 4 bioclimatic variables: bio4, bio5, bio12 and bio15

#+ setup, include=FALSE
knitr::opts_chunk$set(collapse = TRUE, warning=FALSE, echo=FALSE, message=FALSE, fig.width=8, fig.height=6)

#' ## Load and prepare data

rm(list=ls()); invisible(gc())

# Set working directory
filedir <- "/home/matt/Documents/futBats/data/"

# Load packages
library(tidyverse); library(patchwork); library(magrittr); library(sf)
library(scico)

# Specify colour scheme
bluered <- scico(n=255, palette="roma")
bluewhitered <-  scico(n=255, palette="romaO")

# Obtain world map
outline <- rgeos::gPolygonize(rgeos::gNode(as(rworldmap::getMap(resolution = "high"), "SpatialLines")))
outline <- rgeos::gUnaryUnion(outline)
outline <- sf::st_as_sf(outline)

# Load MEM models

m <- get(load(paste0(filedir, "bio4_bio5_bio12_bio15_model_output_GAM_poisson_Eco_block.RData")))

names(m[[1]])

#' ## Analyse model accuracy using AUC



#' ## Analyse variable importance

# Try varImp on GBM
#varImp

#' ## Analyse variable response

#library(caret)
#+ echo=F
#par(mfrow=c(2,2))
#varResp <- lapply(m, function(x) plot(x$mod, se=F))
#varResp <- lapply(varResp, function(y){
#  dat <- lapply(y, function(z){
#    data.frame(var=z$xlab, value=z$x, fit=z$fit)
#  })
#  bind_rows(dat)
#})
#varResp <- bind_rows(varResp)

#+ echo=T
#varResp %>% mutate(value=round(value, 0)) %>%
#  group_by(var, value) %>% summarise(fit=mean(fit)) %>%
#  ggplot() + geom_path(aes(x=value, y=fit)) + facet_wrap(.~var, scales="free")

# Load MEM data

#list.files(pattern=".csv.xz")

bat_gam_1995 <- read.csv(paste0(filedir, "sr_predicted_mem_GAM_poisson_eco_1995.csv.xz"))
bat_gam_1995$model <- "GAM"
#bat_gam_1995$year <- 1995
colnames(bat_gam_1995) <- c("x", "y", paste0("block", 1:10), "EWEMBI_1995", "model", "taxa")
bat_gam_2080 <- read.csv(paste0(filedir, "sr_predicted_mem_GAM_poisson_eco_2080.csv.xz"))
bat_gam_2080$model <- "GAM"
#bat_gam_2080$year <- 2080
bat_gam <- left_join(bat_gam_1995, bat_gam_2080); rm(bat_gam_1995, bat_gam_2080)

bat_gbm_1995 <- read.csv(paste0(filedir, "sr_predicted_mem_GBM_poisson_eco_1995.csv.xz"))
bat_gbm_1995$model <- "GBM"
#bat_gbm_1995$year <- 1995
colnames(bat_gbm_1995) <- c("x", "y", paste0("block", 1:10), "EWEMBI_1995", "model", "taxa")
bat_gbm_2080 <- read.csv(paste0(filedir, "sr_predicted_mem_GBM_poisson_eco_2080.csv.xz"))
bat_gbm_2080$model <- "GBM"
#bat_gbm_2080$year <- 2080
bat_gbm <- left_join(bat_gbm_1995, bat_gbm_2080); rm(bat_gbm_1995, bat_gbm_2080)

bat_dat <- bind_rows(bat_gam, bat_gbm); rm(bat_gam, bat_gbm)

#' ## Look at richness patterns of MEMs

bat_dat %<>% select(x,y,model, EWEMBI_1995, GFDL.ESM2M_rcp26_2080:MIROC5_rcp85_2080) %>% gather(var, value, -c(x, y, model)) %>% 
  separate(var, into=c("gcm", "scenario", "year"), sep="_", fill="left")
bat_dat %<>% filter(year %in% c(1995, 2080)) %>% group_by(x,y,scenario,year) %>%
  summarise(mn=mean(value, na.rm=T)) %>% select(-c(year)) %>% spread(scenario, mn)
colnames(bat_dat) <- c("x", "y", "EWEMBI 1995", "2080 RCP2.6", "2080 RCP6.0", "2080 RCP8.5")

#' ## Plot richness maps
#+ fig.width=12, fig.height=6
bat_dat %>% gather(scen, value, -c(x,y)) %>% 
  mutate(scen=factor(scen, levels=c("EWEMBI 1995", "2080 RCP2.6", "2080 RCP6.0", "2080 RCP8.5"))) %>% 
  ggplot() + geom_tile(aes(x=x,y=y,fill=value)) +
  facet_wrap(.~scen) + scale_fill_gradientn(name="SR", colours = bluered, na.value="transparent") + 
  geom_sf(data=outline, fill=NA, color = "black") +
  coord_sf(xlim=c(-180,180), ylim=c(-56,84), expand=FALSE, ndiscr = FALSE) + theme_minimal() + 
  theme(axis.title=element_blank(), strip.text = element_text(face="bold", size=9))
ggsave("figures/bat_richness_maps_mem.png", width=12, height=6, dpi=600)

#' **Fig. 1.** Richness maps for 1995 and 2080 under three RCPs (RCP2.6, RCP6.0 and RCP8.5).

#' ## Plot maps of change in richness
#+ fig.width=8, fig.height=10
delta_bat <- bat_dat %>% ungroup() %>% mutate_at(vars(-c(x,y,`EWEMBI 1995`)), `-`, y=.$`EWEMBI 1995`) %>% 
  gather(scen, value, -c(x,y,`EWEMBI 1995`))
col_val <- scales::rescale(c(seq(min(delta_bat$value, na.rm=T), 0, 5), seq(0, max(delta_bat$value, na.rm=T), 5)), na.rm=T)
lim <- c(min(delta_bat$value, na.rm=T),  max(delta_bat$value, na.rm=T))
delta_bat %>% ggplot() + geom_tile(aes(x=x,y=y,fill=value)) +
  facet_wrap(.~scen, ncol=1) + 
  scale_fill_gradientn(name="SR", colours = bluewhitered, na.value="transparent",
                       values=col_val, limits=lim) + 
  geom_sf(data=outline, fill=NA, color = "black") +
  coord_sf(xlim=c(-180,180), ylim=c(-56,84), expand=FALSE, ndiscr = FALSE) + theme_minimal() + 
  theme(axis.title=element_blank(), strip.text = element_text(face="bold", size=9))
ggsave("figures/delta_sr_bats_mem.png", width=8, height=10, dpi=600)

#' **Fig. 2.** Maps of the change in species richness for 2080 compared to 1995 under three RCPs (RCP2.6, RCP6.0 and RCP8.5).

#' ## Plot histograms of change in richness
#+ fig.width=8, fig.height=6
delta_bat %>% ggplot() + geom_histogram(aes(x=value), bins=50) +
  facet_wrap(.~scen, nrow=1) + theme_bw() + 
  labs(x="Change in richness", y="Number of cells") + 
  scale_y_continuous(limits=c(0, NA), expand=c(0,200)) + 
  theme(strip.text = element_text(face="bold", size=9),
        strip.background = element_blank())
ggsave("figures/hist_change_sr_bats_mem.png", width=8, height=6, dpi=600)

#' **Fig. 3.** Histograms of the change in species richness for 2080 compared to 1995 under three RCPs (RCP2.6, RCP6.0 and RCP8.5).

#+ fig.width=5, fig.height=8
delta_bat %>% group_by(y, scen) %>% summarise(mn=mean(value)) %>% 
  ggplot() + geom_path(aes(x=mn, y=y, col=scen)) +theme_bw() + 
  labs(x="Change in richness", y="Latitude") + 
  theme(strip.text = element_text(face="bold", size=9),
        strip.background = element_blank())
ggsave("figures/lat_change_sr_bats_mem.png", width=5, height=8, dpi=600)

#' **Fig. 4.** Latitudinal changes in species richness for 2080 compared to 1995 under three RCPs (RCP2.6, RCP6.0 and RCP8.5).
