# ## Load and prepare data

rm(list=ls()); invisible(gc())

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

# Load SDM data

#list.files(pattern=".csv.xz")

bat_gam <- read.csv("data/bat_prob_GAM_dispersal2.csv.xz")
bat_gam$model <- "GAM"
bat_gbm <- read.csv("data/bat_prob_GBM_dispersal2.csv.xz")
bat_gbm$model <- "GBM"

bat_dat <- bind_rows(bat_gam, bat_gbm); rm(bat_gam, bat_gbm); gc()

## Look at richness patterns of SDMs

bat_dat %<>% gather(var, value, -c(x, y, model)) %>% 
  separate(var, into=c("gcm", "scenario", "year"), sep="_", fill="left")
bat_dat %<>% filter(year %in% c(1995, 2080)) %>% group_by(x,y,scenario,year) %>%
  summarise(mn=mean(value, na.rm=T)) %>% select(-year) %>% spread(scenario, mn)
colnames(bat_dat) <- c("x", "y", "EWEMBI 1995", "2080 RCP2.6", "2080 RCP6.0", "2080 RCP8.5")
gc()

#' ## Plot richness maps
bat_dat %>% gather(scen, value, -c(x,y)) %>% 
  mutate(scen=factor(scen, levels=c("EWEMBI 1995", "2080 RCP2.6", "2080 RCP6.0", "2080 RCP8.5"))) %>% 
  ggplot() + geom_tile(aes(x=x,y=y,fill=value)) +
  facet_wrap(.~scen) + scale_fill_gradientn(name="SR", colours = bluered, na.value="transparent") + 
  geom_sf(data=outline, fill=NA, color = "black") +
  coord_sf(xlim=c(-180,180), ylim=c(-56,84), expand=FALSE, ndiscr = FALSE) + theme_minimal() + 
  theme(axis.title=element_blank(), strip.text = element_text(face="bold", size=9))
ggsave("figures/bat_richness_maps_sdm.png", width=12, height=6, dpi=600)

delta_bat <- bat_dat %>% ungroup() %>% mutate_at(vars(-c(x,y,`EWEMBI 1995`)), `-`, y=.$`EWEMBI 1995`) %>% 
  gather(scen, value, -c(x,y,`EWEMBI 1995`)); rm(bat_dat); gc()
col_val <- scales::rescale(c(seq(min(delta_bat$value, na.rm=T), 0, 5), seq(0, max(delta_bat$value, na.rm=T), 5)), na.rm=T)
lim <- c(min(delta_bat$value, na.rm=T),  max(delta_bat$value, na.rm=T))
delta_bat %>% ggplot() + geom_tile(aes(x=x,y=y,fill=value)) +
  facet_wrap(.~scen, ncol=1) + 
  scale_fill_gradientn(name="SR", colours = bluewhitered, na.value="transparent",
                       values=col_val, limits=lim) + 
  geom_sf(data=outline, fill=NA, color = "black") +
  coord_sf(xlim=c(-180,180), ylim=c(-56,84), expand=FALSE, ndiscr = FALSE) + theme_minimal() + 
  theme(axis.title=element_blank(), strip.text = element_text(face="bold", size=9))
ggsave("figures/delta_sr_bats_sdm.png", width=8, height=10, dpi=600)

delta_bat %>% ggplot() + geom_histogram(aes(x=value), bins=50) +
  facet_wrap(.~scen, nrow=1) + theme_bw() + 
  labs(x="Change in richness", y="Number of cells") + 
  scale_y_continuous(limits=c(0, NA), expand=c(0,200)) + 
  theme(strip.text = element_text(face="bold", size=9),
        strip.background = element_blank())
ggsave("figures/hist_change_sr_bats_sdm.png", width=8, height=6, dpi=600)

delta_bat %>% group_by(y, scen) %>% summarise(mn=mean(value)) %>% 
  ggplot() + geom_path(aes(x=mn, y=y, col=scen)) +theme_bw() + 
  labs(x="Change in richness", y="Latitude") + 
  theme(strip.text = element_text(face="bold", size=9),
        strip.background = element_blank())
ggsave("figures/lat_change_sr_bats_sdm.png", width=5, height=8, dpi=600)

