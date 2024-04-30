## Script name: NP International and all-species UK Ringing data cleaning 01
##
## Purpose of script: Clean NP International Ringing data
##
## Author: Dr. James Waterman
##
## Date Started: 2024-02-17
##
## Copyright (c) James Waterman, 2024
## Email: james.o.waterman@gmail.com
##
## Notes:
##   
## packages ----
require(tidyverse)
require(magrittr)
require(data.table)
require(glmmTMB)
require(lme4)
require(DHARMa)
require(hrbrthemes)
require(janitor)
require(effects)
require(car)
require(ggmap)
require(snakecase)

## options ----
options(scipen = 6, digits = 8, max.print = 20000)
ggplot2::theme_set(theme_light())
rstudioapi::writeRStudioPreference("data_viewer_max_columns", 1000L)
extrafont::loadfonts()
# pdf(".pdf", width = 11.7, height = 8.3, paper = "a4r", bg = "white")

## functions ----
diagnostics.func <- function(x){
    par(mfrow=c(2,2))
    simres <- DHARMa::simulateResiduals(fittedModel=x , n=250)
    DHARMa::plotQQunif(simres)
    DHARMa::plotResiduals(simres, rank=T, quantreg=T)
    DHARMa::testZeroInflation(simres)
    DHARMa::testDispersion(simres)
}

# define a function to fix the bbox to be in EPSG:3857 (not needed at present)
ggmap_bbox <- function(map) {
    if (!inherits(map, "ggmap")) stop("map must be a ggmap object")
    # Extract the bounding box (in lat/lon) from the ggmap to a numeric vector, 
    # and set the names to what sf::st_bbox expects:
    map_bbox <- setNames(unlist(attr(map, "bb")), 
                         c("ymin", "xmin", "ymax", "xmax"))
    
    # Coonvert the bbox to an sf polygon, transform it to 3857, 
    # and convert back to a bbox (convoluted, but it works)
    bbox_3857 <- st_bbox(st_transform(st_as_sfc(st_bbox(map_bbox, crs = 4326)), 3857))
    
    # Overwrite the bbox of the ggmap object with the transformed coordinates 
    attr(map, "bb")$ll.lat <- bbox_3857["ymin"]
    attr(map, "bb")$ll.lon <- bbox_3857["xmin"]
    attr(map, "bb")$ur.lat <- bbox_3857["ymax"]
    attr(map, "bb")$ur.lon <- bbox_3857["xmax"]
    map
}

# import international ringing data ----
ring_int <- readr::read_csv("~/MEGA/Bat Conservation Trust/Nathusius' Pipistrelle Project/data/input/international_ring_recaptures.csv")
dplyr::glimpse(ring_int)

# tidy international ringing data (variable names; date and datetime format) ----
ring_int <- janitor::clean_names(ring_int)

ring_int %<>%
    dplyr::rename(capture_lat = x_start,
                  capture_lon = y_start,
                  note = summary) %>%
    dplyr::mutate(across(c(ringing_date, recapture_date),
                         ~ lubridate::parse_date_time(., orders = c("dmy"))),
                  species = "Nathusisus' pipistrelle",
                  note = snakecase::to_mixed_case(note),
                  ring_number = stringr::word(ring_number, -1)) %>%
    dplyr::select(species:ringing_date, capture_lat, capture_lon, x_end:note)


# pivot long (tidy data are 1 row = 1 observation) ----
df <- ring_int %>%
    dplyr::select(species, sex, ring_number, x_end, y_end, recapture_date, note) %>%
    dplyr::mutate(capture_instance = 2) %>%
    dplyr::rename(capture_lat = x_end,
                  capture_lon = y_end,
                  capture_date = recapture_date)

ring_int %<>%
    dplyr::select(-x_end, -y_end, -recapture_date) %>%
    dplyr::mutate(capture_instance = 1) %>%
    dplyr::rename(capture_date = ringing_date)

ring_int_long <- dplyr::bind_rows(ring_int, df) %>% dplyr::mutate(region = "international")

# EDA plot of international capture-recapture locations ----
# top = northerly lat, bottom = southerly lat, left = westerly lon, right = easterly lon

# ring_int_bbox <- c(top = 59.349107, bottom = 48.321177, left =  -13.765270, right = 38.934531)
# int_ring_map_terrain <- ggmap::get_map(location = ring_int_bbox, maptype = "stamen_terrain", source = "stadia",
#                                        zoom = 6, scale = 1)

int_ring_map_terrain <- readRDS("~/MEGA/Bat Conservation Trust/Nathusius' Pipistrelle Project/rds/int_ring_map_terrain.rds")

NP_international_ringing_map <- 
    ggmap::ggmap(int_ring_map_terrain)+
    geom_point(data = ring_int_long, mapping = aes(x = capture_lat, y = capture_lon, color = ring_number), alpha = 0.9)+
    geom_line(data = ring_int_long, mapping = aes(x = capture_lat, y = capture_lon, color = ring_number, group = ring_number))+
    labs(title = "Nathusius' pipistrelle International Ringing Data",
         x = "Latitude", y = "Longitude",
         colour = "Ring ID")+
    theme_ipsum()

# ggsave(filename = "NP_international_ringing_map.svg",
#        plot = NP_international_ringing_map,
#        width = 29.7, height = 21.0, units = 'cm',
#        scale = 1, dpi = 600,
#        bg = "white",
#        limitsize = FALSE)
# 
# ggsave(filename = "NP_international_ringing_map.png",
#        plot = NP_international_ringing_map,
#        width = 29.7, height = 21.0, units = 'cm',
#        scale = 1, dpi = 600,
#        bg = "white",
#        limitsize = FALSE)
# 
# ggsave(filename = "NP_international_ringing_map.pdf",
#        plot = NP_international_ringing_map,
#        width = 29.7, height = 21.0, units = 'cm',
#        scale = 1, dpi = 600,
#        bg = "white",
#        limitsize = FALSE)

# save international ringing basemap ----
saveRDS(int_ring_map_terrain, "~/MEGA/Bat Conservation Trust/Nathusius' Pipistrelle Project/rds/int_ring_map_terrain.rds")

# save cleaned NP international ringing data ----
saveRDS(ring_int_long, "~/MEGA/Bat Conservation Trust/Nathusius' Pipistrelle Project/rds/np_internatinal_ringing_01.rds")

# # import clean NNPP survey data (to extract intra-UK ringing recoveries) ----
nnpp_clean_02 <- readRDS("~/MEGA/Bat Conservation Trust/Nathusius' Pipistrelle Project/rds/nnpp_clean_02.rds")

# # extract intra-UK ringing records from NNPP trapping surveys ----
# ring_uk <- nnpp_clean_02 %>% 
#     dplyr::group_by(ring_number) %>% 
#     dplyr::filter(!is.na(ring_number) & n() > 1)
# 
# janitor::tabyl(ring_uk, species)
# janitor::tabyl(ring_uk, ring_number, species)

# save cleaned all-species UK ringing data ----
# saveRDS(ring_uk, "~/MEGA/Bat Conservation Trust/Nathusius' Pipistrelle Project/rds/uk_ringing_01.rds")

# EDA plot of UK capture-recapture locations ----
# top = northerly lat, bottom = southerly lat, left = westerly lon, right = easterly lon

# uk_bbox <- c(top = 61.17, left = -12.33, bottom = 47.10, right = 4.50)
# uk_ring_map_terrain <- ggmap::get_map(location = uk_bbox, maptype = "stamen_terrain", source = "stadia",
#                                        zoom = 8, scale = 1)

# save UK ringing basemap ----
# saveRDS(uk_ring_map_terrain, "~/MEGA/Bat Conservation Trust/Nathusius' Pipistrelle Project/rds/uk_ring_map_terrain.rds")

ring_uk <- readRDS("~/MEGA/Bat Conservation Trust/Nathusius' Pipistrelle Project/rds/uk_ringing_01.rds")
uk_ring_map_terrain <- readRDS("~/MEGA/Bat Conservation Trust/Nathusius' Pipistrelle Project/rds/uk_ring_map_terrain.rds")

NP_UK_ringing_map <- 
    ggmap::ggmap(uk_ring_map_terrain)+
    geom_point(data = subset(ring_uk, species == "Nathusius' pipistrelle"),
               aes(y = trap_lat, x = trap_lon, color = ring_number),
               alpha = 0.9, size = 0.3)+
    geom_line(data = subset(ring_uk, species == "Nathusius' pipistrelle"),
              aes(y = trap_lat, x = trap_lon, color = ring_number, group = ring_number),
              alpha = 0.5, linewidth = 0.3)+
    labs(title = "Nathusius' pipistrelle UK Ringing Data",
         x = "Latitude", y = "Longitude",
         colour = "Ring ID")+
    theme_ipsum()+
    theme(legend.position = "none")

# ggsave(filename = "NP_UK_ringing_map.svg",
#        plot = NP_UK_ringing_map,
#        width = 21, height = 29.7, units = 'cm',
#        scale = 1, dpi = 600,
#        bg = "white",
#        limitsize = FALSE)
# 
# ggsave(filename = "NP_UK_ringing_map.png",
#        plot = NP_UK_ringing_map,
#        width = 21, height = 29.7, units = 'cm',
#        scale = 1, dpi = 600,
#        bg = "white",
#        limitsize = FALSE)
# 
# ggsave(filename = "NP_UK_ringing_map.pdf",
#        plot = NP_UK_ringing_map,
#        width = 21, height = 29.7, units = 'cm',
#        scale = 1, dpi = 600,
#        bg = "white",
#        limitsize = FALSE)
