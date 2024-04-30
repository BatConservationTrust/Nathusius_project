## Script name: NP local bat groups and centres record processing
##
## Purpose of script:
##
## Author: Dr. James Waterman
##
## Date Started: 2024-02-20
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
require(snakecase)
require(here)
require(readxl)
require(reshape)

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

# import disparate data sources into lists ----
file_paths <- list.files(path = here::here("data/input/adhoc_records/excel"), pattern = "\\.xlsx", full.names = TRUE)
file_names <-  gsub(pattern = "\\.xlsx$", replacement = "", x = basename(file_paths))
data_list <- lapply(file_paths, readxl::read_excel)

file_paths <- list.files(path = here::here("data/input/adhoc_records/csv"), pattern = "\\.csv", full.names = TRUE)
file_names <-  gsub(pattern = "\\.csv$", replacement = "", x = basename(file_paths))
data_list2 <- lapply(file_paths, readr::read_csv)

np_local <- c(data_list, data_list2)
rm(data_list, data_list2, file_names, file_paths)

# clean data ----
np_local <- purrr::map(
    np_local, ~ .x %>%
        dplyr::select(
            dplyr::any_of(c("date", "location", "grid", "survey", "roost", "count",
                            "lat", "lon", "note", "species", "ring_number"))) %>%
        dplyr::mutate(across(where(is.double), as.character)))

# merge into single df ----
np_local <- reshape::merge_recurse(np_local)

# clean dates ----
unique(sort(np_local$date))

np_local %<>% 
    dplyr::mutate(date1 = lubridate::parse_date_time(date, orders = c("dmy", "ymd"))) %>% 
    dplyr::mutate(date2 = dplyr::case_when(is.na(date1) ~ janitor::excel_numeric_to_date(as.numeric(as.character(date)), date_system = "modern"),
                                           TRUE ~ date1)) %>% 
    dplyr::select(-date, -date1) %>%
    dplyr::rename(date = date2)

np_local %>% dplyr::filter(is.na(date)) %>% View()

# clean counts ----
unique(sort(np_local$count))
np_local %<>% 
    dplyr::rename(count_raw = count) %>% 
    dplyr::mutate(count = stringr::str_replace_all(count_raw, "[:alpha:]", ""),
                  count = stringr::str_replace_all(count, "[:space:]", ""),
                  count = stringr::str_replace_all(count, "[:punct:]", ""),
                  count = stringr::str_extract(count, "[:digit:]+"),
                  count = dplyr::if_else(count == "5515", "55", count),
                  count = as.integer(count))

unique(sort(np_local$count))

np_local %>% dplyr::filter(is.na(count)) %>% View()

# extracting counts from notes may take too long to do pre-deadline (leave for now - focus on locations) ----
# 
# keywords <- c("recorded", "pass", "registration", "call", "heard", "count", "passes", "recording",
#               "recordings", "anabat", "registrations", "files", "adult", "calls")
# 
# df <- np_local %>% 
#     dplyr::mutate(count_flag = dplyr::case_when(is.na(count) & stringr::str_detect(note, paste0(keywords, collapse = "|")) ~ "flag",
#                   TRUE ~ NA_character_))
# 
# df %>% dplyr::filter(is.na(count) & count_flag == "flag") %>% View()
# 
# unique(sort(df$count2))

# convert grid to lat lon ----
np_local$row_id <-  seq.int(nrow(np_local))
# data.table::fwrite(np_local, "np_local.csv")
    
grid_locs <- readr::read_csv("data/input/np_local_grid_conversion.csv")                  

np_local <- dplyr::left_join(np_local, grid_locs, by = c("row_id", "grid"))

unique(sort(np_local$date))

np_local %<>% 
    dplyr::mutate(lat.y = ifelse(is.na(lat.y), lat.x, lat.y),
                  lon.y = ifelse(is.na(lon.y), lon.x, lon.y)) %>% 
    dplyr::rename(lat = lat.y,
                  lon = lon.y) %>% 
    dplyr::select(-lat.x, -lon.x, -row_id) %>% 
    dplyr::mutate(lat = dplyr::case_when(is.na(lat) & grid == "SW53M" ~ "50.165288",
                                         is.na(lat) & grid == "SCC321896" ~ "54.274398",
                                         is.na(lat) & location == "Triton Knoll Wind Farm" ~ "56.451685",
                                         is.na(lat) & location == "Fishtoft" ~ "52.960922",
                                         is.na(lat) & location == "Crooksfoot Reservoir" ~ "54.674064",
                                         TRUE ~ lat),
                  lon = dplyr::case_when(is.na(lon) & grid == "SW53M" ~ "-5.4309610",
                                         is.na(lon) & grid == "SCC321896" ~ "-4.5801096",
                                         is.na(lon) & location == "Triton Knoll Wind Farm" ~ "3.563109",
                                         is.na(lon) & location == "Fishtoft" ~ "0.030541",
                                         is.na(lon) & location == "Crooksfoot Reservoir" ~ "-1.332074",
                                         TRUE ~ lon),
                  year = lubridate::year(date),
                  mutate(across(c(lat, lon), as.numeric)))

np_local %<>% dplyr::distinct()

# grid = SW53M (50.165288 , -5.4309610)
# Triton Knoll Wind Farm (56.451685, 3.563109)
# Fishtoft (52.960922, 0.030541)
# Crooksfoot Reservoir (54.674064, -1.332074)
# grid = SCC321896 (54.274398 , -4.5801096)

# rough plot of local records data ----

# save cleaned NP UK adhoc data ----
saveRDS(np_local, "~/MEGA/Bat Conservation Trust/NNPP/rds/np_uk_adhoc_data_01.rds")

# extract any ringed bat data ----
np_uk_adhoc_ringing_data <- np_local %>% 
  dplyr::filter(!is.na(ring_number)) %>% 
  dplyr::distinct(ring_number, date, .keep_all = TRUE)

# add post-draft extra data ----
np_uk_adhoc_ringing_data <- np_uk_adhoc_ringing_data %>% 
  tibble::add_row(ring_number = c("J6163", "J6163"),
                  location = c("Long pits, Dungeness", "Bedfont Lakes Country Park"),
                  note = c("Added late, some info. absent", "Added late, some info. absent"),
                  date = c(as.Date("2016-09-30"), as.Date("2023-09-10")),
                  lat = c(50.930782, 51.444133),
                  lon = c(0.95837119, -0.44841841),
                  year = c(2016, 2023))

# save cleaned NP UK adhoc ringing data ----
saveRDS(np_uk_adhoc_ringing_data, "~/MEGA/Bat Conservation Trust/NNPP/rds/np_uk_adhoc_ringing_data.rds")


