## Script name: NNPP data cleaning 01
##
## Purpose of script: Clean and combine disparate NP data sources
##
## Author: Dr. James Waterman
##
## Date Started: 2024-02-13
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
require(rnrfa)
require(anytime)
library(kde1d)

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

format_cells <- function(df, rows ,cols, value = c("italics", "bold", "strikethrough")){
  
  # select the correct markup
  # one * for italics, two ** for bold
  map <- setNames(c("*", "**", "~~"), c("italics", "bold", "strikethrough"))
  markup <- map[value]  
  
  for (r in rows){
    for(c in cols){
      
      # Make sure values are not factors
      df[[c]] <- as.character( df[[c]])
      
      # Update formatting
      df[r, c] <- paste0(markup, df[r, c], markup)
    }
  }
  
  return(df)
}
# import NNPP capture data ----
nnpp <- readr::read_csv("~/MEGA/Bat Conservation Trust/NNPP/data/input/NNPP complete database export_090224.csv")
dplyr::glimpse(nnpp)

# tidy NNPP capture data (variable names; date and datetime format) ----
nnpp <- janitor::clean_names(nnpp)

# n 18269 ----

# identify KEY fields with missing data ----

# bat_group_id = all present
# bat_group_name = all present
# site_id = absent * 1
# site_name = absent * 1 (same as site_id)
# survey_id = absent * 6 (1 same site_id & site_name)
# survey_date_start = absent * 6 (same as survey_id); 22 records have time but no date || use date from bat_time
# survey_date_end = absent * 6 (same as survey_date_start); 22 records have time but no date || use date from bat_time
# start_temp range(0 - 27); absent = 3158 of 18269 (17.29%)
# end_temp range(0 - 23); absent = 5031 of 18269 (27.54%)
# both start and end temp absent = 3052
# bat_record_id - absent * 79
# bat_number = absent * 129; 1 value = zero (error?)
# bat_time = absent * 79
# species = absent * 79; species = "bat_sp" * 6 (presumably unidentified at species level)
# sex = absent * 123; sex = unknown * 210; 1 value = zero (error?)
# age = absent * 213; age = unknown * 1508; 99 values = zero (error?)
# female_reproductive_status = absent * 207 of 9222 female records; unknown * 2563 of 9222 female records; 132 values = zero (error?)
# male_testes_size = absent * 1526 of 8714 male records; unknown * 1634 of 8714 male records; 49 values = zero (error?)
# male_epididymis_colour = absent * 152 of 8714 male records; unknown * 3792 of 8714 male records; 51 values = zero (error?) 
# bat_weight (range 2.1 - 43.0 excluding oddities) = absent * 4945 of 18269 (27.1%); 1 value = 577 (typo?); 62 values = 0 (error?)
# forearm (range 21.5 - 57.2 excluding oddities) = absent * 4803 of 18269 (26.29%); 1 value = 317 (typo?); 73 values = 0 (error?) mostly same as bat_weight
# ring_number = absent * 14918 of 18269 (81.66%); 595 values = 0 (error?)
# trap_id = absent * 712
# trap_gps = absent * 1064 of 18269 (5.82%); unknown * 260
# equipment_used = absent * 3150; unknown * 651; 8 values = 0 (error?)

# mark records with potentially problematic missing data ----
# continue cleaning while investigating and querying these data with BCT staff
nnpp %<>% 
  dplyr::mutate(absent_key_data = dplyr::case_when(
    is.na(site_id) ~ TRUE, 
    is.na(site_name) ~ TRUE,
    is.na(survey_id) ~ TRUE,
    is.na(survey_date_start) ~ TRUE,
    is.na(survey_date_end) ~ TRUE,
    is.na(species) ~ TRUE,
    is.na(trap_gps) ~ TRUE, TRUE ~ FALSE
  ))

# nnpp %>% dplyr::filter(absent_key_data == TRUE) %>% View()

# extract problematic records for examination (empty bat_record_id are ZERO bat surveys - valuable) ----
nnpp_absences <- nnpp %>% dplyr::filter(absent_key_data == TRUE & !is.na(bat_record_id))

# save and export problematic records for examination ----
saveRDS(nnpp_absences, "~/MEGA/Bat Conservation Trust/NNPP/rds/nnpp_absences.rds")

# extract survey_ids where ZERO bats were trapped (according to PB should - 78) -----
# use this vector/df of survey_id later for calculation/correction of capture rate
nnpp_zero_bat_surveys_confirmed <- nnpp %>% dplyr::filter(is.na(bat_record_id))

# proceed without these "problematic record" data for cleaning framework ----
# nnpp1 <- nnpp %>% dplyr::filter(absent_key_data == FALSE)

# bind zero bat survey records to main df ----
# nnpp1 <- dplyr::bind_rows(nnpp1, nnpp_zero_bat_surveys_confirmed)
nnpp %>% dplyr::filter(is.na(bat_record_id))
nnpp %<>% dplyr::mutate(zero_bat_survey = dplyr::case_when(is.na(bat_record_id) ~ TRUE, TRUE ~ FALSE))

# clean var names, date and datetime formats ----
nnpp %<>%  
  dplyr::rename(start_temp_C = start_temp,
                end_temp_C = end_temp,
                bat_mass_g = bat_weight,
                forearm_length_cm = forearm,
                fifth_finger_forearm_ratio = fifth_finger) %>%
  tidyr::separate(survey_date_start, into = c("survey_date_start", "survey_time_start"), sep = " ", remove = T) %>% 
  tidyr::separate(survey_date_end, into = c("survey_date_end", "survey_time_end"), sep = " ", remove = T) %>%
  tidyr::separate(bat_time, into = c("bat_trap_date", "bat_trap_time"), sep = " ", remove = T) %>% 
  dplyr::mutate(across(c(survey_date_start, survey_date_end, bat_trap_date),
                       ~ lubridate::parse_date_time(., orders = c("ymd", "dmy")))) %>% 
  dplyr::mutate(year = lubridate::year(survey_date_start),
                month = lubridate::month(survey_date_start),
                day = lubridate::day(survey_date_start))

# create survey start, end, and bat time datetime objects (time zone = UTC avoids issues with daylight savings) ----
nnpp %<>% 
  dplyr::mutate(survey_datetime_start = as.POSIXct(paste(survey_date_start, survey_time_start, sep = " "), tz = "UTC", format = "%Y-%m-%d %H:%M"),
                survey_datetime_end = as.POSIXct(paste(survey_date_end, survey_time_end, sep = " "), tz = "UTC", format = "%Y-%m-%d %H:%M"),
                bat_trap_datetime = as.POSIXct(paste(bat_trap_date, bat_trap_time, sep = " "), tz = "UTC", format = "%Y-%m-%d %H:%M"))%>% 
  dplyr::select(bat_group_id:survey_id, year, month, day, survey_date_start, survey_time_start,
                survey_datetime_start, survey_date_end, survey_time_end, survey_datetime_end,
                bat_trap_date, bat_trap_time, bat_trap_datetime, bat_record_id, bat_number,
                species, ring_number, sex:fifth_finger_forearm_ratio, trap_id:trap_gps, everything())


# convert Ordnance Survey (OS) grid reference to latitude/longitude coordinates ----
# clean jumbled entries pre-conversion (eliminate spaces, deal with odd NA/not given, unknown entries, extract existing lat/lon)
unique(sort(nnpp$trap_gps))
nnpp %<>% 
  dplyr::mutate(trap_gps = stringr::str_replace_all(trap_gps, " ", "")) %>% 
  tidyr::separate(trap_gps, into = c("lon1", "lat1"), sep = ",", remove = FALSE) %>% 
  dplyr::mutate(trap_gps = dplyr::case_when(trap_gps %in% c("N/A", "unknown", "Unknown", "Notgiven") ~ NA, TRUE ~ trap_gps),
                trap_gps = ifelse(stringr::str_detect(trap_gps, ","), NA, trap_gps),
                lon1 = ifelse(stringr::str_detect(lon1, "[:alpha:]"), NA, lon1)) %>% 
  dplyr::rename(trap_osg = trap_gps)

unique(sort(nnpp$trap_osg))
unique(sort(nnpp$lon1))
unique(sort(nnpp$lat1))

nnpp %<>% 
  dplyr::group_by(site_name) %>% 
  tidyr::fill(., trap_osg, .direction = "downup")

# create random survey_id numbers for ZERO BAT surveys ----
max(nnpp$survey_id, na.rm = T)
sum(is.na(nnpp$survey_id))
jw_survey_id <- 1970:1975
nnpp$survey_id[is.na(nnpp$survey_id)] <- jw_survey_id[1:sum(is.na(nnpp$survey_id))]
max(nnpp$survey_id, na.rm = T)
sum(is.na(nnpp$survey_id))
rm(jw_survey_id)

# data.table::fwrite(nnpp, "nnpp.csv")

# import transformed trap coordinates (internal function failing) ----
# https://gridreferencefinder.com/batchConvert/batchConvert.php
nnpp_locs <- readr::read_csv("~/MEGA/Bat Conservation Trust/NNPP/data/input/nnpp_trap_locations.csv")
dplyr::glimpse(nnpp_locs)

nnpp_locs %<>% dplyr::distinct()

nnpp <- dplyr::left_join(nnpp, dplyr::select(nnpp_locs, -trap_osg), by = c("survey_id"))

# n 49421 ----

nnpp %>% dplyr::filter(is.na(trap_lat)) %>% View()

# clean inconsistencies in factor levels ----
unique(sort(nnpp$bat_group_name))
unique(sort(nnpp$site_name))
unique(sort(nnpp$survey_id))
unique(sort(nnpp$bat_record_id))
unique(sort(nnpp$bat_number))

# rename unknown bat species ----
unique(sort(nnpp$species))

nnpp %<>% dplyr::mutate(species = dplyr::case_when(species == "Bat sp." ~ "Unknown", TRUE ~ species))
unique(sort(nnpp$species))

unique(sort(nnpp$ring_number))
# scrub 'odd' ring numbers (H8841 & H8841Recap), (7984 & J7984 same animal), and multiple zeros ----
nnpp %<>% 
  dplyr::mutate(ring_number = dplyr::case_when(ring_number == "H8841Recap" ~ "H8841",
                                               ring_number == "0" ~ NA_character_,
                                               ring_number == "7984" ~ "J7984",
                                               TRUE ~ ring_number))
unique(sort(nnpp$ring_number))

# scrub sex ----
unique(sort(nnpp$sex))
length(nnpp$zero_bat_survey[nnpp$zero_bat_survey == TRUE])

nnpp %<>% 
  dplyr::mutate(sex = dplyr::case_when(
    sex == "0" | sex == "Unknown" | is.na(sex) ~ "Unknown",
    sex == "F" | sex == "female" | sex == "Female" ~ "Female",
    sex == "M" | sex == "Male" ~ "Male",
    TRUE ~ sex))

unique(sort(nnpp$sex))

nnpp %>%
  ggplot(aes(x = sex)) +
  geom_histogram(fill = "#404080", colour = "#e9ecef", alpha = 0.6, position = 'identity', stat = "count") +
  theme_ipsum()

nnpp %>%
  ggplot(aes(x = sex)) +
  geom_histogram(fill = "#404080", colour = "#e9ecef", alpha = 0.6, position = 'identity', stat = "count") +
  theme_ipsum()+
  facet_wrap(~species)+
  theme(axis.text.x=element_text(angle = 45, hjust = 1))

# scrub age ----
unique(sort(nnpp$age))
length(nnpp$zero_bat_survey[nnpp$zero_bat_survey == TRUE])

nnpp %<>% 
  dplyr::mutate(age = dplyr::case_when(
    age == "0" ~ "Unknown",
    age == "A" | age == "adult" ~ "Adult",
    age == "J" ~ "Juvenile",
    TRUE ~ age))

unique(sort(nnpp$age))

nnpp %>%
  ggplot(aes(x = age)) +
  geom_histogram(fill = "#404080", colour = "#e9ecef", alpha = 0.6, position = 'identity', stat = "count") +
  theme_ipsum()

nnpp %>%
  ggplot(aes(x = age)) +
  geom_histogram(fill = "#404080", colour = "#e9ecef", alpha = 0.6, position = 'identity', stat = "count") +
  theme_ipsum()+
  facet_wrap(~species)+
  theme(axis.text.x=element_text(angle = 45, hjust = 1))

# scrub male testes size ----
unique(sort(nnpp$male_testes_size))

nnpp %<>% 
  dplyr::mutate(male_testes_size = dplyr::case_when(
    male_testes_size == "0" | male_testes_size == "Unknown" | male_testes_size == "None" | male_testes_size == "unknown" | male_testes_size == "\"\"" ~ "Unknown",
    male_testes_size == "large" | male_testes_size == "Large" ~ "Large",
    male_testes_size == "medium" | male_testes_size == "Medium" ~ "Medium",
    male_testes_size == "small" | male_testes_size == "Small" ~ "Small",
    TRUE ~ male_testes_size))

unique(sort(nnpp$male_testes_size))

nnpp %>%
  ggplot(aes(x = male_testes_size)) +
  geom_histogram(fill = "#404080", colour = "#e9ecef", alpha = 0.6, position = 'identity', stat = "count") +
  theme_ipsum()

nnpp %>%
  ggplot(aes(x = male_testes_size)) +
  geom_histogram(fill = "#404080", colour = "#e9ecef", alpha = 0.6, position = 'identity', stat = "count") +
  theme_ipsum()+
  facet_wrap(~species)+
  theme(axis.text.x=element_text(angle = 45, hjust = 1))

# scrub male epididymis colour ----
unique(sort(nnpp$male_epididymis_colour))

nnpp %<>% 
  dplyr::mutate(male_epididymis_colour = dplyr::case_when(
    male_epididymis_colour == "dark" ~ "Dark",
    male_epididymis_colour == "pale" ~ "Pale",
    male_epididymis_colour == "0" | male_epididymis_colour == "\"\"" ~ "Unknown",
    male_epididymis_colour == "black epi" ~ "Black",
    male_epididymis_colour == "patchy epi" ~ "Patchy",
    TRUE ~ male_epididymis_colour))

unique(sort(nnpp$male_epididymis_colour))

nnpp %>%
  ggplot(aes(x = male_epididymis_colour)) +
  geom_histogram(fill = "#404080", colour = "#e9ecef", alpha = 0.6, position = 'identity', stat = "count") +
  theme_ipsum()

nnpp %>%
  ggplot(aes(x = male_epididymis_colour)) +
  geom_histogram(fill = "#404080", colour = "#e9ecef", alpha = 0.6, position = 'identity', stat = "count") +
  theme_ipsum()+
  facet_wrap(~species)+
  theme(axis.text.x=element_text(angle = 45, hjust = 1))

# scrub female reproductive status ----
unique(sort(nnpp$female_reproductive_status))

nnpp %<>% 
  dplyr::mutate(female_reproductive_status = dplyr::case_when(
    female_reproductive_status == "0" | female_reproductive_status == "Unknown" ~ "Unknown",
    female_reproductive_status == "Non parous" | female_reproductive_status == "NonParous" | female_reproductive_status == "nulliparous" ~ "Non-Parous",
    female_reproductive_status == "parous" ~ "Parous",
    female_reproductive_status == "Post lactating" ~ "Post-Lactating",
    TRUE ~ female_reproductive_status))

unique(sort(nnpp$female_reproductive_status))

nnpp %>%
  ggplot(aes(x = female_reproductive_status)) +
  geom_histogram(fill = "#404080", colour = "#e9ecef", alpha = 0.6, position = 'identity', stat = "count") +
  theme_ipsum()

nnpp %>%
  ggplot(aes(x = female_reproductive_status)) +
  geom_histogram(fill = "#404080", colour = "#e9ecef", alpha = 0.6, position = 'identity', stat = "count") +
  theme_ipsum()+
  facet_wrap(~species)+
  theme(axis.text.x=element_text(angle = 45, hjust = 1))

# scrub bat mass ----
unique(sort(nnpp$bat_mass_g))

nnpp %>%
  ggplot(aes(x = bat_mass_g)) +
  geom_histogram(fill = "#404080", colour = "#e9ecef", alpha = 0.6, position = 'identity') +
  theme_ipsum()

nnpp %<>% 
  dplyr::mutate(bat_mass_g = dplyr::case_when(bat_mass_g == 577 ~ 5.7,
                                              bat_mass_g == 0 ~ NA,
                                              TRUE ~ bat_mass_g))
unique(sort(nnpp$bat_mass_g))

nnpp %>%
  ggplot(aes(x = bat_mass_g)) +
  geom_histogram(fill = "#404080", colour = "#e9ecef", alpha = 0.6, position = 'identity') +
  theme_ipsum()

nnpp %>%
  ggplot(aes(x = bat_mass_g)) +
  geom_histogram(fill = "#404080", colour = "#e9ecef", alpha = 0.6, position = 'identity') +
  theme_ipsum()+
  facet_wrap(~species)

# scrub forearm length ----
unique(sort(nnpp$forearm_length_cm))

nnpp %>%
  ggplot(aes(x = forearm_length_cm)) +
  geom_histogram(fill = "#404080", colour = "#e9ecef", alpha = 0.6, position = 'identity') +
  theme_ipsum()

nnpp %<>% 
  dplyr::mutate(forearm_length_cm = dplyr::case_when(forearm_length_cm == 317 ~ 31.7,  # consistent with soprano range
                                                     forearm_length_cm == 0 ~ NA,
                                                     TRUE ~ forearm_length_cm))

unique(sort(nnpp$forearm_length_cm))

nnpp %>%
  ggplot(aes(x = forearm_length_cm)) +
  geom_histogram(fill = "#404080", colour = "#e9ecef", alpha = 0.6, position = 'identity') +
  theme_ipsum()

nnpp %>%
  ggplot(aes(x = forearm_length_cm)) +
  geom_histogram(fill = "#404080", colour = "#e9ecef", alpha = 0.6, position = 'identity') +
  theme_ipsum()+
  facet_wrap(~species)

# scrub fifth_finger_forearm_ratio (NP only) ----
unique(sort(nnpp$fifth_finger_forearm_ratio))

nnpp %>%
  ggplot(aes(x = fifth_finger_forearm_ratio)) +
  geom_histogram(fill = "#404080", colour = "#e9ecef", alpha = 0.6, position = 'identity') +
  theme_ipsum()

nnpp %<>% 
  dplyr::mutate(fifth_finger_forearm_ratio = dplyr::case_when(fifth_finger_forearm_ratio == 0 ~ NA,
                                                              TRUE ~ fifth_finger_forearm_ratio))

unique(sort(nnpp$fifth_finger_forearm_ratio))

nnpp %>%
  ggplot(aes(x = fifth_finger_forearm_ratio)) +
  geom_histogram(fill = "#404080", colour = "#e9ecef", alpha = 0.6, position = 'identity') +
  theme_ipsum()

unique(sort(nnpp$trap_id))
unique(sort(nnpp$trap_number))

# create grounded_bat var ----
unique(sort(nnpp$survey_organiser))
unique(sort(nnpp$equipment_used))

nnpp %<>% 
  dplyr::mutate(grounded_bat = dplyr::case_when(
    survey_organiser %in% c("N/A - grounded bat", "N/A - grounded bat record", "N/A - grounded bat release") ~ TRUE,
    equipment_used %in% c("Bat found on wind turbine") ~ TRUE,
    TRUE ~ FALSE))

# check sensible temperatures ----
unique(sort(nnpp$start_temp_C))
janitor::tabyl(nnpp, start_temp_C, month)

nnpp %>% 
  ggplot(., aes(x = month, y = start_temp_C, colour = factor(year)), group = factor(year))+
  geom_point(aes(colour = factor(year)), position = position_dodge(0.4))+
  geom_smooth(method = "loess", span = 1.5, position = position_dodge(0.4))+
  stat_summary(fun.data = "mean_cl_boot", position = position_dodge(0.4))+
  facet_wrap(~ year)

nnpp %<>% 
  dplyr::mutate(start_temp_C = dplyr::case_when(start_temp_C == 0 ~ NA, TRUE ~ start_temp_C))

unique(sort(nnpp$end_temp_C))
janitor::tabyl(nnpp, end_temp_C, month)

nnpp %<>% 
  dplyr::mutate(end_temp_C = dplyr::case_when(end_temp_C == 0 ~ NA, TRUE ~ end_temp_C))

nnpp %>% 
  ggplot(., aes(x = month, y = end_temp_C, colour = factor(year)), group = factor(year))+
  geom_point(aes(colour = factor(year)), position = position_dodge(0.4))+
  geom_smooth(method = "loess", span = 1.5, position = position_dodge(0.4))+
  stat_summary(fun.data = "mean_cl_boot", position = position_dodge(0.4))+
  facet_wrap(~ year)

# identifies error in year (2002); investigate ----
# bat_trap_date correct - alter others
nnpp %<>% 
  dplyr::mutate(
    survey_date_start = dplyr::case_when(
      lubridate::year(survey_date_start) == 2002 ~ `year<-`(survey_date_start , 2023),TRUE ~ survey_date_start),
    survey_date_end = dplyr::case_when(
      lubridate::year(survey_date_end) == 2002 ~ `year<-`(survey_date_end , 2023), TRUE ~ survey_date_end),
    year = dplyr::case_when(year == 2002 ~ 2023, TRUE ~ year))

nnpp %>% 
  ggplot(., aes(x = month, y = start_temp_C, colour = factor(year)), group = factor(year))+
  geom_point(aes(colour = factor(year)), position = position_dodge(0.4))+
  geom_smooth(method = "loess", span = 1.5, position = position_dodge(0.4))+
  stat_summary(fun.data = "mean_cl_boot", position = position_dodge(0.4))+
  facet_wrap(~ year)

# homogenise and codify cloud cover ----
unique(sort(nnpp$cloud_cover))

nnpp %<>% 
  dplyr::mutate(cloud_cover = dplyr::case_when(
    cloud_cover == "Clear (clear-1/3)" | cloud_cover == "Clear(clear-1/3)" | cloud_cover == "Clear(clear)-1/3)" ~ "Clear",
    cloud_cover == "Patchy (1/3-2/3)" | cloud_cover == "Patchy(1/3-2/3)" ~ "Patchy",
    cloud_cover == "Full (2/3-complete)" | cloud_cover == "Full(2/3-complete)" ~ "Full",
    cloud_cover == 0 | is.na(cloud_cover) ~ "Unknown",
    TRUE ~ cloud_cover))

nnpp %>%
  ggplot(aes(x = cloud_cover)) +
  geom_histogram(fill = "#404080", colour = "#e9ecef", alpha = 0.6, position = 'identity', stat = "count") +
  theme_ipsum()

nnpp %>%
  ggplot(aes(x = cloud_cover)) +
  geom_histogram(fill = "#404080", colour = "#e9ecef", alpha = 0.6, position = 'identity', stat = "count") +
  theme_ipsum()+
  facet_wrap(~year)

# homogenise and codify wind ----
unique(sort(nnpp$wind))

nnpp %<>% 
  dplyr::mutate(wind = dplyr::case_when(
    wind == "light breeze" | wind == "Light breeze" | wind == "Light Breeze" ~ "Light",
    wind == 0 | is.na(wind) ~ "Unknown",
    TRUE ~ wind))

nnpp %>%
  ggplot(aes(x = wind)) +
  geom_histogram(fill = "#404080", colour = "#e9ecef", alpha = 0.6, position = 'identity', stat = "count") +
  theme_ipsum()

nnpp %>%
  ggplot(aes(x = wind)) +
  geom_histogram(fill = "#404080", colour = "#e9ecef", alpha = 0.6, position = 'identity', stat = "count") +
  theme_ipsum()+
  facet_wrap(~year)+
  theme(axis.text.x=element_text(angle = 45, hjust = 1))

# homogenise and codify rain ----
unique(sort(nnpp$rain))

nnpp %<>% 
  dplyr::mutate(rain = dplyr::case_when(
    rain == "Constant drizzle" ~ "Drizzle",
    is.na(rain) | rain == "0" ~ "Unknown",
    TRUE ~ rain))

nnpp %>%
  ggplot(aes(x = rain)) +
  geom_histogram(fill = "#404080", colour = "#e9ecef", alpha = 0.6, position = 'identity', stat = "count") +
  theme_ipsum()

nnpp %>%
  ggplot(aes(x = rain)) +
  geom_histogram(fill = "#404080", colour = "#e9ecef", alpha = 0.6, position = 'identity', stat = "count") +
  theme_ipsum()+
  facet_wrap(~year)+
  theme(axis.text.x=element_text(angle = 45, hjust = 1))

# homogenise and codify moon ----
unique(sort(nnpp$moon))

nnpp %<>% 
  dplyr::mutate(moon = dplyr::case_when(
    moon == "0" & is.na(start_temp_C) ~ "Unknown",
    is.na(moon) ~ "Unknown",
    moon == "Full" ~ "1",
    moon == "Half" ~ "0.25",
    moon == "New" ~ "0",
    moon == "Quarter" ~ "0.25",
    moon == "Three quarter" | moon == "Three quarters" | moon == "Three quartters" ~ "0.75",
    TRUE ~ moon))

unique(sort(nnpp$moon))

nnpp %>%
  ggplot(aes(x = moon)) +
  geom_histogram(fill = "#404080", colour = "#e9ecef", alpha = 0.6, position = 'identity', stat = "count") +
  theme_ipsum()

nnpp %>%
  ggplot(aes(x = moon)) +
  geom_histogram(fill = "#404080", colour = "#e9ecef", alpha = 0.6, position = 'identity', stat = "count") +
  theme_ipsum()+
  facet_wrap(~year)+
  theme(axis.text.x=element_text(angle = 45, hjust = 1))

# create moon phase var ----
nnpp %<>% 
  dplyr::mutate(moon_phase = dplyr::case_when(
    moon == "0" ~ "New",
    moon == "0.25" ~ "Wax-Wane Crescent",
    moon == "0.5" ~ "First-Last Quarter",
    moon == "0.75" ~ "Wax-Wane Gibbous",
    moon == "1" ~ "Full",
    TRUE ~ moon))

unique(sort(nnpp$moon_phase))

# check no zero bat surveys lost in processing (should = 78) ----
length(nnpp$zero_bat_survey[nnpp$zero_bat_survey == TRUE])

# create simpler equipment var for comparison of trapping method ----
unique(sort(nnpp$equipment_used))

nnpp %<>% 
  dplyr::mutate(equipment_simple = dplyr::case_when(
    stringr::str_detect(equipment_used, "(?i)net") ~ "Net",
    stringr::str_detect(equipment_used, "(?i)mist") ~ "Net",
    stringr::str_detect(equipment_used, "(?i)6m") ~ "Net",
    stringr::str_detect(equipment_used, "(?i)box") ~ "Bat Box",
    stringr::str_detect(equipment_used, "(?i)schwegler") ~ "Bat Box",
    stringr::str_detect(equipment_used, "(?i)2f-dfp") ~ "Bat Box",
    stringr::str_detect(equipment_used, "(?i)harp") ~ "Harp Trap",
    stringr::str_detect(equipment_used, "(?i)trap") ~ "Harp Trap",
    stringr::str_detect(equipment_used, "(?i)bank") ~ "Harp Trap",
    stringr::str_detect(equipment_used, "(?i)HT") ~ "Harp Trap",
    stringr::str_detect(equipment_used, "(?i)triple") ~ "Harp Trap",
    stringr::str_detect(equipment_used, "(?i)double") ~ "Harp Trap",
    stringr::str_detect(equipment_used, "(?i)at100") ~ "Harp Trap",
    stringr::str_detect(equipment_used, "(?i)none") ~ "None (Found)",
    stringr::str_detect(equipment_used, "(?i)found") ~ "None (Found)",
    equipment_used == "0" | is.na(equipment_used) ~ "Unknown",
    TRUE ~ "Unknown"))

unique(sort(nnpp$equipment_simple))

nnpp %>%
  ggplot(aes(x = equipment_simple)) +
  geom_histogram(fill = "#404080", colour = "#e9ecef", alpha = 0.6, position = 'identity', stat = "count") +
  theme_ipsum()

nnpp %>%
  ggplot(aes(x = equipment_simple)) +
  geom_histogram(fill = "#404080", colour = "#e9ecef", alpha = 0.6, position = 'identity', stat = "count") +
  theme_ipsum()+
  facet_wrap(~year)+
  theme(axis.text.x=element_text(angle = 45, hjust = 1))

# homogenise and codify habitat ----
unique(sort(nnpp$habitat))

nnpp %<>% 
  dplyr::mutate(habitat_simple = dplyr::case_when(
    habitat == "Built up areas and gardens" ~ "Built-up areas and gardens",
    habitat == "Broadleaved and mixed woodland" | habitat == "Broadleaved or mixed woodland" | habitat == "Mixed Broadleaf" | habitat == "Mixed Broadleaf Woodland"  ~ "Broadleaf/Mixed Woodland",
    habitat == "Mixed Woodland / Water edge" | habitat == "Mixed woodland / river" | habitat == "Mixed woodland / lake side"  ~ "Mixed Woodland/Water",
    habitat == "Willlow trees" | habitat == "Mature Willow" ~ "Willow",
    habitat == "Willow scrub" ~ "Willow Scrub",
    habitat == "Large willow near shore" ~ "Willow/Water",
    habitat == "Broadleaf Woodland/water edge" ~ "Broadleaf/Water",
    habitat == "Mixed Broadleaf/River bank" ~ "Mixed Broadleaf/Water",
    habitat == "Lake edge" | habitat == "Lake surrounded by parkland and wood" ~ "Lake",
    habitat == "Coastal" ~ "Coast",
    habitat == "Mixed Broadleaf / Wetland" | habitat == "Wetland / broadleaf trees" ~ "Mixed Broadleaf/Wetland",
    habitat == "Conifer trees / Water edge"  ~ "Coniferous Woodland/Water",
    habitat == "Coniferous woodland"  ~ "Coniferous Woodland",
    habitat == "Broad leaf hedgerow/River bank"  ~ "Broadleaf Hedgerow/Water",
    habitat == "Built-up areas and gardens" ~ "Urban",
    habitat == "Improved grassland" | habitat == "Semi-natural grassland" ~ "Grassland",
    habitat == "Mixed scrub" ~ "Scrub",
    habitat == "Mixed woodland" ~ "Mixed Woodland",
    habitat == "Woodland edge/ Mature Trees" ~ "Woodland Edge",
    habitat == "Wetland / Grassland" ~ "Grassland/Wetland",
    habitat == 0 | is.na(habitat) | habitat == "unknown" ~ "Unknown",
    TRUE ~ habitat))

unique(sort(nnpp$habitat_simple))

nnpp %>%
  ggplot(aes(x = habitat_simple)) +
  geom_histogram(fill = "#404080", colour = "#e9ecef", alpha = 0.6, position = 'identity', stat = "count") +
  theme_ipsum()+
  theme(axis.text.x=element_text(angle = 45, hjust = 1))

nnpp %>%
  ggplot(aes(x = habitat_simple)) +
  geom_histogram(fill = "#404080", colour = "#e9ecef", alpha = 0.6, position = 'identity', stat = "count") +
  theme_ipsum()+
  facet_wrap(~year)+
  theme(axis.text.x=element_text(angle = 45, hjust = 1))

# check no zero bat surveys lost in processing (should = 78) ----
length(nnpp$zero_bat_survey[nnpp$zero_bat_survey == TRUE])

# survey count = 1500 ----

# save cleaned NNPP survey data V1 ----
saveRDS(nnpp, "~/MEGA/Bat Conservation Trust/NNPP/rds/nnpp.rds")

# export nnpp (clean) as .csv for visualisation in GIS (WAS NNPP1 IN V01 > nnpp_clean_03 >> nnpp2) ----
data.table::fwrite(nnpp, "~/MEGA/Bat Conservation Trust/NNPP/gis/nnpp.csv",
                   sep = ",", na = "", dec = ".", row.names = FALSE, col.names = TRUE,
                   dateTimeAs = "ISO")

# clean erroneous trap locations ----

# issues with survey_id: (easting, northing, lat, lon)
# 685 Druridge Bay CP-Ladyburn Lake Stepping Stones - correct: 426900 600200 55.295182 -1.5779225
# 1758 Druridge Bay Ladyburn Lake Visitor Centre - correct: 427340 600080 55.29408 -1.5710052
# 678 Druridge Bay CP-Ladyburn Lake Stepping Stones - correct: 426900 600200 55.295182 -1.5779225
# 676 Druridge Bay CP-Ladyburn Lake Stepping Stones - correct: 426900 600200 55.295182 -1.5779225
# 1856 Filby Broad - correct: 645830 313150 52.66029 1.6337067
# 1075 Titchwell RSPB Reserve - correct: easting northing 52.962588 0.606348
# 1058 Lily Broad - correct: 645520 314750 52.674785 1.6303235
# 1396 Filby Broad - correct: 645830 313150 52.66029 1.6337067
# 1397 Filby Broad - correct: 645830 313150 52.66029 1.6337067

nnpp %<>% 
  dplyr::mutate(
    trap_x = dplyr::case_when(
      survey_id == 685 & trap_x != 426900 ~ 426900,
      survey_id == 1758 & trap_x != 427340 ~ 427340,
      survey_id == 678 & trap_x != 426900  ~ 426900,
      survey_id == 676 & trap_x != 426900 ~ 426900,
      survey_id == 1856 & trap_x != 645830 ~ 645830,
      survey_id == 1075  ~ NA_real_,
      survey_id == 1058 & trap_x != 645520 ~ 645520,
      survey_id == 1396 & trap_x != 645830 ~ 645830,
      survey_id == 1397 & trap_x != 645830 ~ 645830,
      TRUE ~ trap_x),
    
    trap_y = dplyr::case_when(
      survey_id == 685 & trap_y != 600200 ~ 600200,
      survey_id == 1758 & trap_y != 600080 ~ 600080,
      survey_id == 678 & trap_y != 600200  ~ 600200,
      survey_id == 676 & trap_y != 600200 ~ 600200,
      survey_id == 1856 & trap_y != 313150  ~ 313150 ,
      survey_id == 1075  ~ NA_real_,
      survey_id == 1058 & trap_y != 314750 ~ 314750,
      survey_id == 1396 & trap_y != 313150 ~ 313150,
      survey_id == 1397 & trap_y != 313150 ~ 313150,
      TRUE ~ trap_y),
    
    trap_lat = dplyr::case_when(
      survey_id == 685 & trap_lat != 55.295182 ~ 55.295182,
      survey_id == 1758 & trap_lat != 55.29408 ~ 55.29408,
      survey_id == 678 & trap_lat != 55.295182  ~ 55.295182,
      survey_id == 676 & trap_lat != 55.295182 ~ 55.295182,
      survey_id == 1856 & trap_lat != 52.66029 ~ 52.66029,
      survey_id == 1075  ~ 52.962588,
      survey_id == 1058 & trap_lat != 52.674785 ~ 52.674785,
      survey_id == 1396 & trap_lat != 52.66029 ~ 52.66029,
      survey_id == 1397 & trap_lat != 52.66029 ~ 52.66029,
      TRUE ~ trap_lat),
    
    trap_lon = dplyr::case_when(
      survey_id == 685 & trap_lon != -1.5779225 ~ -1.5779225,
      survey_id == 1758 & trap_lon != -1.5710052 ~ -1.5710052,
      survey_id == 678 & trap_lon != -1.5779225  ~ -1.5779225,
      survey_id == 676 & trap_lon != -1.5779225 ~ -1.5779225,
      survey_id == 1856 & trap_lon != 1.6337067 ~ 1.6337067,
      survey_id == 1075  ~ 0.606348,
      survey_id == 1058 & trap_lon != 1.6303235 ~ 1.6303235,
      survey_id == 1396 & trap_lon != 1.6337067 ~ 1.6337067,
      survey_id == 1397 & trap_lon != 1.6337067 ~ 1.6337067,
      TRUE ~ trap_lon))

# save cleaned NNPP survey data V2 ----
saveRDS(nnpp, "~/MEGA/Bat Conservation Trust/NNPP/rds/nnpp.rds")

# check no zero bat surveys lost in processing (should = 78) ----
length(nnpp$zero_bat_survey[nnpp$zero_bat_survey == TRUE])

# import bat species names ----
bat_names <- readr::read_csv("~/MEGA/Bat Conservation Trust/NNPP/data/input/bat_species_names.csv")
saveRDS(bat_names, "~/MEGA/Bat Conservation Trust/NNPP/rds/bat_names.rds")

# bind Linnaean names to nnpp ----
nnpp <- 
  dplyr::left_join(nnpp, bat_names, by = "species") %>% 
  dplyr::mutate(species = stringr::str_replace(species, "sp.", "spp."))


# create var to indicate if data are in/from the official NNPP ----

# drop non-NNPP survey data (dropping 'found' bats for this df only - retain for distribution etc...)
# drop bat box finds (occurred during day and not strictly part of NNPP surveying)
nnpp %<>% 
  dplyr::mutate(data_origin = dplyr::case_when(year < 2014 ~ "pre_survey",
                                               equipment_simple == "Bat Box" ~ "adhoc",
                                               stringr::str_detect(site_name, pattern = "non-project") ~ "non_survey",
                                               TRUE ~ "nnpp_survey"))

# n 49421 ----

# save cleaned NNPP survey data V3 ----
saveRDS(nnpp, "~/MEGA/Bat Conservation Trust/NNPP/rds/nnpp.rds")

###############################################################################-
# • • • National Nathusius' Pipistrelle Project - SURVEY DATA ONLY #############
###############################################################################-

# extract NNPP survey data ONLY (902 rows dropped) ----
nnpp_survey <- nnpp %>% dplyr::filter(data_origin == "nnpp_survey")

# n = 48519 ----
# survey count = 1434 ----

# calculate trapping effort per survey ----

### Rationale: 

# Each survey_id is unique, thus 'survey' is (or would appear to be) the observational unit.
# However, at least 2 key variables (related to trapping effort) vary between surveys: number of traps deployed and for how long?
# 1. Deployment duration is a problem with these data because many survey records (228 (~ 16%)) omit this information ----

nnpp_survey %>%
  dplyr::distinct(survey_id, survey_time_start, survey_time_end) %>% 
  dplyr::mutate(duration_info = dplyr::case_when(!is.na(survey_time_start) & !is.na(survey_time_end) ~ "y",
                                                 TRUE ~ "n")) %>%
  janitor::tabyl(., duration_info)

# calculate survey duration where possible ----
# correct some obviously erroneous timestamps ----
# e.g. surveys starting at 10:00/11:00 should be 22:00/23:00
nnpp_survey_durations <- nnpp_survey %>%
  dplyr::distinct(survey_id, survey_datetime_start, survey_datetime_end, survey_time_start, survey_time_end) %>% 
  dplyr::mutate(duration_info = dplyr::case_when(!is.na(survey_time_start) & !is.na(survey_time_end) ~ "y",
                                                 TRUE ~ "n")) %>% 
  dplyr::mutate(survey_datetime_start2 = dplyr::case_when(survey_time_start < "14:00" ~ survey_datetime_start + lubridate::hours(12),
                                                          TRUE ~ survey_datetime_start),
                survey_datetime_end2 = dplyr::case_when(survey_time_start < "14:00" ~ survey_datetime_end + lubridate::hours(12),
                                                        TRUE ~ survey_datetime_end),
                dt_flag = dplyr::case_when(survey_datetime_end2 < survey_datetime_start2 ~ "flag"),
                survey_datetime_end3 = dplyr::case_when(dt_flag == "flag" ~ survey_datetime_end2 + lubridate::hours(24),
                                                        TRUE ~ survey_datetime_end2),
                survey_datetime_start2 = dplyr::case_when(survey_id == 1955 ~ survey_datetime_start2 - lubridate::hours(24),
                                                          TRUE ~ survey_datetime_start2),
                survey_duration = dplyr::case_when(duration_info == "y" ~ difftime(survey_datetime_end3, survey_datetime_start2,
                                                                                   units = "mins"))) %>% 
  dplyr::rename(survey_datetime_start_raw = survey_datetime_start,
                survey_datetime_end_raw = survey_datetime_end,
                survey_time_start_raw = survey_time_start,
                survey_time_end_raw = survey_time_end,
                survey_datetime_start_clean = survey_datetime_start2,
                survey_datetime_end_clean = survey_datetime_end3) %>% 
  dplyr::select(-duration_info, -dt_flag, -survey_datetime_end2) %>% 
  dplyr::mutate(survey_duration_inferred = dplyr::case_when(is.na(survey_duration) ~ "yes",
                                                            TRUE ~ "no"))

janitor::tabyl(nnpp_survey_durations$survey_duration_inferred)

# draw a pseudo-random sample from complete survey durations ----
# draw in a way that matches the distribution of actual survey durations 
nnpp_survey_durations %>% 
  dplyr::filter(!is.na(survey_duration)) %>%  
  ggplot(., aes(x = survey_duration))+
  geom_histogram(binwidth = 15, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  theme_ipsum()

nnpp_survey_durations %>% 
  dplyr::filter(!is.na(survey_duration)) %>%  
  ggplot(., aes(x = 1, y = survey_duration))+
  geom_violin(fill="#69b3a2", color="#e9ecef", alpha=0.9)+
  geom_boxplot(width=0.1, color="grey", alpha=0.2)+
  theme_ipsum()+
  coord_flip()

# distribution approximates normal
summary(nnpp_survey_durations$survey_duration)

fit <- kde1d::kde1d(as.numeric(nnpp_survey_durations$survey_duration))
plot(fit)
simulations <- kde1d::rkde1d(228, fit)
simulations
plot(simulations)
sims_df <- as.data.frame(simulations)

ggplot(sims_df, aes(x = simulations))+
  geom_histogram(binwidth = 15, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  theme_ipsum()

simulations <- as.integer(simulations)
simulations <- lubridate::dminutes(simulations)

# insert inferred survey durations in to actual survey durations df ----
ind <- which(is.na(nnpp_survey_durations$survey_duration))
nnpp_survey_durations[ind, "survey_duration"] <- simulations

nnpp_survey_durations %>% 
  ggplot(., aes(x = survey_duration))+
  geom_histogram(binwidth = 15, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  theme_ipsum()

# clean workspace
rm(ind, simulations, fit, sims_df)
gc()

# 2. Differing number of traps deployed per survey alters survey effort ----

# extract number of traps deployed for each survey
nnpp_traps_per_survey <- nnpp_survey %>%
  dplyr::group_by(survey_id, zero_bat_survey) %>% 
  dplyr::summarise(trap_count_per_survey = length(unique(trap_id)))

# establish mean trap per survey
nnpp_traps_per_survey %>% 
  ggplot(., aes(x = trap_count_per_survey))+
  geom_histogram(binwidth = 1, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  theme_ipsum()

# apply mean trap number (mean = 2) to zero bat count surveys
nnpp_traps_per_survey %<>% 
  dplyr::mutate(trap_count_per_survey = dplyr::case_when(zero_bat_survey == TRUE ~ 2, TRUE ~ trap_count_per_survey))

# bind to survey durations df
nnpp_survey_durations <- dplyr::left_join(nnpp_survey_durations, dplyr::select(nnpp_traps_per_survey, -zero_bat_survey), by = "survey_id")

# calculate survey effort in minutes 
nnpp_survey_durations %<>% dplyr::mutate(survey_effort_mins = survey_duration * trap_count_per_survey)

# # create survey effort master key ----
# nnpp_survey_effort_master_01 <- nnpp_survey_durations
# 
# # save nnpp survey effort master key
# saveRDS(nnpp_survey_effort_master_01 , "~/MEGA/Bat Conservation Trust/Nathusius' Pipistrelle Project/rds/nnpp_survey_effort_master_01.rds")

# 1. continued... Time errors carry over to bat_trap_time etc. fix these ----
nnpp_bat_trap_datetimes_clean_01 <- nnpp_survey %>% dplyr::distinct(survey_id, bat_trap_datetime, bat_trap_time)

nnpp_bat_trap_datetimes_clean_01 <- dplyr::left_join(nnpp_bat_trap_datetimes_clean_01,
                                                     dplyr::select(nnpp_survey_durations, survey_id, survey_datetime_start_clean),
                                                     by = "survey_id")

nnpp_bat_trap_datetimes_clean_01 %<>%
  dplyr::mutate(bat_time_flag = dplyr::case_when(bat_trap_datetime < survey_datetime_start_clean ~ "flag"))

nnpp_bat_trap_datetimes_clean_01 %<>%
  dplyr::mutate(bat_time_flag2 =
                  dplyr::case_when(bat_time_flag == "flag" &
                                     lubridate::date(bat_trap_datetime) == lubridate::date(survey_datetime_start_clean) ~ "flag"))

nnpp_bat_trap_datetimes_clean_01 %<>%
  dplyr::mutate(bat_trap_datetime2 = dplyr::case_when(bat_time_flag2 == "flag" ~ bat_trap_datetime + lubridate::hours(24),
                                                      TRUE ~ bat_trap_datetime))
nnpp_bat_trap_datetimes_clean_01 %<>%
  dplyr::mutate(bat_trap_datetime2 = dplyr::case_when(survey_id == 1635 ~ bat_trap_datetime %m+% months(2),
                                                      survey_id == 165 ~ bat_trap_datetime + lubridate::days(5),
                                                      TRUE ~ bat_trap_datetime2))

nnpp_bat_trap_datetimes_clean_01 %<>%
  dplyr::mutate(bat_time_flag3 = dplyr::case_when(bat_trap_datetime2 < survey_datetime_start_clean ~ "flag"))

nnpp_bat_trap_datetimes_clean_01 %<>%
  dplyr::select(survey_id, bat_trap_datetime, bat_trap_datetime2) %>% 
  dplyr::rename(bat_trap_datetime_clean = bat_trap_datetime2)

# # save corrected bat_trap_datetimes
# saveRDS(nnpp_bat_trap_datetimes_clean_01, "~/MEGA/Bat Conservation Trust/Nathusius' Pipistrelle Project/rds/nnpp_bat_trap_datetimes_clean_01.rds")
# 
# # clean workspace
# rm(nnpp_survey_durations, nnpp_traps_per_survey)
# gc()

# merge all to create clean master nnpp survey df ----
janitor::compare_df_cols(nnpp_survey, nnpp_survey_durations, nnpp_bat_trap_datetimes_clean_01)

nnpp_survey_durations %<>% dplyr::ungroup()

nnpp_survey <- dplyr::left_join(nnpp_survey,
                                dplyr::select(nnpp_survey_durations, -survey_duration, -site_name),
                                by = "survey_id")

nnpp_survey <- dplyr::left_join(nnpp_survey,
                                 nnpp_bat_trap_datetimes_clean_01,
                                 by = c("survey_id", "bat_trap_datetime"))

# data.table::fwrite(nnpp_survey, "nnpp_survey.csv")

# reconcile and rename corrected and erroneous records (datetimes) ----
nnpp_survey %<>% 
  dplyr::rename(survey_date_start_raw = survey_date_start,
                survey_date_end_raw = survey_date_end,
                bat_trap_date_raw = bat_trap_date,
                bat_trap_time_raw = bat_trap_time,
                bat_trap_datetime_raw = bat_trap_datetime) %>% 
  dplyr::mutate(survey_date_start_clean = 
                  dplyr::case_when(is.na(survey_datetime_start_clean) ~ survey_date_start_raw,
                                   TRUE ~ lubridate::date(survey_datetime_start_clean)),
                
                survey_time_start_clean = hms::as_hms(survey_datetime_start_clean),
                
                survey_date_end_clean = 
                  dplyr::case_when(is.na(survey_datetime_end_clean) ~ survey_date_end_raw,
                                   TRUE ~ lubridate::date(survey_datetime_end_clean)),
                
                survey_time_end_clean = hms::as_hms(survey_datetime_end_clean),
                
                bat_trap_date_clean = 
                  dplyr::case_when(is.na(bat_trap_datetime_clean) ~ bat_trap_date_raw,
                                   TRUE ~ lubridate::date(bat_trap_datetime_clean)),
                
                bat_trap_time_clean = hms::as_hms(bat_trap_datetime_clean),
                year = lubridate::year(survey_date_start_clean),
                month = lubridate::month(survey_date_start_clean),
                day = lubridate::day(survey_date_start_clean),
                month_fc = as.factor(month),
                survey_effort_mins_num = as.numeric(survey_effort_mins)) %>% 
  dplyr::select(survey_id, trap_count_per_survey, survey_effort_mins, bat_group_id, bat_group_name,
                site_id, site_name, year, month, day,
                survey_datetime_start_clean,
                survey_date_start_clean,
                survey_time_start_clean,
                survey_datetime_end_clean,
                survey_date_end_clean,
                survey_time_end_clean,
                bat_trap_datetime_clean,
                bat_trap_date_clean,
                bat_trap_time_clean,
                species, species_name, sex, age, female_reproductive_status, ring_number,
                bat_record_id:habitat_simple, survey_duration_inferred,
                survey_datetime_start_raw,
                survey_date_start_raw,
                survey_time_start_raw,
                survey_datetime_end_raw,
                survey_date_end_raw,
                survey_time_end_raw,
                bat_trap_datetime_raw,
                bat_trap_date_raw,
                bat_trap_time_raw, everything()) %>% 
  dplyr::select(-survey_time_start, - survey_datetime_start, -survey_time_end, -survey_datetime_end)

# save nnpp survey master ----
saveRDS(nnpp_survey, "~/MEGA/Bat Conservation Trust/NNPP/rds/nnpp_survey.rds")

# calculate number of bats captured per survey (Pipistrelle species ONLY) ---- 
# ensure all surveys are retained (need those in which zero NP were trapped)

# extract Nathusisus' pipistrelle (NP) data only

length(unique(nnpp_survey$survey_id))

temp_np <- nnpp_survey %>%
  dplyr::group_by(survey_id) %>% 
  dplyr::filter(row_number() == 1| zero_bat_survey == TRUE | species == "Nathusius' pipistrelle")

duplicates <- temp_np %>% 
  group_by(bat_record_id) %>% 
  filter(n() > 1) %>% 
  dplyr::select(survey_id, species, zero_bat_survey, bat_record_id) %>% 
  ungroup()

rm(duplicates)

length(unique(temp_np$survey_id))

temp_np %>% janitor::tabyl(. ,species, zero_bat_survey)

temp_np %<>% dplyr::select(survey_id, species, zero_bat_survey)

temp_np %<>% 
  dplyr::group_by(survey_id, species) %>% 
  dplyr::summarise(np_captures = n())

temp_np %<>% 
  dplyr::group_by(survey_id) %>%
  dplyr::mutate(np_captures2 = dplyr::case_when(n() == 1 & species != "Nathusius' pipistrelle" ~ 0,
                                                n() == 1 & species == "Nathusius' pipistrelle" ~ np_captures,
                                                n() > 1 & species != "Nathusius' pipistrelle" ~ 0,
                                                is.na(species) ~ 0,
                                                TRUE ~ np_captures)) %>% 
  dplyr::group_by(survey_id) %>% 
  dplyr::mutate(np_captures3 = sum(np_captures2)) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(survey_id, np_captures3) %>% 
  dplyr::rename(np_captures = np_captures3) %>% 
  dplyr::distinct()

temp_np %<>% 
  dplyr::left_join(., dplyr::select(nnpp_survey, survey_id, zero_bat_survey), by = "survey_id") %>%
  dplyr::distinct() %>% 
  dplyr::select(survey_id, zero_bat_survey, np_captures)

sum(temp_np$np_captures)

# extract Soprano pipistrelle (SP) data only
temp_sp <- nnpp_survey %>%
  dplyr::group_by(survey_id) %>% 
  dplyr::filter(row_number() == 1| zero_bat_survey == TRUE | species == "Soprano pipistrelle")

duplicates <- temp_sp %>% 
  group_by(bat_record_id) %>% 
  filter(n() > 1) %>% 
  dplyr::select(survey_id, species, zero_bat_survey, bat_record_id) %>% 
  ungroup()

rm(duplicates)

length(unique(temp_sp$survey_id))

temp_sp %>% janitor::tabyl(. ,species, zero_bat_survey)

temp_sp %<>% dplyr::select(survey_id, species, zero_bat_survey)

temp_sp %<>% 
  dplyr::group_by(survey_id, species) %>% 
  dplyr::summarise(sp_captures = n())

temp_sp %<>% 
  dplyr::group_by(survey_id) %>%
  dplyr::mutate(sp_captures2 = dplyr::case_when(n() == 1 & species != "Soprano pipistrelle" ~ 0,
                                                n() == 1 & species == "Soprano pipistrelle" ~ sp_captures,
                                                n() > 1 & species != "Soprano pipistrelle" ~ 0,
                                                is.na(species) ~ 0,
                                                TRUE ~ sp_captures)) %>% 
  dplyr::group_by(survey_id) %>% 
  dplyr::mutate(sp_captures3 = sum(sp_captures2)) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(survey_id, sp_captures3) %>% 
  dplyr::rename(sp_captures = sp_captures3) %>% 
  dplyr::distinct()

temp_sp %<>% 
  dplyr::left_join(., dplyr::select(nnpp_survey, survey_id, zero_bat_survey), by = "survey_id") %>%
  dplyr::distinct() %>% 
  dplyr::select(survey_id, zero_bat_survey, sp_captures)

sum(temp_sp$sp_captures)

# extract Common pipistrelle (CP) data only
temp_cp <- nnpp_survey %>%
  dplyr::group_by(survey_id) %>% 
  dplyr::filter(row_number() == 1| zero_bat_survey == TRUE | species == "Common pipistrelle")

duplicates <- temp_cp %>% 
  group_by(bat_record_id) %>% 
  filter(n() > 1) %>% 
  dplyr::select(survey_id, species, zero_bat_survey, bat_record_id) %>% 
  ungroup()

rm(duplicates)

length(unique(temp_cp$survey_id))

temp_cp %>% janitor::tabyl(. ,species, zero_bat_survey)

temp_cp %<>% dplyr::select(survey_id, species, zero_bat_survey)

temp_cp %<>% 
  dplyr::group_by(survey_id, species) %>% 
  dplyr::summarise(cp_captures = n())

temp_cp %<>% 
  dplyr::group_by(survey_id) %>%
  dplyr::mutate(cp_captures2 = dplyr::case_when(n() == 1 & species != "Common pipistrelle" ~ 0,
                                                n() == 1 & species == "Common pipistrelle" ~ cp_captures,
                                                n() > 1 & species != "Common pipistrelle" ~ 0,
                                                is.na(species) ~ 0,
                                                TRUE ~ cp_captures)) %>% 
  dplyr::group_by(survey_id) %>% 
  dplyr::mutate(cp_captures3 = sum(cp_captures2)) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(survey_id, cp_captures3) %>% 
  dplyr::rename(cp_captures = cp_captures3) %>% 
  dplyr::distinct()

temp_cp %<>% 
  dplyr::left_join(., dplyr::select(nnpp_survey, survey_id, zero_bat_survey), by = "survey_id") %>%
  dplyr::distinct() %>% 
  dplyr::select(survey_id, zero_bat_survey, cp_captures)

sum(temp_cp$cp_captures)

nnpp_survey_pipistrelle <- temp_np %>% 
  dplyr::left_join(temp_sp, by = c("survey_id", "zero_bat_survey")) %>% 
  dplyr::left_join(temp_cp, by = c("survey_id", "zero_bat_survey"))

rm(temp_np, temp_sp, temp_cp)

# create a survey characteristics key to merge with 'nnpp_survey_pipistrelle' ----
nnpp_survey_characteristics <- nnpp_survey %>% 
  dplyr::select(survey_id, survey_duration_inferred, trap_count_per_survey, survey_effort_mins, zero_bat_survey,
                site_name, site_id, bat_group_id, year, month, day, survey_datetime_start_clean, survey_datetime_end_clean,
                survey_date_start_clean, survey_time_start_clean, survey_date_end_clean, survey_time_end_clean,
                start_temp_C, end_temp_C, cloud_cover, wind, rain, moon, moon_phase) %>% 
  dplyr::distinct()

# create SITE lat and lon ----
survey_lat_lon_centroid <- nnpp_survey %>%
  dplyr::ungroup() %>% 
  dplyr::select(survey_id, trap_lat, trap_lon) %>% 
  dplyr::distinct() %>% 
  dplyr::group_by(survey_id) %>% 
  dplyr::summarise(site_lat = mean(trap_lat),
                   site_lon = mean(trap_lon))

nnpp_survey_characteristics <- dplyr::left_join(nnpp_survey_characteristics,
                                                survey_lat_lon_centroid,
                                                by = "survey_id") %>% 
  dplyr::ungroup()

# bind survey characteristics (and SITE lat lon) to pipistrelle survey data ----
nnpp_survey_pipistrelle <- dplyr::left_join(nnpp_survey_pipistrelle,
                                            dplyr::select(nnpp_survey_characteristics, -zero_bat_survey),
                                            by = "survey_id")

# bind survey characteristics (and SITE lat lon) to full NNPP survey data ----
janitor::compare_df_cols(nnpp_survey, nnpp_survey_characteristics)

nnpp_survey <- dplyr::left_join(nnpp_survey,
                                dplyr::select(nnpp_survey_characteristics, survey_id, site_lat, site_lon),
                                by = "survey_id")

# save survey characterstics ----
saveRDS(nnpp_survey_characteristics, "~/MEGA/Bat Conservation Trust/NNPP/rds/nnpp_survey_characteristics.rds")

# save nnpp_survey ----
saveRDS(nnpp_survey, "~/MEGA/Bat Conservation Trust/NNPP/rds/nnpp_survey.rds")

# save nnpp_survey_pipistrelle ----
saveRDS(nnpp_survey_pipistrelle, "~/MEGA/Bat Conservation Trust/NNPP/rds/nnpp_survey_pipistrelle.rds")

rm(survey_lat_lon_centroid)
gc()

# create var indicating if survey site is within/out 10 km of coast ----

# import UK shapefile ----
GBR <- sf::read_sf("~/MEGA/Bat Conservation Trust/NNPP/gis/gadm41_GBR.gpkg")
class(GBR)
head(GBR)
sf::st_crs(GBR)
# set GBR CRS to 27700 (OS Grid) ----
GBR_27700 <- GBR %>% sf::st_transform(crs = 27700)

# import Jersey shapefile ----
JEY <- sf::read_sf("~/MEGA/Bat Conservation Trust/NNPP/gis/gadm41_JEY.gpkg")
class(JEY)
head(JEY)
sf::st_crs(JEY)
# set Jersey CRS to 27700 (OS Grid) ----
JEY_27700 <- JEY %>% sf::st_transform(crs = 27700)

# bind GBR and Jersey shapefiles ----
UK_27700 <- rbind(GBR_27700, JEY_27700)

rm(GBR, GBR_27700, JEY, JEY_27700)
gc()

# examine surveys (pipistrelle only) that lack location data (only 2 and they make sense) ----
nnpp_survey_pipistrelle %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(is.na(site_lat)) %>% 
  dplyr::select(survey_id, site_name) %>% 
  dplyr::distinct()

# remove surveys lacking location data before conversion to sf object ----
nnpp_survey_pipistrelle_sf <- nnpp_survey_pipistrelle %>% dplyr::filter(!is.na(site_lon))

# process pipistrelle data for spatial query ----
nnpp_survey_pipistrelle_sf <- sf::st_as_sf(nnpp_survey_pipistrelle_sf, coords = c("site_lon", "site_lat"), remove = F)
class(nnpp_survey_pipistrelle_sf)
head(nnpp_survey_pipistrelle_sf)
sf::st_crs(nnpp_survey_pipistrelle_sf)

# set survey location CRS to 27700 (OS Grid) ----
nnpp_survey_pipistrelle_sf_4326 <- nnpp_survey_pipistrelle_sf %>% sf::st_set_crs(4326)
nnpp_survey_pipistrelle_sf_27700 <- nnpp_survey_pipistrelle_sf_4326 %>% sf::st_transform(crs = 27700)

# check UK shapefile and survey site data have same CRS ----
sf::st_crs(UK_27700)
sf::st_crs(nnpp_survey_pipistrelle_sf_27700)

# create a 10 km buffer inside UK coast ----
buffer <- sf::st_buffer(UK_27700, -10000)

ggplot()+
  geom_sf(data = UK_27700)+
  geom_sf(data = buffer, fill = "#74A089", alpha = 0.4)+
  geom_sf(data = nnpp_survey_pipistrelle_sf_27700, size = 0.75, colour = "firebrick")

temp_buffer_intersection <- sf::st_intersection(nnpp_survey_pipistrelle_sf_27700, buffer)

temp_buffer_intersection %<>% dplyr::mutate(gis = "intersect?") 

nnpp_survey_pipistrelle <- dplyr::left_join(nnpp_survey_pipistrelle,
                                            dplyr::select(temp_buffer_intersection, survey_id, gis),
                                            by = "survey_id")
nnpp_survey_pipistrelle %<>%
  dplyr::distinct() %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(coastal = dplyr::case_when(gis == "intersect?" ~ "inland",
                                           is.na(site_lat) ~ NA_character_,
                                           TRUE ~ "coastal"))

nnpp_survey_characteristics <- dplyr::left_join(nnpp_survey_characteristics,
                                                dplyr::select(temp_buffer_intersection, survey_id, gis),
                                                by = "survey_id")

nnpp_survey <- dplyr::left_join(nnpp_survey,
                                dplyr::select(temp_buffer_intersection, survey_id, gis),
                                by = "survey_id")

# save survey characterstics ----
saveRDS(nnpp_survey_characteristics, "~/MEGA/Bat Conservation Trust/NNPP/rds/nnpp_survey_characteristics.rds")

# save nnpp_survey ----
saveRDS(nnpp_survey, "~/MEGA/Bat Conservation Trust/NNPP/rds/nnpp_survey.rds")

# save nnpp_survey_pipistrelle ----
saveRDS(nnpp_survey_pipistrelle, "~/MEGA/Bat Conservation Trust/NNPP/rds/nnpp_survey_pipistrelle.rds")

rm(temp_buffer_intersection, nnpp_survey_pipistrelle_sf, nnpp_survey_pipistrelle_sf_4326)
gc()

######################## GTH ###################################################












# calculate scaled body mass index for Nathusius' pipistrelle ----
nnpp_sbmi <- nnpp_survey %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(!is.na(bat_mass_g)) %>% 
  dplyr::select(survey_id, bat_mass_g, forearm_length_cm, species, sex, age, female_reproductive_status, bat_record_id)

nrow(nnpp_sbmi); length(unique(nnpp_sbmi$bat_record_id))

nnpp_sbmi <- dplyr::left_join(nnpp_sbmi,
                              dplyr::select(nnpp_survey_pipistrelle, survey_id, site_lat, site_lon, coastal),
                              by = "survey_id") %>% 
  dplyr::distinct()

p <- nnpp_sbmi %>%
  dplyr::filter(species == "Nathusius' pipistrelle") %>% 
  ggplot(., aes(x = forearm_length_cm, y = bat_mass_g))+
  geom_point()+
  stat_smooth(method = lm, alpha = 0.3)+
  ggpubr::stat_regline_equation(aes(label = ..rr.label..),
                                colour="black",
                                label.x = -Inf, label.y = Inf,
                                vjust = 1.5, hjust = -0.1)
plotly::ggplotly(p)

p1 <- nnpp_sbmi %>%
  dplyr::filter(species == "Common pipistrelle") %>% 
  ggplot(., aes(x = forearm_length_cm, y = bat_mass_g))+
  geom_point()+
  stat_smooth(method = lm, alpha = 0.3)+
  ggpubr::stat_regline_equation(aes(label = ..rr.label..),
                                colour="black",
                                label.x = -Inf, label.y = Inf,
                                vjust = 1.5, hjust = -0.1)
plotly::ggplotly(p1)

p2 <- nnpp_sbmi %>%
  dplyr::filter(species == "Soprano pipistrelle") %>% 
  ggplot(., aes(x = forearm_length_cm, y = bat_mass_g))+
  geom_point()+
  stat_smooth(method = lm, alpha = 0.3)+
  ggpubr::stat_regline_equation(aes(label = ..rr.label..),
                                colour="black",
                                label.x = -Inf, label.y = Inf,
                                vjust = 1.5, hjust = -0.1)
plotly::ggplotly(p2)

nnpp_sbmi2 <- nnpp_sbmi %>% dplyr::filter(!(bat_record_id %in% c(18916, 18713, 18929, 14005, 6608)))

nnpp_sbmi2 %>%
  dplyr::filter(species == "Nathusius' pipistrelle") %>% 
  ggplot(., aes(x = forearm_length_cm, y = bat_mass_g))+
  geom_point()+
  stat_smooth(method = lm, alpha = 0.3)+
  ggpubr::stat_regline_equation(aes(label = ..rr.label..),
                                colour="black",
                                label.x = -Inf, label.y = Inf,
                                vjust = 1.5, hjust = -0.1)
nnpp_sbmi2 %>%
  ggplot(., aes(x = forearm_length_cm, y = bat_mass_g, colour = species))+
  geom_point()+
  stat_smooth(method = lm, alpha = 0.3)+
  ggpubr::stat_regline_equation(aes(label = ..rr.label..),
                                colour="black",
                                label.x = -Inf, label.y = Inf,
                                vjust = 1.5, hjust = -0.1)+
  facet_wrap(~species)

nnpp_sbmi2 %>%
  dplyr::filter(species %in% c("Nathusius' pipistrelle", "Soprano pipistrelle", "Common pipistrelle")) %>% 
  ggplot(., aes(x = forearm_length_cm, y = bat_mass_g, colour = species))+
  geom_point()+
  stat_smooth(method = lm, alpha = 0.3)+
  ggpubr::stat_regline_equation(aes(label = ..rr.label..),
                                colour="black",
                                label.x = -Inf, label.y = Inf,
                                vjust = 1.5, hjust = -0.1)+
  facet_wrap(~species, ncol = 2)

# conduct standardised major axis (SMA) regression on morphometrics

# Fit a line to the natural logs of M and L
# The slope estimate is used as the power function (bSMA) in the scaled mass index (SMI) calculation.

# extract species and metrics for analysis (NP) ----
np_sbmi <- nnpp_sbmi2 %>% dplyr::filter(species == "Nathusius' pipistrelle")

np_l_sma <- smatr::sma(bat_mass_g ~ forearm_length_cm, log="xy", data = np_sbmi, method = "SMA", alpha = 0.05)
plot(np_l_sma, which = "default", )
np_l_sma
np_l_sma_out <- np_l_sma$groupsummary
np_l_bSMA <- as.numeric(np_l_sma_out$Slope)

# calculate mean of L
np_l_meanL <- mean(np_sbmi$forearm_length_cm, na.rm = T)

# calculate scaled mass index (SMI) for each individual 
np_sbmi <- np_sbmi %>% dplyr::mutate(smi = bat_mass_g * (np_l_meanL/forearm_length_cm)^np_l_bSMA)

# extract species and metrics for analysis (CP) ----
cp_sbmi <- nnpp_sbmi2 %>% dplyr::filter(species == "Common pipistrelle")

cp_l_sma <- smatr::sma(bat_mass_g ~ forearm_length_cm, log="xy", data = cp_sbmi, method = "SMA", alpha = 0.05)
plot(cp_l_sma, which = "default", )
cp_l_sma
cp_l_sma_out <- cp_l_sma$groupsummary
cp_l_bSMA <- as.numeric(cp_l_sma_out$Slope)

# calculate mean of L
cp_l_meanL <- mean(cp_sbmi$forearm_length_cm, na.rm = T)

# calculate scaled mass index (SMI) for each individual 
cp_sbmi <- cp_sbmi %>% dplyr::mutate(smi = bat_mass_g * (cp_l_meanL/forearm_length_cm)^cp_l_bSMA)

# extract species and metrics for analysis (SP) ----
sp_sbmi <- nnpp_sbmi2 %>% dplyr::filter(species == "Soprano pipistrelle")

sp_l_sma <- smatr::sma(bat_mass_g ~ forearm_length_cm, log="xy", data = sp_sbmi, method = "SMA", alpha = 0.05)
plot(sp_l_sma, which = "default", )
sp_l_sma
sp_l_sma_out <- sp_l_sma$groupsummary
sp_l_bSMA <- as.numeric(sp_l_sma_out$Slope)

# calculate mean of L
sp_l_meanL <- mean(sp_sbmi$forearm_length_cm, na.rm = T)

# calculate scaled mass index (SMI) for each individual 
sp_sbmi <- sp_sbmi %>% dplyr::mutate(smi = bat_mass_g * (sp_l_meanL/forearm_length_cm)^sp_l_bSMA)

# bind 3 species sbmi ----
np_sp_cp_smbi <- dplyr::bind_rows(np_sbmi, cp_sbmi, sp_sbmi)

saveRDS(np_sp_cp_smbi, "~/MEGA/Bat Conservation Trust/Nathusius' Pipistrelle Project/rds/np_sp_cp_smbi_01.rds")

np_sp_cp_smbi %>% 
  dplyr::filter(sex != "Unknown") %>%
  ggplot(., aes(x = sex, y = smi))+
  geom_jitter(width = 0.1, aes(colour = sex), alpha = 0.7, size = 0.7)+
  stat_summary(aes(x = sex, y = smi), size = 0.6, 
               fun=mean,
               fun.min = function(x) mean(x) - sd(x), 
               fun.max = function(x) mean(x) + sd(x), 
               geom = "pointrange")+
  labs(x = "Sex", y = "Scaled Mass Index (g) Mean ± SD")+
  scale_colour_manual(values = c("#46ACC8", "#E2D200"),
                      breaks = c("Female", "Male"))+
  coord_cartesian(ylim = c(0, 15))+
  theme_ipsum()+
  theme(legend.position = "none",
        axis.text.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))+
  facet_wrap(~ species)

