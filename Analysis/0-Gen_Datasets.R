# ------------------------------------------------
# Generate base data 
# Chris Hoover, Jan 2022
# ------------------------------------------------

library(tidyverse)
library(tidycensus)
library(lubridate)
library(data.table)

# Data ---------------------------
# Master zip-code level file with zip-level population, tests, cases, hospitalizations, deaths, and vaccinations. Amalgamation of many data sources constructed in CMH_Daily.Rmd

zip_master <- fread("//mnt/projects/connect-izb/resources/general/master_zip_by_date_vem_from_linelist.tsv") %>% 
  rename("vempop" = totalpop)

# doses by zip 
vax_zip <- fread("//mnt/projects/connect-izb/rmd_generate/vaccinations/vax_by_agegrp_zip_date.csv") %>% 
  filter(!is.na(zip)) %>% 
  mutate(zip = as.character(zip)) %>% 
  group_by(zip,date) %>% 
  summarise(vax = sum(vax)) %>% 
  ungroup()

vax1_zip <- fread("//mnt/projects/connect-izb/rmd_generate/vaccinations/vax1_by_agegrp_zip_date.csv") %>% 
  filter(!is.na(zip)) %>% 
  mutate(zip = as.character(zip)) %>% 
  group_by(zip,date) %>% 
  summarise(vax1 = sum(vax)) %>% 
  ungroup()

#vax_zip %>% group_by(date) %>% summarise(vax = sum(vax)) %>% ggplot(aes(x = date, y = vax))  +geom_rect(xmin = policy_date-28, xmax = policy_date+21, ymin = -Inf, ymax = Inf, alpha = 0.4, fill = "red")+geom_line()

# Vaccine equity metric
vem <- fread("//mnt/projects/epi/covid19/R/resources/progData/vem_allvariables_2021.csv") %>% 
  dplyr::select(zip, vem, vemquartile, county, totalpop) %>% 
  mutate(zip = as.character(zip))

vem_octiles <- quantile(vem$vem, seq(0,1,by = 1/8))
vem_20ths <- quantile(vem$vem, seq(0,1,by = 1/20))

# Vaccines delivered
vax_dlvr <- fread("https://data.chhs.ca.gov/dataset/ead44d40-fd63-4f9f-950a-3b0111074de8/resource/4ebd81ad-1ab9-4a3d-aeda-8a0f06cfee7b/download/covid19vaccinesshippeddeliveredonhandbyzip.csv") %>% mutate(zip = as.character(zip))

# Tests conducted
tests <- fread("//mnt/projects/connect-izb/rmd_generate/geocoding/zip_calc.tsv") %>% 
  dplyr::select(date, "zip" = location, tests_new)

# ACS 2021 population by ZCTA
acs20pop <- get_acs(geography = "zcta", 
                    table = "B01001", 
                    year = 2020,
                    survey = "acs5")

acs_vars <- load_variables(year = 2020, dataset = "acs5") %>% 
  filter(grepl("B01001_", name))

zcta_pop <- acs20pop %>% 
  mutate(zip = str_remove(NAME, "ZCTA5 ")) %>% 
  dplyr::select(-c(GEOID, NAME)) %>% 
  filter(zip %in% vem$zip) %>% 
  left_join(acs_vars, by = c("variable" = "name")) %>% 
  filter(!variable %in% c("B01001_002", "B01001_026")) %>% 
  mutate(AGEGRP = case_when(grepl("Under 5|5 to 9|10 to 14|15 to 17", label) ~ "u18",
                            grepl("18 and 19|20 years|21 years|22 to 24|25 to 29|30 to 34|35 to 39|40 to 44|45 to 49", label) ~ "18-49",
                            grepl("50 to 54|55 to 59|60 and 61|62 to 64", label) ~ "50-64",
                            grepl("65 and 66|67 to 69|70 to 74|75 to 79|80 to 84|85 years", label) ~ "65+",
                            TRUE ~ "TOT")) %>% 
  group_by(zip, AGEGRP) %>% 
  summarise(pop = sum(estimate)) %>% 
  pivot_wider(values_from = pop, names_from = AGEGRP, names_prefix = "pop")

# Put all together 
zip_fin <- zip_master %>% mutate(zip = as.character(zip)) %>% #filter(zip %in% vem$zip) %>% dplyr::select(-pop) %>% 
  left_join(tests %>% mutate(zip = as.character(zip)), by = c("zip", "date")) %>% 
  left_join(vax_dlvr %>% dplyr::select(date, zip, doses_delivered), by = c("zip", "date")) %>% 
  left_join(vax_zip, by = c("zip", "date")) %>% 
  left_join(vax1_zip, by = c("zip", "date")) %>% 
  left_join(zcta_pop, by = "zip") %>% 
  mutate(doses_delivered = replace_na(doses_delivered, 0),
         vax = replace_na(vax, 0),
         vax1 = replace_na(vax1, 0),
         vemq1 = if_else(vemquartile == 1, 1, 0),
         vemoctile = cut(vem, vem_octiles, label = F),
         vem20ths = cut(vem, vem_20ths, label = F))

saveRDS(zip_fin, here::here("data/zip_fin_vem.rds"))

# County population from zip level dataset  
cnty_pop <- zip_fin %>% 
  filter(county!= "", !is.na(county), date == min(date)) %>% 
  group_by(county) %>% 
  summarise(pop = sum(vempop))

saveRDS(cnty_pop, here::here("data/county_pop.RDS"))

# County population from zip level dataset  
cnty_vem_pop <- zip_fin %>% 
  filter(county!= "", !is.na(county), date == min(date)) %>% 
  group_by(county, vemquartile) %>% 
  summarise(pop = sum(vempop))

saveRDS(cnty_vem_pop, here::here("data/county_vem_pop.RDS"))


# Vaccines delivered 
vax_dlvr_fin <- vax_dlvr %>% 
  left_join(vem %>% dplyr::select(zip,vemquartile, totalpop), by = "zip") %>% 
  left_join(zcta_pop, by = "zip") %>% 
  mutate(week = floor_date(date, "week")) %>% 
  group_by(week, zip) %>% 
  summarise(county = first(county),
            vemquartile = first(vemquartile),
            doses_delivered = sum(doses_delivered),
            pop = first(totalpop)) %>% 
  ungroup() %>% 
  mutate(dd_100k = doses_delivered/pop*1e5)
  
saveRDS(vax_dlvr_fin, here::here("data/vax_delivered.RDS"))

# Summarise to month
# monthly_zip_sum <- zip_fin %>% 
#   #st_drop_geometry() %>% 
#   filter(!is.na(vemquartile)) %>% 
#   group_by(zip,moyr) %>% 
#   summarise(date = first(date),
#             county = first(county),
#             region = first(region),
#             pop = first(totalpop),
#             vemq1 = mean(vemq1), # should be 0 or 1, but use this just to check
#             vemquartile = first(vemquartile),
#             vemoctile = first(vemoctile),
#             vem20th = first(vem20ths),
#             tests = sum(tests_new),
#             vax = sum(vax),
#             cases = sum(cases_new),
#             hosps = sum(hosp),
#             deaths = sum(deaths_new)) %>% 
#   ungroup() %>% 
#   mutate(testp100k = tests/pop*1e5,
#          vaxp100k = vax/pop*1e5,
#          casep100k = cases/pop*1e5,
#          hospsp100k = hosps/pop*1e5,
#          deathp100k = deaths/pop*1e5) %>% 
#   arrange(zip, date)
# 
# saveRDS(monthly_zip_sum, here::here("data/monthly_zip_sum.RDS"))
# 
# Summarise to week
weekly_zip_sum <- zip_fin %>% 
  filter(date >= ymd(20200301) & date <= ymd(20220101) & !is.na(vemquartile)) %>% 
  mutate(wkyr = floor_date(date, "week")+1, # Want mondays rather than sundays
         t_wk = as.numeric(wkyr - ymd(20210301))/7) %>% 
  group_by(zip,t_wk) %>% 
  summarise(date = first(wkyr),
            county = first(county),
            region = first(region),
            pop = first(vempop),
            pop50up = first(`pop50-64`)+first(`pop65+`),
            pop65up = first(`pop65+`),
            vemq1 = mean(vemq1), # should be 0 or 1, but use this just to check
            vemquartile = first(vemquartile),
            vemoctile = first(vemoctile),
            vem20th = first(vem20ths),
            vemscore = first(vem),
            tests = sum(tests_new),
            vax = sum(vax),
            vax1 = sum(vax1),
            doses_dd = sum(doses_delivered),
            cases = sum(cases),
            hosps = sum(hosp),
            deaths = sum(mort)) %>% 
  ungroup() %>% 
  group_by(zip) %>% 
  arrange(date) %>% 
  mutate(cases_cum = cumsum(cases),
         vax_cum = cumsum(vax),
         vax1_cum = cumsum(vax1),
         prop_case = cases_cum/pop,
         logit_prop_case = log(prop_case/(1-prop_case)),
         hosps_cum = cumsum(hosps),
         deaths_cum = cumsum(deaths)) %>% 
  ungroup() %>% 
  mutate(per50up = pop50up/pop*100,
         per65up = pop65up/pop*100,
         testp100k = tests/pop*1e5,
         vaxp100k = vax/pop*1e5,
         vax1p100k = vax1/pop*1e5,
         cumvaxp100k = vax_cum/pop*1e5,
         cumvax1p100k = vax1_cum/pop*1e5,
         cumcasep100k = cases_cum/pop*1e5,
         suscp100k = (pop-vax_cum/2-cases_cum)/pop*1e5,
         cumhospp100k = hosps_cum/pop*1e5,
         cumdeathp100k = deaths_cum/pop*1e5,
         casep100k = cases/pop*1e5,
         hospsp100k = hosps/pop*1e5,
         deathp100k = deaths/pop*1e5) %>% 
  arrange(zip, date)

saveRDS(weekly_zip_sum, here::here("data/weekly_zip_sum.RDS"))
