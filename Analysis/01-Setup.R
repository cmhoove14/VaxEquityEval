source(here::here("Utils/Plot_Utils.R"))

# Key values ------------------
policy_date <- ymd(20210301)
its_start_date <- ymd(20201214)# First week vaccines really started ramping up as opposed to ymd(20210103) which is just beginning of new year
its_end_date <- ymd(20211101) # Prior to widespread booster rollout, omicron emergence, and ~ 8 months after policy when vaxes from policy may be waning ymd(20210830)
nboots <- 10000

# Data ---------------------
zip_fin <- readRDS(here::here("data/zip_fin_vem.rds")) %>% 
  filter(!is.na(vemquartile)) %>% 
  group_by(zip) %>% 
  arrange(date) %>% 
  mutate(cases_cum = cumsum(cases),
         prop_case = cases_cum/popTOT,
         hosps_cum = cumsum(hosp),
         logit_prop_case = log(prop_case/(1-prop_case))) %>% 
  ungroup()

# view(zip_fin %>% filter(prop_case > 1)) 

# Zip code 90071 in downtown LA has pop=126 but cases > 400 for some reason, so filter it out

# zip_fin <- zip_fin %>% filter(zip != "90071")

#excess_cases_zips <- unique(zip_fin %>% filter(date <= ymd(20220101)) %>% filter(prop_case > 1) %>% pull(zip))
#no_cases_zips <- unique(zip_fin %>% filter(date >= ymd(20210103) & date <= ymd(20220101) & prop_case <= 0) %>% pull(zip))
no_pop_zips <- unique(zip_fin %>% filter(vempop <= 0) %>% pull(zip))

no_case_zips <- zip_fin %>% 
  filter(date < its_end_date) %>% 
  group_by(zip) %>% 
  summarise(tot_cases = sum(cases)) %>% 
  filter(tot_cases <= 0) %>% 
  pull(zip)

zip_fin <- zip_fin %>% 
  filter(!zip %in% c(no_pop_zips,no_case_zips))#c(excess_cases_zips, no_cases_zips, no_pop_zips))

weekly_zip_sum <- readRDS(here::here("data/weekly_zip_sum.RDS")) %>% 
  filter(!zip %in% c(no_pop_zips, no_case_zips))#c(excess_cases_zips, no_cases_zips, no_pop_zips))

# County population from zip level dataset  
cnty_pop <- zip_fin %>% 
  filter(county!= "", !is.na(county), date == min(date)) %>% 
  group_by(county) %>% 
  summarise(pop = sum(vempop))

# County population from zip level dataset  
cnty_vem_pop <- zip_fin %>% 
  filter(county!= "", !is.na(county), date == min(date)) %>% 
  group_by(county, vemquartile) %>% 
  summarise(pop = sum(vempop))

vem_pop <- cnty_vem_pop %>% 
  group_by(vemquartile) %>% 
  summarise(pop = sum(pop))

ca_pop <- sum(cnty_pop$pop)
