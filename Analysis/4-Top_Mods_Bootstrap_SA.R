library(tidyverse)
library(lubridate)
#library(data.table)
library(splines)
library(parallel)
library(doParallel)
library(foreach)
library(fastglm)
library(tictoc)

source(here::here("Analysis/01-Setup.R"))
mod_prfrm <- readRDS(here::here("data/mods_mean_mses.rds"))

mod_prfrm_cases_mean <- rowMeans(mod_prfrm$case_mses)
mod_prfrm_hosps_mean <- rowMeans(mod_prfrm$hosp_mses)
mod_prfrm_morts_mean <- rowMeans(mod_prfrm$mort_mses)

mod_prfrm_cases_min <- min(mod_prfrm_cases_mean)
mod_prfrm_hosps_min <- min(mod_prfrm_hosps_mean)
mod_prfrm_morts_min <- min(mod_prfrm_morts_mean)

cases_top_mods <- mod_prfrm$forms[order(mod_prfrm_cases_mean)][2:10] # Best mods already evaluated, so look at 2-10 here
hosps_top_mods <- mod_prfrm$forms[order(mod_prfrm_hosps_mean)][2:10]
morts_top_mods <- mod_prfrm$forms[order(mod_prfrm_morts_mean)][2:10]

topboots <- 100

cases_top_fill <- expand_grid(mod = cases_top_mods, iter = c(1:topboots)) %>% mutate(outcome = "cases", obs = NA_real_, pred0 = NA_real_)
hosps_top_fill <- expand_grid(mod = hosps_top_mods, iter = c(1:topboots)) %>% mutate(outcome = "hosps", obs = NA_real_, pred0 = NA_real_)
morts_top_fill <- expand_grid(mod = morts_top_mods, iter = c(1:topboots)) %>% mutate(outcome = "deaths", obs = NA_real_, pred0 = NA_real_)

its_df_quarts <- weekly_zip_sum %>% 
  filter(county != "") %>% 
  filter(date >= its_start_date & date < its_end_date) %>%  
  mutate(int = if_else(date >= policy_date & vemquartile == 1, 1, 0)) %>% 
  dplyr::select(zip,county,vemquartile,pop,per50up,per65up,date,t_wk,int,cases,hosps,deaths,tests,testp100k,cumcasep100k,cumvaxp100k)

# zips_exclude <- its_df_quarts %>% 
#   group_by(zip) %>% 
#   summarise(tot_cases = sum(cases), pop = first(pop)) %>% 
#   filter(pop > 499, tot_cases == 0) %>% 
#   pull(zip)

# its_df_quarts_use <- its_df_quarts %>% filter(!zip %in% zips_exclude) %>% filter(pop > 0)

zs <- unique(its_df_quarts$zip)

# Run bootstrapping 

for(b in 1:(topboots*length(cases_top_mods))){
  #tic()
  case_mod <- paste0("cases ~ ", cases_top_fill$mod[b])
  hosp_mod <- paste0("hosps ~ ", hosps_top_fill$mod[b])
  mort_mod <- paste0("deaths ~ ", morts_top_fill$mod[b])
  
  z_samps <- sample(zs, size = length(zs), replace = T)
  
  b_df <- bind_rows(lapply(z_samps, function(z){
    its_df_quarts[its_df_quarts$zip == z,]
  }))
  
  # Cases model and cases avoided  
  cases_mm <- model.matrix(as.formula(case_mod), data = b_df)
  cases_mm0 <- cases_mm
  cases_mm0[,which(grepl("int", colnames(cases_mm)))] <- 0
  
  c_mod <- fastglm(x = cases_mm,# model.matrix(as.formula(case_mod), data = b_df)
                   y = b_df$cases,
                   family = poisson(link = "log"),
                   offset = log(b_df$pop),
                   method = 3)
  
  b_df$pred_c0 <- exp(drop(cases_mm0 %*% c_mod$coefficients) + log(b_df$pop))
  
  # Hosps model and hosps avoided  
  hosps_mm <- model.matrix(as.formula(hosp_mod), data = b_df)
  hosps_mm0 <- hosps_mm
  hosps_mm0[,which(grepl("int", colnames(hosps_mm)))] <- 0
  
  h_mod <- fastglm(x = hosps_mm, 
                   y = b_df$hosps,
                   family = poisson(link = "log"),
                   offset = log(b_df$pop),
                   method = 3)
  
  b_df$pred_h0 <- exp(drop(hosps_mm0 %*% h_mod$coefficients) + log(b_df$pop))
  
  # Deaths model and deaths avoided  
  deaths_mm <- model.matrix(as.formula(mort_mod), data = b_df)
  deaths_mm0 <- deaths_mm
  deaths_mm0[,which(grepl("int", colnames(deaths_mm)))] <- 0
  
  d_mod <- fastglm(x = deaths_mm, 
                   y = b_df$deaths,
                   family = poisson(link = "log"),
                   offset = log(b_df$pop),
                   method = 3)
  
  b_df$pred_d0 <- exp(drop(deaths_mm0 %*% d_mod$coefficients) + log(b_df$pop))
  
  
  # Estimate VEM Q1 observed, counterfactual outcomes
  b_df_vemq1 <- b_df %>% 
    filter(vemquartile == 1, date >= policy_date)
    
  cases_top_fill$obs[b] <- sum(b_df_vemq1$cases)
  hosps_top_fill$obs[b] <- sum(b_df_vemq1$hosps)
  morts_top_fill$obs[b] <- sum(b_df_vemq1$deaths)
  
  cases_top_fill$pred0[b] <- sum(b_df_vemq1$pred_c0)
  hosps_top_fill$pred0[b] <- sum(b_df_vemq1$pred_h0)
  morts_top_fill$pred0[b] <- sum(b_df_vemq1$pred_d0)
  
  #toc()
  cat(b, ", ")
}

all_filled <- bind_rows(cases_top_fill, hosps_top_fill, morts_top_fill) %>% 
  mutate(avoided = pred0 - obs)

saveRDS(all_filled, here::here(paste0("data/bootstrap_top_mods_SA", topboots,".rds")))
