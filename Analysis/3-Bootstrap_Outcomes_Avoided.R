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
mods <- readRDS(here::here("data/best_mods.rds"))

case_mod <- paste0("cases ~ ", mods$cases)
hosp_mod <- paste0("hosps ~ ", mods$hosps)
mort_mod <- paste0("deaths ~ ", mods$deaths)

dev <- F

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

  RNGkind("L'Ecuyer-CMRG")
  set.seed(430)

   if(dev) nboots <- 100 ; cat("Running", nboots, "bootstrap simulations\n")
  
# Run bootstrapping in parallel  \
  
clooster <- parallel::makeCluster(4)


parallel::clusterExport(clooster, c("its_df_quarts", "zs", "nboots", "case_mod", "hosp_mod", "mort_mod"))
parallel::clusterEvalQ(clooster, {
  library(tidyverse)
  library(splines)
  library(fastglm)
})

registerDoParallel(clooster)

if(dev) tictoc::tic()

its_quart_boots <- foreach(n = 1:nboots) %dopar% { #bind_rows(parallel::parLapply(clooster, 1:nboots, function(b){
#  its_quart_boots <- lapply(1:nboots, function(b){
    
    #tic()
    z_samps <- sample(zs, size = length(zs), replace = T)
    
    b_df <- bind_rows(lapply(z_samps, function(z){
      its_df_quarts[its_df_quarts$zip == z,]
    }))
    
    b_pop = sum(b_df %>% filter(date == min(date)) %>% pull(pop))
    
    # Cases model and cases avoided  
    cases_mm <- model.matrix(as.formula(case_mod), data = b_df)
    cases_mm0 <- cases_mm
    cases_mm0[,which(grepl("int", colnames(cases_mm)))] <- 0
    
    c_mod <- fastglm(x = cases_mm,# model.matrix(as.formula(case_mod), data = b_df)
                     y = b_df$cases,
                     family = poisson(link = "log"),
                     offset = log(b_df$pop),
                     method = 3)
    
    b_df$pred_c1 <- c_mod$fitted.values
    b_df$pred_c0 <- exp(drop(cases_mm0 %*% c_mod$coefficients) + log(b_df$pop))
    
    # c_mod2 <- glm(as.formula(case_mod),
    #               family = poisson(link = "log"),
    #               offset = log(pop),
    #               data = b_df)
    # 
    # View(cbind(c_mod$coefficients, c_mod2$coefficients))
    # View(cbind(cases_mm[,which(grepl("bs", colnames(cases_mm)))][,12], bs(b_df$t_wk, knots = seq(min(b_df$t_wk), max(b_df$t_wk), by = 3))[,12]))
    # View(cbind(c_mod$fitted.values, 
    #            exp(predict(c_mod2)), 
    #            exp(drop(cases_mm %*% c_mod$coefficients) + log(b_df$pop))))
    # View(cbind(exp(predict(c_mod2, newdata = b_df %>% mutate(int = 0))), 
    #            exp(drop(cases_mm0 %*% c_mod$coefficients) + log(b_df$pop))))
    # 
    # b_df$pred_c1 <- exp(predict(c_mod2))
    # b_df$pred_c0 <- exp(predict(c_mod2, newdata = b_df %>% mutate(int = 0)))
    
    # Hosps model and hosps avoided  
    hosps_mm <- model.matrix(as.formula(hosp_mod), data = b_df)
    hosps_mm0 <- hosps_mm
    hosps_mm0[,which(grepl("int", colnames(hosps_mm)))] <- 0
    
    h_mod <- fastglm(x = hosps_mm, 
                     y = b_df$hosps,
                     family = poisson(link = "log"),
                     offset = log(b_df$pop),
                     method = 3)
    
    b_df$pred_h1 <- h_mod$fitted.values
    b_df$pred_h0 <- exp(drop(hosps_mm0 %*% h_mod$coefficients) + log(b_df$pop))
    
    # h_mod <- glm(as.formula(hosp_mod),
    #              family = poisson(link = "log"),
    #              offset = log(pop),
    #              data = b_df)
    # 
    # b_df$pred_h1 <- exp(predict(h_mod))
    # b_df$pred_h0 <- exp(predict(h_mod, newdata = b_df %>% mutate(int = 0)))
    

    # Deaths model and deaths avoided  
    deaths_mm <- model.matrix(as.formula(mort_mod), data = b_df)
    deaths_mm0 <- deaths_mm
    deaths_mm0[,which(grepl("int", colnames(deaths_mm)))] <- 0
    
    d_mod <- fastglm(x = deaths_mm, 
                     y = b_df$deaths,
                     family = poisson(link = "log"),
                     offset = log(b_df$pop),
                     method = 3)
    
    b_df$pred_d1 <- d_mod$fitted.values
    b_df$pred_d0 <- exp(drop(deaths_mm0 %*% d_mod$coefficients) + log(b_df$pop))
    
    # d_mod <- glm(as.formula(mort_mod),
    #              family = poisson(link = "log"),
    #              offset = log(pop),
    #              data = b_df)
    # 
    # b_df$pred_d1 <- exp(predict(d_mod))
    # b_df$pred_d0 <- exp(predict(d_mod, newdata = b_df %>% mutate(int = 0)))
    
    out <- b_df %>% 
             group_by(date, county, vemquartile) %>% 
             summarise(c_obs = sum(cases),
                       h_obs = sum(hosps),
                       d_obs = sum(deaths),
                       c_mod = sum(pred_c1),
                       h_mod = sum(pred_h1),
                       d_mod = sum(pred_d1),
                       c_avd_obs = sum(pred_c0-cases),
                       h_avd_obs = sum(pred_h0-hosps),
                       d_avd_obs = sum(pred_d0-deaths),
                       c_avd = sum(pred_c0-pred_c1),
                       h_avd = sum(pred_h0-pred_h1),
                       d_avd = sum(pred_d0-pred_d1)) #%>% 
             #dplyr::select(zip:date,pred_c1:pred_d0) %>% 
             #mutate(iter = b)
    
    #toc()
    
    return(out)
  
  }#)

if(dev) tictoc::toc()
  
# its_quart_boots <- bind_rows(lapply(1:nboots, function(b){
#   
#   z_samps <- sample(zs, size = length(zs), replace = T)
#   
#   b_df <- bind_rows(lapply(z_samps, function(z){
#     its_df_quarts[its_df_quarts$zip == z,]
#   }))
#   
#   b_pop = sum(b_df %>% filter(date == min(date)) %>% pull(pop))
#   
# # Cases model and cases avoided  
#   c_mod <- glm(cases ~ county + bs(t_wk, knots = seq(min(its_df_quarts$t_wk), max(its_df_quarts$t_wk), by = 2)) + factor(vemquartile, levels = c("4", "3", "2", "1")) + tests_p100k + int + int:bs(t_wk, knots = seq(min(its_df_quarts$t_wk), max(its_df_quarts$t_wk), by = 2)),
#                family = poisson(link = "log"),
#                offset = log(pop),
#                data = b_df)
#   
#   b_df$pred_c1 <- exp(predict(c_mod))
#   b_df$pred_c0 <- exp(predict(c_mod, newdata = b_df %>% mutate(int = 0)))
#   
#   #c_est <- sum(b_df %>% filter(int == 1) %>% mutate(c_avd = pred_c0 - pred_c1) %>% pull(c_avd)) 
# 
# # Hosps model and hosps avoided  
#   h_mod <- glm(hosps ~ county + bs(t_wk, knots = seq(min(its_df_quarts$t_wk), max(its_df_quarts$t_wk), by = 2)) + factor(vemquartile, levels = c("4", "3", "2", "1")) + int + int:bs(t_wk, knots = seq(min(its_df_quarts$t_wk), max(its_df_quarts$t_wk), by = 2)),
#                family = poisson(link = "log"),
#                offset = log(pop),
#                data = b_df)
#   
#   b_df$pred_h1 <- exp(predict(h_mod))
#   b_df$pred_h0 <- exp(predict(h_mod, newdata = b_df %>% mutate(int = 0)))
#   
#   #h_est <- sum(b_df %>% filter(int == 1) %>% mutate(h_avd = pred_h0 - pred_h1) %>% pull(h_avd)) 
#   
# # Deaths model and deaths avoided  
#   d_mod <- glm(deaths ~ county + bs(t_wk, knots = seq(min(its_df_quarts$t_wk), max(its_df_quarts$t_wk), by = 2)) + factor(vemquartile, levels = c("4", "3", "2", "1")) + int + int:bs(t_wk, knots = seq(min(its_df_quarts$t_wk), max(its_df_quarts$t_wk), by = 2)),
#                family = poisson(link = "log"),
#                offset = log(pop),
#                data = b_df)
#   
#   b_df$pred_d1 <- exp(predict(d_mod))
#   b_df$pred_d0 <- exp(predict(d_mod, newdata = b_df %>% mutate(int = 0)))
#   
#   #d_est <- sum(b_df %>% filter(int == 1) %>% mutate(d_avd = pred_d0 - pred_d1) %>% pull(d_avd)) 
#   
#   
#   # b_est <- b_df %>% filter(int == 1) %>% 
#   #   mutate(pred1 = b_pred1,
#   #          pred0 = b_pred0,
#   #          IRR = pred1/pred0,
#   #          IRD = pred1-pred0,
#   #          wt = pop/b_pop,
#   #          IRR_wt = IRR*wt,
#   #          boot = b)
#   # 
#   # cat(b, "... ")
#   
#   return(b_df %>% 
#            group_by(date) %>% 
#            summarise(c_avd = sum(pred_c0-pred_c1),
#                      h_avd = sum(pred_h0-pred_h1),
#                      d_avd = sum(pred_d0-pred_d1)) %>% 
#            #dplyr::select(zip:date,pred_c1:pred_d0) %>% 
#            mutate(iter = b))
# }))

parallel::stopCluster(clooster)

for(l in 1:length(its_quart_boots)){
  its_quart_boots[[l]] <- its_quart_boots[[l]] %>% mutate(iter = l)
}

saveRDS(bind_rows(its_quart_boots), here::here(paste0("data/bootstrapped_outcomes_avoided_its", nboots,".rds")))
