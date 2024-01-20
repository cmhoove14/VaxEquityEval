library(tidyverse)
library(lubridate)
library(data.table)
library(splines)
library(origami)
library(fastglm)
library(doParallel)

source(here::here("Analysis/01-Setup.R"))

n_rounds <- 10

log_trns <- function(v){
  out <- log(v)
  out[is.infinite(out)] <- NA_real_
  
  out
}

set.seed(430)

its_df_quarts <- weekly_zip_sum %>% 
  filter(county != "") %>% 
  filter(date >= its_start_date & date < its_end_date) %>%  
  mutate(int = if_else(date >= policy_date & vemquartile == 1, 1, 0),
         across(c(testp100k,cumcasep100k,cumvaxp100k),
                log_trns,
                .names = "log_{.col}"),
         across(c(testp100k,cumcasep100k,cumvaxp100k),
                .fn = function(x){scale(x)[,1]},
                .names = "scale_{.col}")) %>% 
  dplyr::select(zip,county,vemquartile,pop,per50up,per65up,date,t_wk,int,cases,hosps,deaths,tests,testp100k,cumcasep100k,cumvaxp100k,log_testp100k,log_cumcasep100k,log_cumvaxp100k,scale_testp100k,scale_cumcasep100k,scale_cumvaxp100k)

# zips_exclude <- its_df_quarts %>% 
#   group_by(zip) %>% 
#   summarise(tot_cases = sum(cases), pop = first(pop)) %>% 
#   filter(pop > 499, tot_cases == 0) %>% 
#   pull(zip)

its_df_quarts_largeCs <- its_df_quarts %>% #filter(!zip %in% zips_exclude) %>% 
  filter(county %in% cnty_pop$county[which(cnty_pop$pop > 100000)])

# Evaluate case models -------------------
# Candidate case GLMs
prop_forms_grid <- cross_df(list(base = "county + bs(t_wk, knots = seq(min(its_df_quarts$t_wk), max(its_df_quarts$t_wk), by = 3)) + factor(vemquartile, levels = c('4', '3', '2', '1')) + int",
                             int = c("", 
                                     "int:county", 
                                     "int:factor(vemquartile, levels = c('4', '3', '2', '1'))",
                                     "int:bs(t_wk, knots = seq(min(its_df_quarts$t_wk), max(its_df_quarts$t_wk), by = 3))"),
                             case = c("",
                                      #"scale_cumcasep100k",
                                      #"log_cumcasep100k",
                                      "cumcasep100k"),
                             vax = c("",
                                     #"scale_cumvaxp100k",
                                     #"log_cumvaxp100k",
                                     "cumvaxp100k"),
                             test = c("", 
                                      #"scale_testp100k",
                                      "testp100k"),
                             pop = c("", "per50up"))) %>% 
  arrange(base, int)

prop_formulas <- apply(prop_forms_grid, 1, function(i){
  str_c(i[i != ""], collapse = "+")
})


# Function to evaluate each formula in each fold
cv_glm <- function(fold, data, reg_form) {
  # get name and index of outcome variable from regression formula
  out_var <- as.character(unlist(str_split(reg_form, " "))[1])
  #out_var_ind <- as.numeric(which(colnames(data) == out_var))
  
  # split up data into training and validation sets
  train_data <- training(data)
    train_mm <- model.matrix(as.formula(reg_form), data = train_data)
  valid_data <- validation(data)
    valid_mm <- model.matrix(as.formula(reg_form), data = valid_data)
  
  # fit linear model on training set and predict on validation set
    mod <- fastglm(x = train_mm,
                   y = train_data[[out_var]],
                   family = poisson(link = "log"),
                   offset = log(train_data$pop),
                   method = 3)
  
  # Remove coefficients for groups not present in validation set
  coef_valid <- mod$coefficients
  coef_valid <- coef_valid[-which(!names(mod$coefficients) %in% colnames(valid_mm))]
  
  preds <- as.numeric(exp(drop(valid_mm %*% coef_valid) + log(valid_data$pop)))
  
  out_data <- valid_data
  out_data$preds <- preds
  out_data$SE <- as.numeric((out_data$preds - out_data[[out_var]])^2)
  
  # mod <- glm(as.formula(reg_form), 
  #            data = train_data,
  #            offset = log(pop),
  #            family = poisson(link = "log"))
  # preds <- predict(mod, newdata = valid_data, type = "response")

  return(out_data)
}

# Run cases model evaluation in parallel  \
clooster <- parallel::makeCluster(4)

RNGkind("L'Ecuyer-CMRG")
set.seed(430)

parallel::clusterExport(clooster, c("its_df_quarts", "its_df_quarts_largeCs", "n_rounds", "cv_glm", "prop_formulas"))
parallel::clusterEvalQ(clooster, {
  library(tidyverse)
  library(splines)
  library(fastglm)
  library(origami)
})

registerDoParallel(clooster)

case_mses <- foreach(i = 1:n_rounds, .combine = cbind) %dopar% { 
  
#case_mses <- bind_cols(lapply(1:n_rounds, function(i){ #nonParallel run version
  # Conduct cross validation across all candidate glms
  case_folds <- make_folds(its_df_quarts_largeCs %>% filter(!is.na(cases)), cluster_ids = its_df_quarts_largeCs %>% filter(!is.na(cases)) %>% pull(zip))
  
  # cv_test <- cross_validate(cv_fun = cv_glm, folds = case_folds, data = its_df_quarts_largeCs,
  #                           reg_form = paste0("cases ~ ", prop_formulas[[1]]))
  
  
  case_evals <- lapply(prop_formulas, function(l){
    cross_validate(cv_fun = cv_glm, folds = case_folds, data = its_df_quarts_largeCs %>% filter(!is.na(cases)),
                   reg_form = paste0("cases ~ ", l))
  })
  
  case_mses <- unlist(lapply(case_evals, function(l){
    mean(as.numeric(l$SE)[which(as.Date(l$date) >= policy_date)])
  }))

  #cat("Cases iteration", i, "complete\n")
  
  return(case_mses)
}#))

cat("Cases model evaluation complete\n")

#case_mses

case_modselect <- which.min(rowMeans(case_mses))

# case_preds <- as_tibble(case_evals[[case_modselect]][1:18])
# 
# case_mse_ts <- case_preds %>% 
#   mutate(SE = as.numeric(SE)) %>% 
#   group_by(date, vemquartile) %>% 
#   summarise(SE_mean = mean(SE),
#             SE_med = median(SE),
#             SE_q25 = quantile(SE, 0.25),
#             SE_q75 = quantile(SE, 0.75))
# 
# case_mse_ts %>% 
#   mutate(date = as.Date(date)) %>% 
#   filter(date >= policy_date) %>% 
#   ggplot(aes(x = date, y = SE_med, col = vemquartile, fill = vemquartile)) +
#   geom_ribbon(aes(ymin = SE_q25, ymax = SE_q75), alpha = 0.3) +
#   geom_line() +
#   scale_color_hpiq() +
#   scale_fill_hpiq() +
#   theme_classic() +
#   theme_here()

# Evaluate hospitalization models -------------------
# Candidate hospitalization GLMs
# Conduct cross validation across all candidate glms

hosp_mses <- foreach(i = 1:n_rounds, .combine = cbind) %dopar% { 
#hosp_mses <- bind_cols(lapply(1:n_rounds, function(i){
  hosp_folds <- make_folds(its_df_quarts_largeCs %>% filter(!is.na(hosps)), cluster_ids = its_df_quarts_largeCs %>% filter(!is.na(hosps)) %>% pull(zip))

  hosp_evals <- lapply(prop_formulas, function(l){
    cross_validate(cv_fun = cv_glm, folds = hosp_folds, data = its_df_quarts_largeCs %>% filter(!is.na(hosps)),
                   reg_form = paste0("hosps ~ ", l))
  })
  
  hosp_mses <- unlist(lapply(hosp_evals, function(l){
    mean(as.numeric(l$SE)[which(as.Date(l$date) >= policy_date)])
  }))
  
  #cat("Hospitalizations iteration", i, "complete\n")
  
  return(hosp_mses)
}#))

cat("Hospitalizations model evaluation complete\n")

hosp_modselect <- which.min(rowMeans(hosp_mses))

# hosp_preds <- as_tibble(hosp_evals[[hosp_modselect]][1:18])

# Evaluate mortality models -------------------
# Candidate mortality GLMs same as hospitalization glms
# Conduct cross validation across all candidate glms
mort_mses <- foreach(i = 1:n_rounds, .combine = cbind) %dopar% { 
#mort_mses <- bind_cols(lapply(1:n_rounds, function(i){
  
  mort_folds <- make_folds(its_df_quarts_largeCs %>% filter(!is.na(deaths)), cluster_ids = its_df_quarts_largeCs %>% filter(!is.na(deaths)) %>% pull(zip))
  
  mort_evals <- lapply(prop_formulas, function(l){
    cross_validate(cv_fun = cv_glm, folds = mort_folds, data = its_df_quarts_largeCs %>% filter(!is.na(deaths)),
                   reg_form = paste0("deaths ~ ", l))
  })
  
  mort_mses <- unlist(lapply(mort_evals, function(l){
    mean(as.numeric(l$SE)[which(as.Date(l$date) >= policy_date)])
  }))
  
  #cat("Deaths iteration", i, "complete\n")
  
  return(mort_mses)
}#))

parallel::stopCluster(clooster)

cat("Mortalities model evaluation complete\n")

mort_modselect <- which.min(rowMeans(mort_mses))

# mort_preds <- as_tibble(mort_evals[[mort_modselect]][1:18])

out_mods <- list("cases" = prop_formulas[case_modselect],
                 "hosps" = prop_formulas[hosp_modselect],
                 "deaths" = prop_formulas[mort_modselect])

saveRDS(out_mods, here::here("data/best_mods.rds"))
# saveRDS(list(case_preds, hosp_preds, mort_preds), 
#         here::here("data/best_mods_cv_results.rds"))
saveRDS(list("forms" = prop_formulas,
             "case_mses" = case_mses,
             "hosp_mses" = hosp_mses,
             "mort_mses" = mort_mses),
        here::here("data/mods_mean_mses.rds"))
