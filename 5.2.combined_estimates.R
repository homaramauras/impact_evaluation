# ------------------------------------------------------------- #
# Combined matching estimates ----
# ------------------------------------------------------------- #

# libraries
library(haven)
library(dplyr)
library(cem)
library(ggplot2)
library(tibble)
library(patchwork)
library(cobalt)
library(MatchIt)
library(boot)
library(tidyr)
library(Matching)
library(rgenoud)
library(tibble)
library(knitr)
library(gt)
library(officer)
library(flextable)


# loading data
provida_df <- read_dta("/Users/homi/Library/CloudStorage/Box-Box/PROVIDA/01_data/02_build/Final_PROVIDA_data.dta")

# creating col_21 dummies
provida_df <- provida_df %>%
  mutate(col_21 = as.factor(col_21)) %>%
  bind_cols(model.matrix(~ col_21 - 1, data = .)) %>%
  dplyr::select(-col_21)

# drop missing values for all outcomes
provida_df <- provida_df %>%
  filter(!is.na(ipbcie_stand), !is.na(rent_charge), !is.na(dec_extent))

# helper function
format_att_result <- function(method, outcome, att, se, p_val, stars, mt, mc) {
  tibble(
    Method = method,
    Outcome = outcome,
    ATT = att,
    SE = se,
    P_Value = p_val,
    Stars = stars,
    Matched_Treated = mt,
    Matched_Control = mc,
    Total_Matched = mt + mc
  )
}

# ------------------------------------------------------------- #
# CEM function execution ----
# ------------------------------------------------------------- #

run_cem <- function(data, outcome, match_vars, label) {
  cem_df <- data %>% dplyr::select(treat, all_of(match_vars))
  
  cem_result <- cem(
    treatment = "treat",
    data = cem_df,
    keep.all = TRUE
  )
  
  data <- data %>%
    mutate(
      cem_weight = cem_result$w,
      cem_group  = cem_result$group
    )
  
  matched_data <- data %>% filter(!is.na(cem_weight) & cem_weight > 0)
  
  att <- with(matched_data, 
              weighted.mean(get(outcome)[treat == 1], cem_weight[treat == 1]) -
                weighted.mean(get(outcome)[treat == 0], cem_weight[treat == 0]))
  
  weighted_se <- function(y, weights) {
    wm <- weighted.mean(y, weights)
    variance <- sum(weights * (y - wm)^2) / sum(weights)
    sqrt(variance / length(y))
  }
  
  se <- with(matched_data, {
    se_treated <- weighted_se(get(outcome)[treat == 1], cem_weight[treat == 1])
    se_control <- weighted_se(get(outcome)[treat == 0], cem_weight[treat == 0])
    sqrt(se_treated^2 + se_control^2)
  })
  
  ci_lower <- att - 1.96 * se
  ci_upper <- att + 1.96 * se
  t_stat <- att / se
  p_val <- 2 * (1 - pnorm(abs(t_stat)))
  
  stars <- case_when(
    p_val < 0.001 ~ "***",
    p_val < 0.01  ~ "**",
    p_val < 0.05  ~ "*",
    TRUE          ~ ""
  )
  
  matched_treated  <- sum(cem_result$matched & cem_df$treat == 1)
  matched_control  <- sum(cem_result$matched & cem_df$treat == 0)
  total_matched    <- sum(cem_result$matched)
  
  result_table <- tibble(
    Method = label,
    Outcome = outcome,
    ATT = att,
    SE = se,
    P_Value = p_val,
    Stars = stars,
    Matched_Treated = matched_treated,
    Matched_Control = matched_control,
    Total_Matched = total_matched
  )
  
  return(list(
    results = result_table,
    matched_data = matched_data
  ))
}

# complete cases
provida_df <- provida_df %>% 
  filter(!is.na(ipbcie_stand), !is.na(rent_charge), !is.na(dec_extent))

# colonies fixed effects
col21_dummies <- names(provida_df)[grepl("^col_21", names(provida_df))]

# DL Covariates
match_vars_ipbcie <- c("age_ben", "hh_men_min_21", "hh_women_min_21", "housing_index_stand", col21_dummies)
match_vars_rent    <- c("gender", "hh_men_min_21", "housing_index_stand", col21_dummies)
match_vars_dec     <- c("age_ben", "hh_men_21", "hh_men_min_21", "hh_women_min_21", "housing_index_stand", col21_dummies)

ipbcie_dl_out <- run_cem(provida_df, "ipbcie_stand", match_vars_ipbcie, "CEM (DL Covariates)")
rent_dl_out   <- run_cem(provida_df, "rent_charge", match_vars_rent, "CEM (DL Covariates)")
dec_dl_out    <- run_cem(provida_df, "dec_extent", match_vars_dec, "CEM (DL Covariates)")

# CEM all covariates
outcomes <- c("ipbcie_stand", "rent_charge", "dec_extent")
all_covariates <- setdiff(names(provida_df), c("treat", outcomes, "cem_weight", "cem_group"))

clean_covariates <- all_covariates[ sapply(provida_df[all_covariates], function(x) {
  is_numeric_or_factor <- is.numeric(x) | is.factor(x)
  has_valid_data <- any(!is.na(x)) & all(is.finite(x[!is.na(x)]))
  return(is_numeric_or_factor & has_valid_data)
})]

ipbcie_all_out <- run_cem(provida_df, "ipbcie_stand", clean_covariates, "CEM (All Covariates)")
rent_all_out   <- run_cem(provida_df, "rent_charge", clean_covariates, "CEM (All Covariates)")
dec_all_out    <- run_cem(provida_df, "dec_extent", clean_covariates, "CEM (All Covariates)")

cem_results <- bind_rows(
  ipbcie_dl_out$results,
  rent_dl_out$results,
  dec_dl_out$results,
  ipbcie_all_out$results,
  rent_all_out$results,
  dec_all_out$results
)

# ------------------------------------------------------------- #
# MDM function execution ----
# ------------------------------------------------------------- #

run_mdm <- function(data, outcome, formula, label) {
  matched <- matchit(formula, data = data, method = "nearest", distance = "mahalanobis")
  matched_df <- match.data(matched, data = data)
  
  att <- with(matched_df, weighted.mean(get(outcome)[treat == 1], weights[treat == 1]) -
                weighted.mean(get(outcome)[treat == 0], weights[treat == 0]))
  
  att_fun <- function(data, indices) {
    d <- data[indices, ]
    wm1 <- weighted.mean(d[[outcome]][d$treat == 1], d$weights[d$treat == 1])
    wm0 <- weighted.mean(d[[outcome]][d$treat == 0], d$weights[d$treat == 0])
    return(wm1 - wm0)
  }
  
  boot_out <- boot(data = matched_df, statistic = att_fun, R = 1000)
  se <- sd(boot_out$t)
  
  t_stat <- att / se
  p_val <- 2 * (1 - pnorm(abs(t_stat)))
  stars <- case_when(
    p_val < 0.001 ~ "***",
    p_val < 0.01  ~ "**",
    p_val < 0.05  ~ "*",
    TRUE ~ ""
  )
  
  mt <- sum(matched_df$treat == 1)
  mc <- sum(matched_df$treat == 0)
  
  return(format_att_result(label, outcome, att, se, p_val, stars, mt, mc))
}

# Clean up names to make them safe for formulas
names(provida_df) <- make.names(names(provida_df))

# Recreate dummy name list after sanitizing
col21_dummies <- names(provida_df)[grepl("^col_21", names(provida_df))]

# IPBCIE outcome (poverty)
form_ipbcie <- reformulate(
  c("age_ben", "hh_men_min_21", "hh_women_min_21", "housing_index_stand", col21_dummies),
  response = "treat"
)

# RENT outcome
form_rent <- reformulate(
  c("gender", "hh_men_21", "housing_index_stand", col21_dummies),
  response = "treat"
)

# DEC_EXTENT outcome (empowerment)
form_dec <- reformulate(
  c("age_ben", "hh_men_21", "hh_men_min_21", "hh_women_min_21", "housing_index_stand", col21_dummies),
  response = "treat"
)


# Run MDM
mdm_results <- bind_rows(
  run_mdm(provida_df, "ipbcie_stand", form_ipbcie, "MDM (DL Covariates)"),
  run_mdm(provida_df, "rent_charge", form_rent, "MDM (DL Covariates)"),
  run_mdm(provida_df, "dec_extent", form_dec, "MDM (DL Covariates)")
)

# all covariates
### define list of all usable covariates
outcomes <- c("ipbcie_stand", "rent_charge", "dec_extent", "cem_weight", "cem_group")
all_covariates <- setdiff(names(provida_df), c("treat", outcomes))

# get rid of variables that have *any* NA or non-finite values
clean_covariates <- all_covariates[
  sapply(provida_df[all_covariates], function(x) {
    is.numeric(x) || is.factor(x)
  }) &
    colSums(is.na(provida_df[all_covariates])) == 0 &
    colSums(!is.finite(as.matrix(provida_df[all_covariates]))) == 0
]


form_all <- reformulate(clean_covariates, response = "treat")

mdm_results_all <- bind_rows(
  run_mdm(provida_df, "ipbcie_stand", form_all, "MDM (All Covariates)"),
  run_mdm(provida_df, "rent_charge", form_all, "MDM (All Covariates)"),
  run_mdm(provida_df, "dec_extent", form_all, "MDM (All Covariates)")
)


mdm_results <- bind_rows(mdm_results, mdm_results_all)

# ------------------------------------------------------------- #
# PSM function execution ----
# ------------------------------------------------------------- #

# Helper: Bootstrap ATT with IPW
get_psm_att <- function(data, outcome_var, weight_var = "ipw_att", label) {
  att_fun <- function(data, indices) {
    d <- data[indices, ]
    treated <- weighted.mean(d[[outcome_var]][d$treat == 1], d[[weight_var]][d$treat == 1])
    control <- weighted.mean(d[[outcome_var]][d$treat == 0], d[[weight_var]][d$treat == 0])
    return(treated - control)
  }
  
  boot_out <- boot(data = data, statistic = att_fun, R = 1000)
  att <- boot_out$t0
  se <- sd(boot_out$t)
  t_stat <- att / se
  p_val <- 2 * (1 - pnorm(abs(t_stat)))
  stars <- case_when(
    p_val < 0.001 ~ "***",
    p_val < 0.01 ~ "**",
    p_val < 0.05 ~ "*",
    TRUE ~ ""
  )
  
  tibble(
    Method = label,
    Outcome = outcome_var,
    ATT = att,
    SE = se,
    P_Value = p_val,
    Stars = stars,
    Matched_Treated = sum(data$treat == 1),
    Matched_Control = sum(data$treat == 0),
    Total_Matched = nrow(data)
  )
}

# ------------------------------------------------------------- #
# Generate PSM results for all 3 outcomes
# ------------------------------------------------------------- #

# POVERTY
provida_pov_dl <- provida_df %>%
  dplyr::select(treat, ipbcie_stand, age_ben, hh_men_min_21, hh_women_min_21,
         housing_index_stand, starts_with("col_21"))
form_pov_dl <- reformulate(setdiff(names(provida_pov_dl), c("treat", "ipbcie_stand")), response = "treat")
model_pov_dl <- glm(form_pov_dl, data = provida_pov_dl, family = binomial)
provida_pov_dl$pscore <- predict(model_pov_dl, type = "response")
provida_pov_dl <- provida_pov_dl %>%
  mutate(ipw_att = ifelse(treat == 1, 1, pscore / (1 - pscore)))
psm_pov_dl <- get_psm_att(provida_pov_dl, "ipbcie_stand", "ipw_att", "PSM (DL Covariates)")

provida_pov_all <- provida_df %>%
  dplyr::select(treat, ipbcie_stand, gender, age_ben, dep_21, mun_21, dis_21,
         hh_num_21, hh_men_21, hh_men_min_21, hh_women_21, hh_women_min_21,
         housing_index_stand, starts_with("col_21"))
form_pov_all <- reformulate(setdiff(names(provida_pov_all), c("treat", "ipbcie_stand")), response = "treat")
model_pov_all <- glm(form_pov_all, data = provida_pov_all, family = binomial)
provida_pov_all$pscore <- predict(model_pov_all, type = "response")
provida_pov_all <- provida_pov_all %>%
  mutate(ipw_att = ifelse(treat == 1, 1, pscore / (1 - pscore)))
psm_pov_all <- get_psm_att(provida_pov_all, "ipbcie_stand", "ipw_att", "PSM (All Covariates)")

# RENT
provida_rent_dl <- provida_df %>%
  dplyr::select(treat, rent_charge, gender, hh_men_min_21, housing_index_stand, starts_with("col_21"))
form_rent_dl <- reformulate(setdiff(names(provida_rent_dl), c("treat", "rent_charge")), response = "treat")
model_rent_dl <- glm(form_rent_dl, data = provida_rent_dl, family = binomial)
provida_rent_dl$pscore <- predict(model_rent_dl, type = "response")
provida_rent_dl <- provida_rent_dl %>%
  mutate(ipw_att = ifelse(treat == 1, 1, pscore / (1 - pscore)))
psm_rent_dl <- get_psm_att(provida_rent_dl, "rent_charge", "ipw_att", "PSM (DL Covariates)")

provida_rent_all <- provida_df %>%
  dplyr::select(treat, rent_charge, gender, age_ben, dep_21, mun_21, dis_21,
         hh_num_21, hh_men_21, hh_men_min_21, hh_women_21, hh_women_min_21,
         housing_index_stand, starts_with("col_21"))
form_rent_all <- reformulate(setdiff(names(provida_rent_all), c("treat", "rent_charge")), response = "treat")
model_rent_all <- glm(form_rent_all, data = provida_rent_all, family = binomial)
provida_rent_all$pscore <- predict(model_rent_all, type = "response")
provida_rent_all <- provida_rent_all %>%
  mutate(ipw_att = ifelse(treat == 1, 1, pscore / (1 - pscore)))
psm_rent_all <- get_psm_att(provida_rent_all, "rent_charge", "ipw_att", "PSM (All Covariates)")

# EMPOWER
provida_emp_dl <- provida_df %>%
  dplyr::select(treat, dec_extent, age_ben, hh_men_21, hh_men_min_21, hh_women_min_21,
         housing_index_stand, starts_with("col_21"))
form_emp_dl <- reformulate(setdiff(names(provida_emp_dl), c("treat", "dec_extent")), response = "treat")
model_emp_dl <- glm(form_emp_dl, data = provida_emp_dl, family = binomial)
provida_emp_dl$pscore <- predict(model_emp_dl, type = "response")
provida_emp_dl <- provida_emp_dl %>%
  mutate(ipw_att = ifelse(treat == 1, 1, pscore / (1 - pscore)))
psm_emp_dl <- get_psm_att(provida_emp_dl, "dec_extent", "ipw_att", "PSM (DL Covariates)")

provida_emp_all <- provida_df %>%
  dplyr::select(treat, dec_extent, gender, age_ben, dep_21, mun_21, dis_21,
         hh_num_21, hh_men_21, hh_men_min_21, hh_women_21, hh_women_min_21,
         housing_index_stand, starts_with("col_21"))
form_emp_all <- reformulate(setdiff(names(provida_emp_all), c("treat", "dec_extent")), response = "treat")
model_emp_all <- glm(form_emp_all, data = provida_emp_all, family = binomial)
provida_emp_all$pscore <- predict(model_emp_all, type = "response")
provida_emp_all <- provida_emp_all %>%
  mutate(ipw_att = ifelse(treat == 1, 1, pscore / (1 - pscore)))
psm_emp_all <- get_psm_att(provida_emp_all, "dec_extent", "ipw_att", "PSM (All Covariates)")

# psm results
psm_results <- bind_rows(
  psm_pov_dl, psm_pov_all,
  psm_rent_dl, psm_rent_all,
  psm_emp_dl, psm_emp_all
)

# ------------------------------------------------------------- #
# Genetic matching function execution ----
# ------------------------------------------------------------- #

# Load data
provida_df <- read_dta("/Users/homi/Library/CloudStorage/Box-Box/PROVIDA/01_data/02_build/Final_PROVIDA_data.dta") %>%
  mutate(treat = as.numeric(treat)) %>%
  filter(treat %in% c(0,1)) # ensure binary treatment

# Clean up col_21 dummy names to be syntactically safe
col21_dummies <- model.matrix(~ col_21 - 1, data = provida_df)
colnames(col21_dummies) <- make.names(colnames(col21_dummies))

# Append dummies to main data
provida_df <- cbind(provida_df, col21_dummies)

# Define DL covariates
cov_dl_ipbcie <- c("age_ben", "hh_men_min_21", "hh_women_min_21", "housing_index_stand", colnames(col21_dummies))
cov_dl_rent   <- c("gender", "hh_men_min_21", "housing_index_stand", colnames(col21_dummies))
cov_dl_dec    <- c("age_ben", "hh_men_21", "hh_men_min_21", "hh_women_min_21", "housing_index_stand", colnames(col21_dummies))

# Define all covariates
outcomes <- c("ipbcie_stand", "rent_charge", "dec_extent", "cem_weight", "cem_group")
cov_all <- setdiff(names(provida_df), c("treat", outcomes))
cov_all <- cov_all[sapply(provida_df[cov_all], function(x) is.numeric(x) && all(is.finite(x)))]

# ATT function for bootstrapping
att_fun <- function(outcome) {
  function(data, indices) {
    d <- data[indices, ]
    treated <- weighted.mean(d[[outcome]][d$treat == 1], d$weights[d$treat == 1])
    control <- weighted.mean(d[[outcome]][d$treat == 0], d$weights[d$treat == 0])
    return(treated - control)
  }
}


# Wrapper to run GenMatch
run_genmatch_att <- function(data, outcome, formula, label) {
  m.out <- matchit(formula, data = data, method = "genetic", pop.size = 10)
  matched <- match.data(m.out)
  
  att <- with(matched, weighted.mean(get(outcome)[treat == 1], weights[treat == 1]) -
                weighted.mean(get(outcome)[treat == 0], weights[treat == 0]))
  
  boot_out <- boot(data = matched, statistic = att_fun(outcome), R = 1000)
  se <- sd(boot_out$t)
  p_val <- 2 * (1 - pnorm(abs(att / se)))
  
  stars <- case_when(
    p_val < 0.001 ~ "***",
    p_val < 0.01  ~ "**",
    p_val < 0.05  ~ "*",
    TRUE ~ ""
  )
  
  tibble(
    Method = label,
    Outcome = outcome,
    ATT = att,
    SE = se,
    P_Value = p_val,
    Stars = stars,
    Matched_Treated = sum(matched$treat == 1),
    Matched_Control = sum(matched$treat == 0),
    Total_Matched = nrow(matched)
  )
}

# Run for DL Covariates
gen_dl_results <- bind_rows(
  run_genmatch_att(provida_df %>% drop_na(ipbcie_stand), "ipbcie_stand", as.formula(paste("treat ~", paste(cov_dl_ipbcie, collapse = "+"))), "GenMatch (DL Covariates)"),
  run_genmatch_att(provida_df %>% drop_na(rent_charge), "rent_charge", as.formula(paste("treat ~", paste(cov_dl_rent, collapse = "+"))), "GenMatch (DL Covariates)"),
  run_genmatch_att(provida_df %>% drop_na(dec_extent), "dec_extent", as.formula(paste("treat ~", paste(cov_dl_dec, collapse = "+"))), "GenMatch (DL Covariates)")
)

# Run for All Covariates
form_all <- as.formula(paste("treat ~", paste(cov_all, collapse = "+")))

gen_all_results <- bind_rows(
  run_genmatch_att(provida_df %>% drop_na(ipbcie_stand), "ipbcie_stand", form_all, "GenMatch (All Covariates)"),
  run_genmatch_att(provida_df %>% drop_na(rent_charge), "rent_charge", form_all, "GenMatch (All Covariates)"),
  run_genmatch_att(provida_df %>% drop_na(dec_extent), "dec_extent", form_all, "GenMatch (All Covariates)")
)

# Combine everything
genmatch_results <- bind_rows(gen_dl_results, gen_all_results)

# ------------------------------------------------------------- #
# OLS Results Tables ----
# ------------------------------------------------------------- #

# Poverty models
model_pov_bivar <- lm(ipbcie_stand ~ treat, data = provida_df)
model_pov_multiv <- lm(ipbcie_stand ~ treat + gender + age_ben + hh_num_21 + hh_men_21 +
                         hh_men_min_21 + hh_women_21 + hh_women_min_21 + housing_index_stand +
                         ., data = provida_df[, c("ipbcie_stand", "treat", "gender", "age_ben", "hh_num_21", 
                                                  "hh_men_21", "hh_men_min_21", "hh_women_21", 
                                                  "hh_women_min_21", "housing_index_stand", 
                                                  grep("^col_21", names(provida_df), value = TRUE))])
model_pov_dl <- lm(ipbcie_stand ~ treat + age_ben + hh_men_min_21 + hh_women_min_21 +
                     housing_index_stand + ., data = provida_df[, c("ipbcie_stand", "treat", "age_ben",
                                                                    "hh_men_min_21", "hh_women_min_21",
                                                                    "housing_index_stand", 
                                                                    grep("^col_21", names(provida_df), value = TRUE))])

# Housing Quality models
model_rent_bivar <- lm(rent_charge ~ treat, data = provida_df)
model_rent_multiv <- lm(rent_charge ~ treat + gender + age_ben + hh_num_21 + hh_men_21 +
                          hh_men_min_21 + hh_women_21 + hh_women_min_21 + housing_index_stand +
                          ., data = provida_df[, c("rent_charge", "treat", "gender", "age_ben", "hh_num_21", 
                                                   "hh_men_21", "hh_men_min_21", "hh_women_21", 
                                                   "hh_women_min_21", "housing_index_stand", 
                                                   grep("^col_21", names(provida_df), value = TRUE))])
model_rent_dl <- lm(rent_charge ~ treat + gender + hh_men_21 + housing_index_stand + ., 
                    data = provida_df[, c("rent_charge", "treat", "gender", "hh_men_21",
                                          "housing_index_stand", 
                                          grep("^col_21", names(provida_df), value = TRUE))])

# Empowerment models
model_emp_bivar <- lm(dec_extent ~ treat, data = provida_df)
model_emp_multiv <- lm(dec_extent ~ treat + gender + age_ben + hh_num_21 + hh_men_21 +
                         hh_men_min_21 + hh_women_21 + hh_women_min_21 + housing_index_stand +
                         ., data = provida_df[, c("dec_extent", "treat", "gender", "age_ben", "hh_num_21", 
                                                  "hh_men_21", "hh_men_min_21", "hh_women_21", 
                                                  "hh_women_min_21", "housing_index_stand", 
                                                  grep("^col_21", names(provida_df), value = TRUE))])
model_emp_dl <- lm(dec_extent ~ treat + age_ben + hh_men_21 + hh_men_min_21 + 
                     hh_women_min_21 + housing_index_stand + ., 
                   data = provida_df[, c("dec_extent", "treat", "age_ben", "hh_men_21", 
                                         "hh_men_min_21", "hh_women_min_21", 
                                         "housing_index_stand", 
                                         grep("^col_21", names(provida_df), value = TRUE))])

# Combine all OLS results
ols_results <- bind_rows(
  broom::tidy(model_pov_bivar) %>%
    filter(term == "treat") %>%
    mutate(Method = "OLS – Bivariate", Outcome = "Poverty", Observations = nobs(model_pov_bivar)),
  broom::tidy(model_pov_multiv) %>%
    filter(term == "treat") %>%
    mutate(Method = "OLS – Multivariate", Outcome = "Poverty", Observations = nobs(model_pov_multiv)),
  broom::tidy(model_pov_dl) %>%
    filter(term == "treat") %>%
    mutate(Method = "OLS – DL Covariates", Outcome = "Poverty", Observations = nobs(model_pov_dl)),
  
  broom::tidy(model_rent_bivar) %>%
    filter(term == "treat") %>%
    mutate(Method = "OLS – Bivariate", Outcome = "Housing Quality", Observations = nobs(model_rent_bivar)),
  broom::tidy(model_rent_multiv) %>%
    filter(term == "treat") %>%
    mutate(Method = "OLS – Multivariate", Outcome = "Housing Quality", Observations = nobs(model_rent_multiv)),
  broom::tidy(model_rent_dl) %>%
    filter(term == "treat") %>%
    mutate(Method = "OLS – DL Covariates", Outcome = "Housing Quality", Observations = nobs(model_rent_dl)),
  
  broom::tidy(model_emp_bivar) %>%
    filter(term == "treat") %>%
    mutate(Method = "OLS – Bivariate", Outcome = "Empowerment", Observations = nobs(model_emp_bivar)),
  broom::tidy(model_emp_multiv) %>%
    filter(term == "treat") %>%
    mutate(Method = "OLS – Multivariate", Outcome = "Empowerment", Observations = nobs(model_emp_multiv)),
  broom::tidy(model_emp_dl) %>%
    filter(term == "treat") %>%
    mutate(Method = "OLS – DL Covariates", Outcome = "Empowerment", Observations = nobs(model_emp_dl))
) %>%
  transmute(
    Method,
    Outcome,
    ATT = estimate,
    SE = std.error,
    P_Value = p.value,
    Stars = case_when(
      P_Value < 0.001 ~ "***",
      P_Value < 0.01  ~ "**",
      P_Value < 0.05  ~ "*",
      TRUE ~ ""
    ),
    Observations
  ) %>%
  arrange(Outcome, Method)

# ------------------------------------------------------------- #
# ATT Results Tables ----
# ------------------------------------------------------------- #

att_table <- bind_rows(
  cem_results,
  mdm_results,
  psm_results,
  genmatch_results,
  ols_results)

att_table <- att_table %>%
  mutate(Outcome = case_when(
    Outcome == "ipbcie_stand" ~ "Poverty",
    Outcome == "rent_charge"  ~ "Housing Quality",
    Outcome == "dec_extent"   ~ "Empowerment",
    TRUE ~ Outcome
  ))

att_table <- att_table %>%
  mutate(
    Covariate_Set = case_when(
      grepl("DL Covariates", Method)   ~ "DL Covariates",
      grepl("All Covariates", Method)  ~ "All Covariates",
      grepl("OLS", Method)             ~ "OLS",
      TRUE                             ~ "Other"
    )
  ) %>%
  arrange(Covariate_Set, Outcome, Method)


att_table <- att_table %>% dplyr::select(-Covariate_Set)

# table viz
### basic
print(att_table, digits = 4, n = Inf)

### markdown table
kable(att_table, digits = 4, format = "markdown")

### tibble with stars
att_table_pretty <- att_table %>%
  mutate(
    ATT_pretty = sprintf("%.4f%s", ATT, Stars),
    SE_pretty = sprintf("(%.4f)", SE)
  ) %>%
  unite("ATT (SE)", ATT_pretty, SE_pretty, sep = " ") %>%
  dplyr::select(Method, Outcome, `ATT (SE)`, Matched_Treated, Matched_Control, Total_Matched)
print(att_table_pretty, n = Inf)


### gt table
gt_table <- att_table %>%
  mutate(
    ATT = round(ATT, 2),
    SE = round(SE, 3),
    P_Value = round(P_Value, 3)
  ) %>%
  gt() %>%
  tab_header(title = "ATT Estimates by Method and Outcome") %>%
  cols_label(
    Method = "Method",
    Outcome = "Outcome Variable",
    ATT = "ATT",
    SE = "Standard Error",
    P_Value = "P-Value",
    Stars = "Significance",
    Matched_Treated = "Treated",
    Matched_Control = "Control",
    Total_Matched = "Total"
  ) %>%
  sub_missing(everything(), missing_text = "") %>%
  tab_options(table.font.size = "small") %>%
  tab_source_note(
    md("_Significance levels: *****p < 0.001**, ****p < 0.01**, ***p < 0.05**_")
  )

gtsave(gt_table, "att_table.png")


# ------------------------------------------------------------- #
# GT Table Generation with Reordered Method & Outcome ----
# ------------------------------------------------------------- #

# Define desired orderings
outcome_order <- c("Poverty", "Housing Quality", "Empowerment")
method_order <- c("PSM", "MDM", "CEM", "GenMatch")

# Create subset: DL Covariates
att_dl <- att_table %>%
  filter(grepl("DL Covariates", Method)) %>%
  mutate(
    Outcome = factor(Outcome, levels = outcome_order),
    Method = factor(case_when(
      grepl("PSM", Method) ~ "PSM",
      grepl("MDM", Method) ~ "MDM",
      grepl("CEM", Method) ~ "CEM",
      grepl("GenMatch", Method) ~ "GenMatch",
      TRUE ~ Method
    ), levels = method_order)
  ) %>%
  arrange(Outcome, Method)

# Create subset: All Covariates
att_all <- att_table %>%
  filter(grepl("All Covariates", Method)) %>%
  mutate(
    Outcome = factor(Outcome, levels = outcome_order),
    Method = factor(case_when(
      grepl("PSM", Method) ~ "PSM",
      grepl("MDM", Method) ~ "MDM",
      grepl("CEM", Method) ~ "CEM",
      grepl("GenMatch", Method) ~ "GenMatch",
      TRUE ~ Method
    ), levels = method_order)
  ) %>%
  arrange(Outcome, Method)

# subset: OLS
### OLS with number of observations
# ------------------------------------------------------------- #
# OLS Results Tables (Bivariate, Multivariate, DL) ----
# ------------------------------------------------------------- #

ols_results <- bind_rows(
  # Poverty
  broom::tidy(lm(ipbcie_stand ~ treat, data = provida_df)) %>%
    filter(term == "treat") %>%
    mutate(Method = "OLS – Bivariate", Outcome = "Poverty", Observations = nobs(lm(ipbcie_stand ~ treat, data = provida_df))),
  
  broom::tidy(lm(ipbcie_stand ~ treat + gender + age_ben + hh_num_21 + hh_men_21 +
                   hh_men_min_21 + hh_women_21 + hh_women_min_21 + housing_index_stand +
                   ., data = provida_df[, c("ipbcie_stand", "treat", "gender", "age_ben", "hh_num_21", 
                                            "hh_men_21", "hh_men_min_21", "hh_women_21", 
                                            "hh_women_min_21", "housing_index_stand", 
                                            grep("^col_21", names(provida_df), value = TRUE))])) %>%
    filter(term == "treat") %>%
    mutate(Method = "OLS – Multivariate", Outcome = "Poverty", Observations = nobs(.)),
  
  broom::tidy(lm(ipbcie_stand ~ treat + age_ben + hh_men_min_21 + hh_women_min_21 +
                   housing_index_stand + ., data = provida_df[, c("ipbcie_stand", "treat", "age_ben",
                                                                  "hh_men_min_21", "hh_women_min_21",
                                                                  "housing_index_stand", 
                                                                  grep("^col_21", names(provida_df), value = TRUE))])) %>%
    filter(term == "treat") %>%
    mutate(Method = "OLS – DL Covariates", Outcome = "Poverty", Observations = nobs(.)),
  
  # Housing Quality
  broom::tidy(lm(rent_charge ~ treat, data = provida_df)) %>%
    filter(term == "treat") %>%
    mutate(Method = "OLS – Bivariate", Outcome = "Housing Quality", Observations = nobs(lm(rent_charge ~ treat, data = provida_df))),
  
  broom::tidy(lm(rent_charge ~ treat + gender + age_ben + hh_num_21 + hh_men_21 +
                   hh_men_min_21 + hh_women_21 + hh_women_min_21 + housing_index_stand +
                   ., data = provida_df[, c("rent_charge", "treat", "gender", "age_ben", "hh_num_21", 
                                            "hh_men_21", "hh_men_min_21", "hh_women_21", 
                                            "hh_women_min_21", "housing_index_stand", 
                                            grep("^col_21", names(provida_df), value = TRUE))])) %>%
    filter(term == "treat") %>%
    mutate(Method = "OLS – Multivariate", Outcome = "Housing Quality", Observations = nobs(.)),
  
  broom::tidy(lm(rent_charge ~ treat + gender + hh_men_21 + housing_index_stand +
                   ., data = provida_df[, c("rent_charge", "treat", "gender", "hh_men_21",
                                            "housing_index_stand", 
                                            grep("^col_21", names(provida_df), value = TRUE))])) %>%
    filter(term == "treat") %>%
    mutate(Method = "OLS – DL Covariates", Outcome = "Housing Quality", Observations = nobs(.)),
  
  # Empowerment
  broom::tidy(lm(dec_extent ~ treat, data = provida_df)) %>%
    filter(term == "treat") %>%
    mutate(Method = "OLS – Bivariate", Outcome = "Empowerment", Observations = nobs(lm(dec_extent ~ treat, data = provida_df))),
  
  broom::tidy(lm(dec_extent ~ treat + gender + age_ben + hh_num_21 + hh_men_21 +
                   hh_men_min_21 + hh_women_21 + hh_women_min_21 + housing_index_stand +
                   ., data = provida_df[, c("dec_extent", "treat", "gender", "age_ben", "hh_num_21", 
                                            "hh_men_21", "hh_men_min_21", "hh_women_21", 
                                            "hh_women_min_21", "housing_index_stand", 
                                            grep("^col_21", names(provida_df), value = TRUE))])) %>%
    filter(term == "treat") %>%
    mutate(Method = "OLS – Multivariate", Outcome = "Empowerment", Observations = nobs(.)),
  
  broom::tidy(lm(dec_extent ~ treat + age_ben + hh_men_21 + hh_men_min_21 + 
                   hh_women_min_21 + housing_index_stand + ., 
                 data = provida_df[, c("dec_extent", "treat", "age_ben", "hh_men_21", 
                                       "hh_men_min_21", "hh_women_min_21", 
                                       "housing_index_stand", 
                                       grep("^col_21", names(provida_df), value = TRUE))])) %>%
    filter(term == "treat") %>%
    mutate(Method = "OLS – DL Covariates", Outcome = "Empowerment", Observations = nobs(.))
) %>%
  transmute(
    Method,
    Outcome,
    ATT = estimate,
    SE = std.error,
    P_Value = p.value,
    Stars = case_when(
      P_Value < 0.001 ~ "***",
      P_Value < 0.01  ~ "**",
      P_Value < 0.05  ~ "*",
      TRUE ~ ""
    ),
    Observations = Observations
  ) %>%
  arrange(Outcome, Method)


make_att_gt <- function(df, title) {
  df %>%
    mutate(
      ATT = round(ATT, 2),
      SE = round(SE, 3),
      P_Value = round(P_Value, 3)
    ) %>%
    gt() %>%
    tab_header(title = title) %>%
    {
      if ("Observations" %in% names(df)) {
        cols_label(
          .,
          Method = "Method",
          Outcome = "Outcome Variable",
          ATT = "ATT",
          SE = "Standard Error",
          P_Value = "P-Value",
          Stars = "Significance",
          Observations = "Observations"
        )
      } else {
        cols_label(
          .,
          Method = "Method",
          Outcome = "Outcome Variable",
          ATT = "ATT",
          SE = "Standard Error",
          P_Value = "P-Value",
          Stars = "Significance",
          Matched_Treated = "Treated",
          Matched_Control = "Control",
          Total_Matched = "Total"
        )
      }
    } %>%
    tab_source_note(
      md("_Significance levels: *****p < 0.001**, ****p < 0.01**, ***p < 0.05**_")
    ) %>%
    tab_options(table.font.size = "small")
}



# Generate GT Tables
gt_dl  <- make_att_gt(att_dl,  "ATT Estimates – DL Covariates")
gt_all <- make_att_gt(att_all, "ATT Estimates – All Covariates")
gt_ols <- make_att_gt(att_ols, "ATT Estimates – OLS Regressions")

# Save
gtsave(gt_dl,  "att_dl.png")
gtsave(gt_all, "att_all.png")
gtsave(gt_ols, "att_ols.png")

# Create flextable
ft_ols <- flextable(att_ols)
ft_dl <- flextable(att_dl)
ft_all <- flextable(att_all)

# Add to Word document
doc_ols <- read_docx() %>%
  body_add_par("OLS ATT Estimates", style = "heading 1") %>%
  body_add_flextable(ft_ols)

doc_dl <- read_docx() %>%
  body_add_par("DL ATT Estimates", style = "heading 1") %>%
  body_add_flextable(ft_dl)

doc_all <- read_docx() %>%
  body_add_par("ALL ATT Estimates", style = "heading 1") %>%
  body_add_flextable(ft_all)

# Save
print(doc_ols, target = "ols_results.docx")
print(doc_dl, target = "dl_results.docx")
print(doc_all, target = "all_results.docx")





