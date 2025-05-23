# CEM: Poverty

# Loading libraries
library(haven)
library(dplyr)
library(cem)
library(tibble)

# Loading and filtering data first
provida_df <- read_dta("/Users/homi/Library/CloudStorage/Box-Box/PROVIDA/01_data/02_build/Final_PROVIDA_data.dta") %>%
  filter(!is.na(ipbcie_stand)) %>%
  mutate(col_21 = as.factor(col_21)) %>%
  bind_cols(model.matrix(~ col_21 - 1, data = .)) %>%
  dplyr::select(-col_21)

# Define Double Lasso covariates
dl_vars <- c("age_ben", "hh_men_min_21", "hh_women_min_21", "housing_index_stand", grep("^col_21", names(provida_df), value = TRUE))

# Define All covariates
exclude <- c("ipbcie_stand", "rent_charge", "dec_extent", "treat", "cem_weight", "cem_group")
all_vars <- setdiff(names(provida_df), exclude)
all_vars <- all_vars[sapply(provida_df[all_vars],
                            function(x) (is.numeric(x) | is.factor(x)) && all(is.finite(x)) && !any(is.na(x)))]

# Run CEM function
run_cem <- function(data, outcome, match_vars, label) {
  cem_result <- cem("treat", data = data %>% 
                      dplyr::select(treat, all_of(match_vars)),
                    keep.all = TRUE)
  
  data <- data %>% 
    mutate(cem_weight = cem_result$w)
  matched <- data %>% 
    filter(!is.na(cem_weight) & cem_weight > 0)
  
  att <- with(matched, 
              weighted.mean(get(outcome)[treat == 1],
                            cem_weight[treat == 1]) -
                weighted.mean(get(outcome)[treat == 0],
                              cem_weight[treat == 0])
              )
  
  weighted_se <- function(y, w) {
    wm <- weighted.mean(y, w)
    sqrt(sum(w * (y - wm)^2) / sum(w) / length(y))
  }
  
  se <- with(matched, {
    sqrt(weighted_se(get(outcome)[treat == 1], cem_weight[treat == 1])^2 +
           weighted_se(get(outcome)[treat == 0], cem_weight[treat == 0])^2)
  })
  
  p_val <- 2 * (1 - pnorm(abs(att / se)))
  stars <- case_when(p_val < 0.001 ~ "***", p_val < 0.01 ~ "**", p_val < 0.05 ~ "*", TRUE ~ "")
  
  tibble(
    Method = label,
    Outcome = outcome,
    ATT = round(att, 3),
    SE = round(se, 3),
    P_Value = round(p_val, 3),
    Stars = stars,
    Matched_Treated = sum(cem_result$matched & data$treat == 1),
    Matched_Control = sum(cem_result$matched & data$treat == 0),
    Total_Matched = sum(cem_result$matched)
  )
}

# Results
dl_result  <- run_cem(provida_df, "ipbcie_stand", dl_vars, "CEM (DL Covariates)")
all_result <- run_cem(provida_df, "ipbcie_stand", all_vars, "CEM (All Covariates)")
bind_rows(dl_result, all_result)
