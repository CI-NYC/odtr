library(tidyverse)
library(lmtp)

# check to see if their relapse dates are before their doses begin
# these people would likely need to be removed

xbot <- readRDS("data/drv/xbot.rds")

all_ids <- xbot$baseline$PATID

ids <- setdiff(all_ids, filter(xbot$dose, is.na(first_dose_visit))$PATID)                     # removing patients with no doses
ids <- setdiff(ids, filter(xbot$relapse, week_of_relapse < xbot$dose$first_dose_visit)$PATID) # removing patients who relapse before receiving a dose

xbot <- map(xbot, \(data) filter(data, PATID %in% ids))

xbot$dose[is.na(xbot$dose)] <- 0

# Is this correct? 
xbot$use[is.na(xbot$use)] <- 1

dose_increase <- 
    pivot_longer(xbot$dose, starts_with("dose.")) |> 
    mutate(dose_increase = as.numeric(value > lag(value, n = 1, default = 0))) |> 
    select(PATID, name, dose_increase) |> 
    pivot_wider(names_from = name, values_from = dose_increase)

relapse_weeks <- c(3, 4, 6, 8, 10, 12, 14, 16, 20)

tmp <- left_join(xbot$baseline, xbot$relapse) |> 
    left_join(dose_increase) |> 
    left_join(xbot$use)

A <- glue::glue("dose.{relapse_weeks[1:(length(relapse_weeks) - 1)]}")
Y <- glue::glue("relapse.{relapse_weeks[2:length(relapse_weeks)]}")
W <- c("age", "hispanic")
L <- lapply(relapse_weeks[1:(length(relapse_weeks) - 1)], 
            \(x) c(glue::glue("anyopioiduse.{x}")))

lmtp_sdr(
    tmp, 
    A, Y, W, L,
    shift = static_binary_on, 
    outcome_type = "survival", 
    folds = 1, 
    learners_outcome = sl3::Lrnr_xgboost$new(), 
    learners_trt = sl3::Lrnr_xgboost$new()
)
