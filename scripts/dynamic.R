library(tidyverse)
library(lmtp)
library(glue)
library(future)

devtools::load_all("odtr")

# source("scripts/SL.lightgbm.R")

combined <- readRDS("data/drv/clean_weeks_with_relapse_wide_080922.rds")
imputed <- readRDS("data/drv/clean_patients_imputed_080922.rds")

dat <- imputed$data
bup <- filter(combined, medicine == "bup")

A <- glue("wk{2:11}.dose_increase_this_week")

# Was there previous time opioid use? 
condA <- bup[, glue("wk{2:11}.use_this_week")] == 1
colnames(condA) <- A

# Was dose under the allowable maximum? 
condB <- bup[, glue("wk{2:11}.dose_this_week")] < 32
colnames(condB) <- A

demog <- c("sex", "age", "xrace")

comorbidities <- c(
    "hwithdraw",
    "alcdisorder",
    "cocdisorder",
    "hasBrainDamage",
    "hasEpilepsy",
    "hasSchiz",
    "hasBipolar",
    "hasAnxPan",
    "hasMajorDep",
    "bamphetamine30_base",
    "bcannabis30_base",
    "bbenzo30_base",
    "ivdrug"
)

process_missing <- function(data) {
    as.data.frame(cbind(who = data$who, sl3::make_sl3_Task(data, names(data)[-1])$data))
}

dat <- 
    filter(dat, medicine == "bup") |> 
    select(who, project, all_of(demog), all_of(comorbidities)) |> 
    mutate(sex = if_else(sex == "female", 1, 0)) |> 
    process_missing()

W <- names(dat)[-1]
L <- lapply(2:11, \(x) c(glue("wk{x-1}.dose_this_week"), glue("wk{x}.use_this_week")))
Y <- glue("wk{3:12}.relapse_this_week")

sl <- c("SL.glm", "SL.lightgbm", "SL.earth", "SL.mean")

dat <- left_join(bup, dat)

baseline <- as.data.frame(model.matrix(reformulate(W), dat)[, -1])
W <- names(baseline)

dat <- cbind(baseline, dat[, c(A, unlist(L), Y)])

dynamic <- dat
dynamic[, A] <- apply(condA & condB, 2, \(x) as.numeric(x), simplify = FALSE)
  
plan(multisession, workers = 10) 

estims <- lmtp_sdr(
    dat, 
    A, Y, W, L,
    shifted = dynamic, 
    outcome_type = "survival", 
    folds = 10, 
    learners_outcome = sl, 
    learners_trt = sl
)

constant <- lmtp_sdr(
    dat, 
    A, Y, W, L,
    shift = static_binary_off, 
    outcome_type = "survival", 
    folds = 10, 
    learners_outcome = sl, 
    learners_trt = sl
)

plan(sequential)

saveRDS(estims, "data/drv/survival-combined-dynamic-030823.rds")
saveRDS(constant, "data/drv/survival-combined-constant-030823.rds")
