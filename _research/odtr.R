library(tidyverse)
library(lmtp)
library(glue)
library(future)

devtools::load_all("odtr")

imputed <- readRDS("data/drv/clean_patients_imputed_080922.rds")
combined <- readRDS("data/drv/clean_weeks_with_relapse_wide_080922.rds")

dat <- imputed$data

bup <- filter(combined, medicine == "bup")

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

A <- glue("wk{2:11}.dose_increase_this_week")
W <- names(dat)[-1]
L <- lapply(2:11, \(x) c(glue("wk{x-1}.dose_this_week"), glue("wk{x}.use_this_week")))
Y <- glue("wk{3:12}.relapse_this_week")

sl <- c("SL.glm", "SL.mean", "SL.lightgbm", "SL.glmnet", "SL.earth")
# sl <- c("SL.glm", "SL.mean", "SL.lightgbm", "SL.earth")
sl <- c("SL.glm", "SL.mean", "SL.lightgbm", "SL.glmnet")

dat <- left_join(bup, dat)

baseline <- as.data.frame(model.matrix(reformulate(W), dat)[, -1])
W <- names(baseline)

dat <- cbind(baseline, dat[, c(A, unlist(L), Y)])

sem <- Npsem$new(W, L, A, Y)
d <- odtr(dat, sem, 1, sl, sl, "binomial", FALSE)

saveRDS(d, "data/drv/optimal-rule-042823.rds")
# d <- readRDS("data/drv/optimal-rule-042823.rds")

shifted <- dat
for (a in A) {
    shifted[[a]] <- d[, a]
}

estims <- lmtp_sdr(
    dat,
    A, Y, W, L,
    shifted = shifted,
    outcome_type = "survival",
    folds = 1,
    learners_outcome = sl,
    learners_trt = sl,
    k = 1
)

saveRDS(estims, "data/drv/survival-combined-odtr-042823.rds")
