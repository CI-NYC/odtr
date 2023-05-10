library(lmtp)
library(glue)
library(devtools)
library(future)
library(mlr3extralearners)
suppressPackageStartupMessages(library(tidyverse))

load_all("odtr")

source("_research/covariates.R")

visits_wide <- readRDS("data/src/clean_weeks_with_relapse_wide_080922.rds")
imputed <- readRDS("data/src/clean_patients_imputed_080922.rds")

tau <- 5
A <- glue("wk{2:tau}.dose_increase_this_week")
L <- lapply(2:tau, function(x) c(glue("wk{x-1}.dose_this_week"), glue("wk{x}.use_this_week")))
Y <- glue("wk{3:(tau + 1)}.relapse_this_week")
learners <- c("glm", "glmnet", "earth", "lightgbm")
sl <- c("SL.glm", "SL.mean", "SL.lightgbm", "SL.glmnet", "SL.earth")
med <- "bup"

task_list <- expand.grid(imp = 1:5, 
                         med = med, 
                         stringsAsFactors = FALSE)

for (node in 1:5) {
    W <- c(demog, comorbidities)
    task <- task_list[node, ]
    
    baseline <- 
        mice::complete(imputed, task$imp)[, W] |> 
        model.matrix(reformulate(W), data = _) |> 
        as.data.frame()
    baseline <- cbind(mice::complete(imputed, task$imp)[, "who", drop = F], baseline[, -1])
    W <- names(baseline[, -1])
    sem <- Npsem$new(W, L, A, Y)
    
    observed <- left_join(visits_wide, baseline)
    observed <- observed[as.character(observed$medicine) == task$med, ]
    
    plan(multisession, workers = 10)
    d <- odtr(as.data.frame(observed), sem, 1, learners, learners, "binomial")
    plan(sequential)
    
    saveRDS(d, glue("data/drv/optimal-rule-{med}-week{tau+1}_{node}.rds"))
}

for (node in 1:5) {
    W <- c(demog, comorbidities)
    task <- task_list[node, ]
    
    baseline <- 
        mice::complete(imputed, task$imp)[, W] |> 
        model.matrix(reformulate(W), data = _) |> 
        as.data.frame()
    baseline <- cbind(mice::complete(imputed, task$imp)[, "who", drop = F], baseline[, -1])
    W <- names(baseline[, -1])
    sem <- Npsem$new(W, L, A, Y)
    
    observed <- left_join(visits_wide, baseline)
    observed <- observed[as.character(observed$medicine) == task$med, ]
    
    d <- readRDS(glue("data/drv/optimal-rule-{med}-week{tau+1}_{node}.rds"))
    
    shifted <- as.data.frame(observed)
    for (a in colnames(d$A_opt)) {
        shifted[, a] <- d$A_opt[, ..a]
    }
    
    sdr <- lmtp_sdr(
        as.data.frame(observed), 
        A, Y, W, L,
        shifted = shifted, 
        outcome_type = "survival", 
        folds = 1, 
        learners_outcome = sl, 
        learners_trt = sl
    )
    
    saveRDS(sdr, glue("data/drv/optimal-rule2-{med}-week{tau+1}-lmtp-{node}.rds"))
}
