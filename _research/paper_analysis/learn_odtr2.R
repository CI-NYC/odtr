library(lmtp)
library(glue)
library(devtools)
library(doFuture)
library(tidyverse)
library(here)
library(glue)
library(odtr)

analysis_root <- "_research/paper_analysis"

# load covariates
source(here(analysis_root, "covariates.R"))

visits <- readRDS(here(analysis_root, "data/src/clean_weeks_with_relapse_wide_080922.rds"))
imputed <- readRDS(here(analysis_root, "data/src/clean_patients_imputed_080922.rds"))

tau <- 5
A <- glue("wk{2:tau}.dose_increase_this_week")
L <- lapply(2:tau, function(x) c(glue("wk{x-1}.dose_this_week"), glue("wk{x}.use_this_week")))
Y <- glue("wk{3:(tau + 1)}.relapse_this_week")

learners1 <- list("glm", "earth", list("nnet", trace = FALSE), "cv_glmnet")
# [glm, rpart, cv_glmnet]
learners2 <- c("cv_glmnet") 
med <- "bup"
folds <- 10

task_list <- expand.grid(imp = 1:5, 
                         med = med, 
                         stringsAsFactors = FALSE)

plan(multisession, workers = 5)

res <- foreach(node = 1:5, 
               .options.future = list(seed = TRUE, 
                                      packages = "mlr3extralearners")) %dofuture% {
    W <- c(demog, comorbidities)
    task <- task_list[node, ]
    
    baseline <- mice::complete(imputed, task$imp)[, W] |> 
        model.matrix(reformulate(W), data = _) |> 
        as.data.frame()
    
    baseline <- cbind(mice::complete(imputed, task$imp)[, "who", drop = F], baseline[, -1])
    W <- names(baseline[, -1])
    
    observed <- left_join(visits, baseline)
    observed <- observed[as.character(observed$medicine) == task$med, ] |> 
        as.data.frame()

    d <- odtr(observed, A, Y, W, L, learners1, learners1, learners2, folds, "binomial")
    
    shifted <- as.data.frame(observed)
    shifted[, A] <- d$A_opt[, A]
    
    sdr <- lmtp_sdr(
        as.data.frame(observed), 
        A, Y, W, L,
        shifted = shifted, 
        outcome_type = "survival", 
        folds = folds, 
        learners_outcome = learners1, 
        learners_trt = learners1
    )
    
    list(odtr = d, 
         theta = sdr)
}

plan(sequential)

saveRDS(res, here(analysis_root, glue("data/drv/odtr_dSL_{med}_{folds}_{learners2}.rds")))
