library(lmtp)
library(glue)
library(devtools)
library(doFuture)
library(tidyverse)
library(here)
library(glue)
# library(odtr)
library(mlr3pipelines)
library(mlr3filters)
library(mlr3extralearners)

load_all()

analysis_root <- "_research/paper_analysis"

# load data ---------------------------------------------------------------

# load covariates
source(here(analysis_root, "covariates.R"))

seed <- readRDS(here(analysis_root, "data/drv/seed.rds"))

visits <- readRDS(here(analysis_root, "data/src/clean_weeks_with_relapse_wide_080922.rds"))
imputed <- readRDS(here(analysis_root, "data/src/clean_patients_imputed_080922.rds"))

# housekeeping ------------------------------------------------------------

tau <- 5
A <- glue("wk{2:tau}.dose_increase_this_week")
L <- lapply(2:tau, function(x) c(glue("wk{x-1}.dose_this_week"), glue("wk{x}.use_this_week")))
Y <- glue("wk{3:(tau + 1)}.relapse_this_week")

learners1 <- list("earth", list("nnet", trace = FALSE), "cv_glmnet", "lightgbm")
learners2 <- "cv_glmnet"
filters <- po("filter", filter = flt("mim"), filter.nfeat = 10)

med <- "bup"

task_list <- expand.grid(imp = 1:5, med = med, stringsAsFactors = FALSE)

node <- 1

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

folds <- origami::make_folds(n = nrow(observed), V = 3, strata_ids = observed$wk6.relapse_this_week)
train <- observed[folds[[1]]$training_set, ]
holdout <- observed[folds[[1]]$validation_set, ]

# estimate ODTR -----------------------------------------------------------

d <- odtr(train, A, Y, W, L, learners1, learners1, learners2, filters, 1, "binomial")

# evaluate ODTR -----------------------------------------------------------

treat <- function(x) ifelse(x > 0, 1, 0)
rules <- matrix(data = NA_real_, ncol = length(A), nrow = nrow(holdout))

rules[holdout[["wk5.relapse_this_week"]] == 0, 4] <- 
    treat(predict(d$decision_fits[[4]][[1]], holdout[holdout[["wk5.relapse_this_week"]] == 0, ]))

rules[holdout[["wk4.relapse_this_week"]] == 0, 3] <- 
    treat(predict(d$decision_fits[[3]][[1]], holdout[holdout[["wk4.relapse_this_week"]] == 0, ]))

rules[holdout[["wk3.relapse_this_week"]] == 0, 2] <- 
    treat(predict(d$decision_fits[[2]][[1]], holdout[holdout[["wk3.relapse_this_week"]] == 0, ]))

rules[, 1] <- treat(predict(d$decision_fits[[1]][[1]], holdout))

shifted <- as.data.frame(holdout)
shifted[, A] <- rules

set.seed(seed)

theta_odtr <- lmtp_sdr(as.data.frame(holdout),
                        A, Y, W, L,
                        shifted = shifted,
                        outcome_type = "survival",
                        folds = 1,
                        learners_outcome = learners1,
                        learners_trt = learners1)

# evaluate Rudolph et al --------------------------------------------------

# Was there previous time opioid use? 
condA <- holdout[, glue("wk{2:tau}.use_this_week")] == 1

# Was dose under the allowable maximum? 
condB <- matrix(FALSE, dim(condA)[1], dim(condA)[2])

condB[holdout$medicine == "bup", ] <- 
    holdout[holdout$medicine == "bup", glue("wk{2:tau}.dose_this_week")] < 32

# Was dose under the dose threshold? 
condC <- matrix(FALSE, dim(condA)[1], dim(condA)[2])
condC[holdout$medicine == "bup", ] <- 
    holdout[holdout$medicine == "bup", glue("wk{2:tau}.dose_this_week")] < 16

A <- glue("wk{2:tau}.dose_increase_this_week")

colnames(condA) <- A
colnames(condB) <- A
colnames(condC) <- A

hybrid <- matrix(data = NA_real_, ncol = length(A), nrow = nrow(holdout))

hybrid[holdout[["wk5.relapse_this_week"]] == 0, 4] <- 
    apply(condC | (condB & condA), 2, function(x) as.numeric(x))[holdout[["wk5.relapse_this_week"]] == 0, 4]

hybrid[holdout[["wk4.relapse_this_week"]] == 0, 3] <- 
    apply(condC | (condB & condA), 2, function(x) as.numeric(x))[holdout[["wk4.relapse_this_week"]] == 0, 3]

hybrid[holdout[["wk3.relapse_this_week"]] == 0, 2] <- 
    apply(condC | (condB & condA), 2, function(x) as.numeric(x))[holdout[["wk3.relapse_this_week"]] == 0, 2]

hybrid[, 1] <- apply(condC | (condB & condA), 2, function(x) as.numeric(x))[, 1]

shifted2 <- as.data.frame(holdout)
shifted2[, A] <- hybrid

theta_rudolph <- lmtp_sdr(as.data.frame(holdout),
                           A, Y, W, L,
                           shifted = shifted2,
                           outcome_type = "survival",
                           folds = 1,
                           learners_outcome = learners1,
                           learners_trt = learners1)

# evaluate constant dose --------------------------------------------------

theta_constant <- lmtp_sdr(as.data.frame(holdout),
                            A, Y, W, L,
                            shift = static_binary_off,
                            outcome_type = "survival",
                            folds = 1,
                            learners_outcome = learners1,
                            learners_trt = learners1) 

# save results ------------------------------------------------------------

res <- list(odtr = d, 
            theta_odtr = theta_odtr, 
            rudolph = hybrid, 
            theta_rudolph = theta_rudolph, 
            theta_constant = theta_constant)

# saveRDS(res, here(analysis_root, glue("data/drv/odtr_dSL_{med}_holdout.rds")))
