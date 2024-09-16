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

# seed <- readRDS(here(analysis_root, "data/drv/seed.rds"))

visits <- readRDS(here(analysis_root, "data/src/clean_weeks_with_relapse_wide_080922.rds"))
imputed <- readRDS(here(analysis_root, "data/src/clean_patients_imputed_080922.rds"))

# housekeeping ------------------------------------------------------------

tau <- 5
A <- glue("wk{2:tau}.dose_increase_this_week")
L <- lapply(2:tau, function(x) c(glue("wk{x-1}.dose_this_week"), glue("wk{x}.use_this_week")))
Y <- glue("wk{3:(tau + 1)}.relapse_this_week")

learners1 <- list("earth", list("nnet", trace = FALSE), "cv_glmnet", "xgboost", "ranger", "glm")
learners1 <- list("earth", list("nnet", trace = FALSE), "xgboost", "ranger", "glm")
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

observed_white <- filter(observed, xrace2 == 0, xrace3 == 0, xrace4 == 0)
observed_black <- filter(observed, xrace2 == 1)

# estimate ODTR -----------------------------------------------------------

W <- setdiff(W, c("xrace2", "xrace3", "xrace4"))
d_white <- odtr(observed_white, A, Y, W, L, learners1, learners1, learners2, filters, 1, "binomial")
d_black <- odtr(observed_black, A, Y, W, L, learners1, learners1, learners2, filters, 1, "binomial")

shifted_white <- as.data.frame(observed_white)
shifted_white[, A] <- d_white$A_opt

theta_white <- lmtp_sdr(as.data.frame(observed_white),
                        A, Y, W, L,
                        shifted = shifted_white,
                        outcome_type = "survival",
                        folds = 1,
                        learners_outcome = learners1,
                        learners_trt = learners1)

shifted_black <- as.data.frame(observed_black)
shifted_black[, A] <- d_black$A_opt

theta_black <- lmtp_sdr(as.data.frame(observed_black),
                        A, Y, W, L,
                        shifted = shifted_black,
                        outcome_type = "survival",
                        folds = 1,
                        learners_outcome = learners1,
                        learners_trt = learners1)

obs_white <- 1 - mean(observed_white[[Y[length(Y)]]], na.rm = T)
obs_black <- 1 - mean(observed_black[[Y[length(Y)]]], na.rm = T)

data.frame(cat = c("Observed", "Observed", "ODTR", "ODTR"), 
           race = c("White", "Black", "White", "Black"),
           theta = c(obs_white, obs_black, theta_white$theta, theta_black$theta),
           conf.low = c(obs_white, obs_black, theta_white$low, theta_black$low), 
           conf.high = c(obs_white, obs_black, theta_white$high, theta_black$high)) |> 
    ggplot(aes(x = race, y = theta, shape = cat)) + 
    geom_point(position = position_dodge(width=0.5), 
               size = 4) + 
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                  width = 0, 
                  position = position_dodge(width=0.5)) + 
    labs(x = "Race", 
         y = "Wk. 6 Survival", 
         shape = NULL) + 
    theme_bw()

se_wb <- sqrt(theta_white$standard_error^2 + theta_black$standard_error^2)
ci_wb <- (theta_white$theta - theta_black$theta) + c(-1, 1)*qnorm(0.975)*se_wb

data.frame(cat = c("Observed", "Observed", "Observed", "ODTR", "ODTR", "ODTR"), 
           race = c("White", "Black", "White - Black", "White", "Black", "White - Black"),
           theta = c(obs_white, obs_black, 
                     obs_white - obs_black,
                     theta_white$theta, theta_black$theta, 
                     theta_white$theta - theta_black$theta),
           conf.low = c(obs_white, obs_black, obs_white - obs_black,
                        theta_white$low, theta_black$low, ci_wb[1]), 
           conf.high = c(obs_white, obs_black, obs_white - obs_black, 
                         theta_white$high, theta_black$high, ci_wb[2])) |> 
    ggplot(aes(x = cat, y = theta, shape = race)) + 
    geom_point(position = position_dodge(width=0.5), 
               size = 4) + 
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                  width = 0, 
                  position = position_dodge(width=0.5)) + 
    labs(x = NULL, 
         y = "Wk. 6 Survival", 
         shape = "Race") + 
    theme_bw()

data.frame(cat = c("Observed", "ODTR"), 
           race = c("White - Black", "White - Black"),
           theta = c(obs_white - obs_black,
                     theta_white$theta - theta_black$theta),
           conf.low = c(obs_white - obs_black,
                        ci_wb[1]), 
           conf.high = c(obs_white - obs_black, ci_wb[2])) |> 
    ggplot(aes(x = cat, y = theta)) + 
    geom_point(size = 4) + 
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                  width = 0) + 
    labs(x = NULL, 
         y = "Difference in Wk. 6 Survival, \nWhite vs. Black beneficiaries") + 
    theme_bw()
