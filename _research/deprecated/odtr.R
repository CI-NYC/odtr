library(tidyverse)
library(lmtp)
library(glue)
library(sl3)

devtools::load_all("pkg")

imputed <- readRDS("data/src/clean•patients•imputed•010422.rds")

# testing on first imputation for now...
bup <- 
    readRDS("data/src/clean•weeks•with•relapse•wide•010422.rds") |> 
    filter(medicine == "bup") |> 
    left_join(filter(mice::complete(imputed, 1), medicine == "bup"))

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

W <- c(demog, comorbidities)
A <- glue("wk{3:11}.dose_increase_this_week")
L <- lapply(3:11, \(x) c(glue("wk{x-1}.dose_this_week"), glue("wk{x}.use_this_week")))
Y <- glue("wk{4:(11 + 1)}.relapse_this_week")


sl <- c("SL.glm", "SL.xgboost", "SL.earth", "SL.mean")

# sem <- Npsem$new(W, L, A, Y)
# d <- odtr(bup, sem, 1, sl, sl, "binomial")
# saveRDS(d, "data/test.rds")

d <- readRDS("data/test.rds")
bup_d <- data.table::copy(bup)

for (a in colnames(d$A_opt)) {
    data.table::set(bup_d, j = a, value = d$A_opt[, a])
}

# lapply(colnames(d$A_opt), \(x) table(bup_d[[x]], bup[[x]], dnn = list("optimal", "observed")))

opt <- lmtp_sdr(
    bup, A, Y, W, L, 
    shifted = bup_d, 
    outcome_type = "survival", 
    folds = 1, 
    learners_outcome = sl, learners_trt = sl
)
