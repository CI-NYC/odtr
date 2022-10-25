library(tidyverse)
library(lmtp)
library(glue)

devtools::load_all("pkg")

imputed <- readRDS("data/src/clean•patients•imputed•010422.rds")
combined <- readRDS("data/src/clean•weeks•with•relapse•wide•010422.rds")

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
    select(who, all_of(demog), all_of(comorbidities)) |> 
    mutate(sex = if_else(sex == "female", 1, 0)) |> 
    process_missing()

W <- names(dat)[-1]
A <- glue("wk{3:11}.dose_increase_this_week")
L <- lapply(3:11, \(x) c(glue("wk{x-1}.dose_this_week"), glue("wk{x}.use_this_week")))
Y <- glue("wk{4:12}.relapse_this_week")

sl <- c("SL.glm", "SL.xgboost", "SL.earth", "SL.mean")

dat <- left_join(bup, dat)

baseline <- as.data.frame(model.matrix(reformulate(W), dat)[, -1])
W <- names(baseline)

dat <- cbind(baseline, dat[, c(A, unlist(L), Y)])

# sem <- Npsem$new(W, L, A, Y)
# d <- odtr(dat, sem, 1, sl, sl, "binomial", TRUE)
# 
# saveRDS(d, "data/drv/optimal-rule.rds")
d <- readRDS("data/drv/optimal-rule-cf.rds")

shifted <- dat
for (a in A) {
    shifted[[a]] <- d[, a]
}

future::plan(future::multisession)
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

future::plan(future::sequential)
saveRDS(estims, "data/drv/survival-combined-odtr.rds")
