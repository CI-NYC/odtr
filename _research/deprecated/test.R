library(sl3)
library(lmtp)
library(glue)
library(tidyverse)

devtools::load_all("pkg")

xbot <- readRDS("data/drv/xbot.rds")$baseline
craving <- readRDS("data/drv/craving.rds")$original
hamd <- readRDS("data/drv/hamd.rds")$original
qlpain <- readRDS("data/drv/qlpain.rds")$original
sows <- readRDS("data/drv/sows.rds")$original

xbot <- mutate(xbot, PATID = paste0("0", PATID))

weeks <- as.character(c(0, 3, 4, 8))

craving <- select(craving, PATID, any_of(weeks))
names(craving) <- c("PATID", paste0("craving_", names(craving)[-1]))
hamd <- select(hamd, PATID, any_of(weeks))
names(hamd) <- c("PATID", paste0("hamd_", names(hamd)[-1]))
qlpain <- select(qlpain, PATID, any_of(weeks))
names(qlpain) <- c("PATID", paste0("qlpain_", names(qlpain)[-1]))
sows <- select(sows, PATID, any_of(weeks))
names(sows) <- c("PATID", paste0("sows_", names(sows)[-1]))

process_missing <- function(data) {
    as.data.frame(cbind(PATID = data$PATID, sl3::make_sl3_Task(data, names(data)[-1])$data))
}

craving <- process_missing(craving)
hamd <- process_missing(hamd)
qlpain <- process_missing(qlpain)
sows <- process_missing(sows)

relapses <- 
    readRDS("data/drv/clean•weeks•with•relapse•wide•010422.rds") |> 
    filter(medicine == "bup", project == "51") |> 
    rename(PATID = who)

tau <- 8
A <- glue("wk{3:tau}.dose_increase_this_week")
# W <- c(demog, comorbidities, "site")


L <- list(
    c("wk2.dose_this_week", "wk3.use_this_week", "craving_3", "delta_craving_3"), 
    c("wk3.dose_this_week", "wk4.use_this_week", "craving_4", "delta_craving_4"), 
    c("wk4.dose_this_week", "wk5.use_this_week"), 
    c("wk5.dose_this_week", "wk6.use_this_week"), 
    c("wk6.dose_this_week", "wk7.use_this_week"), 
    c("wk7.dose_this_week", "wk8.use_this_week", "craving_8", "delta_craving_8")
)

Y <- glue("wk{c(4:tau, 12)}.relapse_this_week")

W <- c("age", "hispanic", "gender")

# W <- c("age", "hispanic", "gender", "heroin_user", "iv_user", "smoker", "amphetamine_user", 
#        "sedative_user", "cannabis_user", "cocaine_crack_user", 
#        paste0("dsm5_", c("alcohol", "amphetamine", "sedatives", "cannabis", "cocaine")), 
#        "baseline_hamd", "baseline_asp_score", "baseline_qlanxdep", "baseline_pain", 
#        "homeless", "past_withdrawal_discomfort")

tmp <- list(xbot, relapses, craving, hamd, qlpain, sows) |> 
    reduce(left_join)

# baseline <- as.data.frame(model.matrix(reformulate(W, intercept = FALSE), xbot$baseline[, W]))

# # one hot encoding
# tmp <- cbind(select(tmp, -all_of(W)), baseline)
# W <- names(baseline)

sem <- Npsem$new(W, A = A, Y = Y, L = L)

sl <- c("SL.glm", "SL.xgboost", "SL.earth")
optimal <- odtr(tmp, sem, 1, sl, sl, "binomial")

shifted <- tmp
for (a in A) {
    shifted[[a]] <- optimal$A_opt[, a]
}

sdr <- lmtp_sdr(
    tmp, 
    A, Y, W, L,
    shifted = shifted, 
    outcome_type = "survival", 
    folds = 1, 
    learners_outcome = sl, 
    learners_trt = sl, 
    k = 1
)

tmle <- lmtp_tmle(
    tmp, 
    A, Y, W, L,
    shifted = shifted, 
    outcome_type = "survival", 
    folds = 1, 
    learners_outcome = sl, 
    learners_trt = sl, 
    k = 1
)
