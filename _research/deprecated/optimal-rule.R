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

weeks <- as.character(c(0, 3, 4, 8, 12, 16, 20))

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

dose_weeks <- weeks
dose_weeks[1] <- 2
for (x in 2:length(dose_weeks)) {
    relapses[[glue("wk{dose_weeks[x]}.increase")]] <- 
        as.numeric(
            relapses[[glue("wk{dose_weeks[x]}.dose_this_week")]] > 
                relapses[[glue("wk{dose_weeks[x - 1]}.dose_this_week")]]
        )
}

relapses <- select(relapses, -ends_with("dose_increase_this_week"))

A <- glue("wk{weeks[2:(length(weeks))]}.increase")
Y <- glue("wk{c(4, 8, 12, 16, 20, 24)}.relapse_this_week")

W <- c("age", "hispanic", "gender")

# W <- c("age", "hispanic", "gender", "heroin_user", "iv_user", "smoker", "amphetamine_user", 
#        "sedative_user", "cannabis_user", "cocaine_crack_user", 
#        paste0("dsm5_", c("alcohol", "amphetamine", "sedatives", "cannabis", "cocaine")), 
#        "baseline_hamd", "baseline_asp_score", "baseline_qlanxdep", "baseline_pain", 
#        "homeless", "past_withdrawal_discomfort")

L <- list(
    c("wk3.use_this_week", "wk3.dose_this_week", "craving_3", "delta_craving_3"), 
    c("wk4.use_this_week", "wk4.dose_this_week", "craving_4", "delta_craving_4"), 
    c("wk8.use_this_week", "wk8.dose_this_week", "craving_8", "delta_craving_8"), 
    c("wk12.use_this_week", "wk12.dose_this_week", "craving_12", "delta_craving_12"), 
    c("wk16.use_this_week", "wk16.dose_this_week", "craving_16", "delta_craving_16"), 
    c("wk20.use_this_week", "wk20.dose_this_week", "craving_20", "delta_craving_20") 
)

# L <- list(
#     c("wk3.use_this_week", "wk3.dose_this_week", "craving_3", "hamd_3", "sows_3"), 
#     c("wk4.use_this_week", "wk4.dose_this_week", "craving_4", "hamd_4", "qlpain_4", "sows_4"), 
#     c("wk8.use_this_week", "wk8.dose_this_week", "craving_8", "hamd_8", "qlpain_8", "sows_8"), 
#     c("wk12.use_this_week", "wk12.dose_this_week", "craving_12", "hamd_12", "qlpain_12", "sows_12"), 
#     c("wk16.use_this_week", "wk16.dose_this_week", "craving_16", "hamd_16", "qlpain_16", "sows_16"), 
#     c("wk20.use_this_week", "wk20.dose_this_week", "craving_20", "hamd_20", "qlpain_20", "sows_20") 
# )

tmp <- list(xbot, relapses, craving, hamd, qlpain, sows) |> reduce(left_join)

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

gcomp <- lmtp_sub(
    tmp, 
    A, Y, W, L,
    shifted = shifted, 
    outcome_type = "survival", 
    folds = 1, 
    learners = sl, 
    k = 1
)

res <- lmtp_sdr(
    tmp, 
    A, Y, W, L,
    shifted = shifted, 
    outcome_type = "survival", 
    folds = 1, 
    learners_outcome = sl, 
    learners_trt = sl, 
    k = 1
)

# curve -------------------------------------------------------------------

opt <- map(3:10, function(x) {
    A <- glue("wk{weeks[2:x]}.increase")
    Y <- glue("wk{weeks[3:(x + 1)]}.relapse_this_week")
    L <- lapply(weeks[2:x], 
                \(x) c(glue("wk{x}.use_this_week"), glue("wk{x}.dose_this_week")))
    lmtp_sdr(
        tmp, 
        A, Y, W, L,
        shifted = shifted,
        outcome_type = "survival", 
        folds = 1, 
        learners_outcome = Sl, 
        learners_trt = Sl
    )
})

on <- map(3:10, function(x) {
    A <- glue("wk{weeks[2:x]}.increase")
    Y <- glue("wk{weeks[3:(x + 1)]}.relapse_this_week")
    L <- lapply(weeks[2:x], 
                \(x) c(glue("wk{x}.use_this_week"), glue("wk{x}.dose_this_week")))
    lmtp_sdr(
        tmp, 
        A, Y, W, L,
        shift = static_binary_on,
        outcome_type = "survival", 
        folds = 1, 
        learners_outcome = Lrnr_xgboost$new(), 
        learners_trt = Lrnr_xgboost$new()
    )
})

off <- map(3:10, function(x) {
    A <- glue("wk{weeks[2:x]}.increase")
    Y <- glue("wk{weeks[3:(x + 1)]}.relapse_this_week")
    L <- lapply(weeks[2:x], 
                \(x) c(glue("wk{x}.use_this_week"), glue("wk{x}.dose_this_week")))
    lmtp_sdr(
        tmp, 
        A, Y, W, L,
        shift = static_binary_off,
        outcome_type = "survival", 
        folds = 1, 
        learners_outcome = Lrnr_xgboost$new(), 
        learners_trt = Lrnr_xgboost$new()
    )
})

obs <- map(3:10, function(x) {
    A <- glue("wk{weeks[2:x]}.increase")
    Y <- glue("wk{weeks[3:(x + 1)]}.relapse_this_week")
    L <- lapply(weeks[2:x], 
                \(x) c(glue("wk{x}.use_this_week"), glue("wk{x}.dose_this_week")))
    lmtp_sdr(
        tmp, 
        A, Y, W, L,
        shift = NULL,
        outcome_type = "survival", 
        folds = 1, 
        learners_outcome = Lrnr_xgboost$new(), 
        learners_trt = Lrnr_xgboost$new()
    )
})
