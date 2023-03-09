library(lmtp)
library(glue)
library(tidyverse)

devtools::load_all("pkg")

process_missing <- function(data) {
    as.data.frame(cbind(PATID = data$PATID, sl3::make_sl3_Task(data, names(data)[-1])$data))
}

imputed <- readRDS("data/src/clean•patients•imputed•010422.rds")

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

dat <- imputed$data
dat <- 
    filter(dat, medicine == "bup", project == "51") |> 
    select(PATID = who, all_of(demog), all_of(comorbidities)) |> 
    mutate(sex = if_else(sex == "female", 1, 0)) |> 
    process_missing()

craving <- readRDS("data/drv/craving.rds")$original
hamd <- readRDS("data/drv/hamd.rds")$original
qlpain <- readRDS("data/drv/qlpain.rds")$original
sows <- readRDS("data/drv/sows.rds")$original

craving <- select(craving, PATID, any_of(c("0", "2", "3", "4", "5", "6", "7", "8")))
names(craving) <- c("PATID", paste0("craving_", names(craving)[-1]))

hamd <- select(hamd, PATID, any_of(c("0", "2", "3", "4", "8")))
names(hamd) <- c("PATID", paste0("hamd_", names(hamd)[-1]))

qlpain <- select(qlpain, PATID, any_of(c("0", "4", "8")))
names(qlpain) <- c("PATID", paste0("qlpain_", names(qlpain)[-1]))

sows <- select(sows, PATID, any_of(c("0", "2", "3", "4", "8")))
names(sows) <- c("PATID", paste0("sows_", names(sows)[-1]))

craving <- process_missing(craving)
hamd <- process_missing(hamd)
qlpain <- process_missing(qlpain)
sows <- process_missing(sows)

qlpain <- 
    model.matrix(reformulate(names(qlpain)[-1]), data = qlpain) |> 
    as.data.frame() |> 
    (\(x) x[, -1])() |> 
    (\(x) cbind(PATID = qlpain$PATID, x))()

relapses <- 
    readRDS("data/drv/clean•weeks•with•relapse•wide•010422.rds") |> 
    filter(medicine == "bup", project == "51") |> 
    rename(PATID = who)

tmp <- 
    list(dat, craving, hamd, qlpain, sows, relapses) |> 
    reduce(left_join) |> 
    mutate(across(starts_with("hamd"), \(x) as.numeric(as.character(x))))

A <- glue("wk{3:8}.dose_increase_this_week")

W <- c(names(dat)[-c(1, 5)], "craving_0", "hamd_0", "qlpain_02", "qlpain_03", "sows_0")

L <- list(
    c("wk2.dose_this_week", "wk3.use_this_week", "craving_3", "delta_craving_3", "hamd_3", "delta_hamd_3", "sows_3", "delta_sows_3"), 
    c("wk3.dose_this_week", "wk4.use_this_week", "craving_4", "delta_craving_4", "hamd_4", "delta_hamd_4", "qlpain_42", "qlpain_43", "delta_qlpain_4", "sows_4", "delta_sows_4"), 
    c("wk4.dose_this_week", "wk5.use_this_week", "craving_5", "delta_craving_5"), 
    c("wk5.dose_this_week", "wk6.use_this_week", "craving_6", "delta_craving_6"), 
    c("wk6.dose_this_week", "wk7.use_this_week", "craving_7", "delta_craving_7"), 
    c("wk7.dose_this_week", "wk8.use_this_week", "craving_8", "delta_craving_8", "hamd_8", "delta_hamd_8", "qlpain_82", "qlpain_83", "delta_qlpain_8", "sows_8", "delta_sows_8")
)

Y <- glue("wk{c(4:8, 12)}.relapse_this_week")

sem <- Npsem$new(W, A = A, Y = Y, L = L)

sl <- c("SL.glm", "SL.xgboost", "SL.earth", "SL.mean")
d <- odtr(tmp, sem, 1, sl, sl, "binomial")

saveRDS(d, "data/drv/optimal-rules-xbot.rds")

shifted <- tmp
for (a in A) {
    shifted[[a]] <- d[, a]
}

sdr <- lmtp_sdr(
    tmp, 
    A, Y, W, L, 
    shifted = shifted, 
    learners_outcome = sl, 
    learners_trt = sl,
    outcome_type = "survival",
    k = 1, folds = 1
)

