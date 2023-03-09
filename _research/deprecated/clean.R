library(tidyverse)

bup_raw <- vroom::vroom("data/src/bupwithbaseline.csv")

# weeks bup was dispensed
good_weeks <- c(0, 1, 2, 3, 4, 6, 8, 10, 12, 14, 16, 20, 24)

bup <- bup_raw |> 
    select(
        -1, -TRT, -GENDER,
        -starts_with("VISNO"), 
        -starts_with("RANDDT."), 
        -starts_with("visitnumber."), 
        -starts_with("relapsenew.y"), 
        -starts_with("relapseold"), 
        -starts_with("relapse."), 
        -starts_with("enddatenew"), 
        -ends_with(as.character(setdiff(0:24, good_weeks))),
        -starts_with("opioiduse."), 
        -starts_with("urineopioiduse"),
        -starts_with("begintimeearly"), 
        -starts_with("date."), 
        -starts_with("EODRUGDT."), 
        -starts_with("EOSTOP51."), 
        -MOMEDANY, -MOMEDBUP, -MOMEDNAL, 
        -re, -ageonsether
    ) |> 
    relocate(gender, .before = AGE) |> 
    mutate(PATID = as.character(PATID), 
           across(c("AgeOnset_NIC", "AgeOnset_HER"), \(x) na_if(x, ".")))

# relapse -----------------------------------------------------------------

relapse <- select(bup, PATID, starts_with("relapsenew"))
names(relapse) <- gsub("new.x", "", names(relapse))

relapse <- select(relapse, PATID, all_of(paste0("relapse.", setdiff(good_weeks, 0:2))))
relapse[is.na(relapse)] <- 1

bup <- select(bup, -starts_with("relapsenew"))

relapse_visit <- rep(NA_real_, nrow(bup))
for (i in 1:nrow(relapse)) {
    x <- match(1, relapse[i, paste0("relapse.", setdiff(good_weeks, 0:2))])
    relapse_visit[i] <- setdiff(good_weeks, 0:2)[x]
}

relapse$week_of_relapse <- relapse_visit
relapse <- relocate(relapse, week_of_relapse, .after = PATID)

# dose --------------------------------------------------------------------

dose <- select(bup, PATID, starts_with("D51DOSE"))
names(dose) <- gsub("D51DOSE", "dose", names(dose))

bup <- select(bup, -starts_with("D51DOSE"))

first_dose <- rep(NA_real_, nrow(bup))
for (i in 1:nrow(dose)) {
    x <- min(which(!is.na(dose[i, paste0("dose.", good_weeks)])))
    x <- ifelse(is.infinite(x), NA_real_, x)
    first_dose[i] <- good_weeks[x]
}

dose$first_dose_visit <- first_dose
dose <- relocate(dose, first_dose_visit, .after = PATID)

# filling in future missing with LOCF
for (x in 2:length(good_weeks)) {
    current <- dose[[paste0("dose.", good_weeks[x])]]
    old <- dose[[paste0("dose.", good_weeks[x - 1])]]
    new <- ifelse(is.na(current) & !is.na(old), old, current)
    data.table::set(dose, j = paste0("dose.", good_weeks[x]), value = new)
}

# # can fill LOCB from dose.3
# for (x in 2:0) {
#     current <- dose[[paste0("dose.", x)]]
#     future <- dose[[paste0("dose.", x + 1)]]
#     new <- ifelse(is.na(current) & !is.na(future), future, current)
#     data.table::set(dose, j = paste0("dose.", x), value = new)
# }

# HAMD --------------------------------------------------------------------

hamd <- select(bup, PATID, all_of(paste0("HAMSCORE.", c(1, 2, 3, 4, 8, 12, 16, 20)))) # |> 
#    mutate(across(starts_with("HAMSCORE"), \(score) as.numeric(!is.na(score))))

names(hamd) <- gsub("HAMSCORE", "hamd", names(hamd))

bup <- select(bup, -starts_with("HAMSCORE"))

# SOWS --------------------------------------------------------------------

sows <- select(bup, PATID, all_of(paste0("sows.", c(1, 2, 3, 4, 8, 12, 16, 20)))) # |> 
#     mutate(across(starts_with("sows"), \(score) as.numeric(!is.na(score))))

bup <- select(bup, -starts_with("sows."))

# craving -----------------------------------------------------------------

craving <- select(bup, PATID, all_of(paste0("craving.", good_weeks[2:(length(good_weeks) - 1)])))

bup <- select(bup, -starts_with("craving."))

# opioid use --------------------------------------------------------------

use <- select(bup, PATID, all_of(paste0("anyopioiduse.", good_weeks)))

bup <- select(bup, -starts_with("anyopioiduse."))

# pain --------------------------------------------------------------------

pain <- select(bup, PATID, all_of(paste0("QLPAIN.", c(4, 8, 12, 16, 20)))) # |> 
#     mutate(across(starts_with("QLPAIN"), \(score) as.numeric(!is.na(score))))

bup <- select(bup, -starts_with("QLPAIN"))

names(pain) <- gsub("QLPAIN", "pain", names(pain))

# baseline ----------------------------------------------------------------

bup <- rename(bup, 
    site = SITE, 
    age = AGE, 
    race = race_new, 
    edu = c_edu, 
    employment = c_employment, 
    primary_drug_cost = PrimDrugCost, 
    heroin_user = Her_ASI_TLFB, 
    iv_user = IV_user, 
    methadone_user = Met_ASI_TLFB, 
    methadone_prescription = CurrMetPres_ASI, 
    smoker = CurrSmoker_ASI, 
    amphetamine_user = Amp_ASI_TLFB, 
    sedative_user = Sed_ASI_TLFB,
    cannabis_user = Can_ASI_TLFB, 
    cocaine_crack_user = Coc_ASI_TLFB, 
    dsm5_alcohol = DSM5_Alc, 
    dsm5_amphetamine = DSM5_Amp, 
    dsm5_sedatives = DSM5_Sed, 
    dsm5_cannabis = DSM5_Can, 
    dsm5_cocaine = DSM5_Coc, 
    age_nicotine = AgeOnset_NIC, 
    age_any_drug = AgeOnset_any, 
    age_heroin = AgeOnset_HER, 
    first_treatment = FirstTrt_fix, 
    past_successful_trt = SuccessTrt, 
    past_successful_trt_bup_met = SuccessTrt_BupMet, 
    past_successful_trt_nal = SuccessTrt_Nal, 
    baseline_hamd = HAMD, 
    baseline_asp_score = ASPSCORE_fix, 
    baseline_qlanxdep = c_QLANXDEP, 
    psych_disorders = PsychHx_fix, 
    chronic_pain_6m = ChronicPain6M, 
    baseline_pain = c_QLPAIN, 
    probation = probation_fix, 
    friends_alc_problems = friends_alc, 
    live_with_alcoholic = live_alcoholic, 
    live_with_druguser = live_druguser,
    past_withdrawal_discomfort = MHOPIWDL,
    baseline_sow = SOWSCORE,
    randomization_timing = RAND_TIMING,
    relapse_24wks = relapse
)

bup <- mutate(bup,
    site = as.character(site), 
    hispanic = as.numeric(hispanic == "Hispanic"), 
    across(c(starts_with("past_successful"), 
             "relapse_24wks", "probation", 
             "chronic_pain_6m"), \(x) as.numeric(x == "Yes"))
)

xbot <- list(baseline = bup, 
             dose = dose, 
             relapse = relapse, 
             hamd = hamd, 
             pain = pain, 
             craving = craving, 
             sows = sows, 
             use = use)

saveRDS(xbot, "data/drv/xbot.rds")
