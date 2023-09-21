demog <- c("sex", "age", "xrace")

treatment_info <- c(
    "project",
    "site",
    "trt",
    "rand_dt",
    "end_of_detox",
    "medicine",
    "switched_meds",
    "med_switch_date",
    "never_initiated"
)

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

outcomes <- c(
    "relapse_overall",
    "relapse_date",
    "relapse_this_week",
    "relapse_12wk_date",
    "relapse_24wk_date",
    "opioiduse12",
    "fusedt12",
    "opioiduse24",
    "fusedt24"
)

visit_data <- c(
    "use_today",
    "when",
    "day_of_intervention",
    "medicine",
    "selfopioid",
    "uopioid",
    "bup_dose",
    "naloxone_dose",
    "met_dose",
    "naltrexone_dose",
    "any_dose"
)

weekly_indicators <- c(
    "use_this_week",
    "relapsed",
    "dose_change_since_last_week",
    "week_of_intervention",
    "max_dose_this_week", 
    "any_dose_this_week", 
    "min_nonzero_dose_this_week", 
    "dose_change_during_week"
)
