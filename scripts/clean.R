library(tidyverse)

bup_raw <- vroom::vroom("data/src/bupwithbaseline.csv")

# D51DOSE.{time} is the dose at that visit? 
# RANDDT.{time} appears to be constant
# VISNObup.{time} appears to be constant
# Using relapsenew.x.{time} for relapse defintions
# HAMSCORE by visits likely has too much missingness

bup <- bup_raw |> 
    select(-1, -starts_with("VISNO"), -starts_with("RANDDT."), 
           -starts_with("sows."), -starts_with("visitnumber."), 
           -starts_with("relapsenew.y"), -starts_with("relapseold"), 
           -starts_with("relapse."), -starts_with("enddatenew"), 
           -TRT) |> 
    mutate(PATID = as.character(PATID), 
           across(c("AgeOnset_NIC", "AgeOnset_HER"), \(x) na_if(x, ".")))

# relapse -----------------------------------------------------------------

relapse <- select(bup, PATID, starts_with("relapsenew"))
names(relapse) <- gsub("new.x", "", names(relapse))

relapse <- select(relapse, PATID, all_of(paste0("relapse.", 3:24)))
relapse[is.na(relapse)] <- 1

bup <- select(bup, -starts_with("relapsenew"))

# dose --------------------------------------------------------------------

dose <- select(bup, PATID, starts_with("D51DOSE"))
names(dose) <- gsub("D51DOSE", "dose", names(dose))

bup <- select(bup, -starts_with("D51DOSE"))

for (x in 1:24) {
    current <- dose[[paste0("dose.", x)]]
    old <- dose[[paste0("dose.", x - 1)]]
    new <- ifelse(is.na(current) & !is.na(old), old, current)
    data.table::set(dose, j = paste0("dose.", x), value = new)
}
