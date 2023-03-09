library(haven)
library(tidyverse)
library(lubridate)

crave_raw <- read_csv("data/src/ctn51vas1.csv")
sows_raw <- read_sas("data/src/sow.sas7bdat")
xbot <- readRDS("data/drv/xbot.rds")

ids <- paste0("0", xbot$baseline$PATID)

randdt <- 
    select(sows_raw, PATID, RANDDT) |> 
    distinct() |> 
    right_join(tibble(PATID = ids))

crave <- 
    mutate(crave_raw, VASASMDT = mdy(VASASMDT)) |> 
    rename(PATID = who)

og <- 
    filter(crave, VISNO %in% str_pad(0:24, width = 2, pad = "0")) |> 
    mutate(VISNO = as.numeric(VISNO)) |> 
    select(PATID, VISNO, VACRVOPI) |> 
    pivot_wider(names_from = VISNO, values_from = VACRVOPI, names_sort = TRUE) |> 
    right_join(tibble(PATID = ids))

redefined <- 
    group_by(crave, PATID, VASASMDT) |> 
    filter(row_number(VASASMDT) == 1) |> 
    ungroup() |> 
    right_join(randdt) |> 
    mutate(weeks = interval(RANDDT, VASASMDT) %/% weeks(1)) |> 
    group_by(PATID, weeks) |> 
    filter(row_number(weeks) == 1) |> 
    ungroup() |> 
    filter(weeks <= 24) |> 
    select(PATID, weeks, VACRVOPI) |> 
    pivot_wider(names_from = weeks, values_from = VACRVOPI, names_sort = TRUE) |> 
    mutate(`0` = ifelse(is.na(`0`), `-1`, `0`)) |> 
    select(-`-1`)

out <- list(
    original = og, 
    redefined = redefined
)

saveRDS(out, "data/drv/craving.rds")

# case_when(
#     VACRVOPI == 0 ~ 0, 
#     VACRVOPI > 0 & VACRVOPI < 6 ~ 1, 
#     VACRVOPI > 5 & VACRVOPI < 26 ~ 2, 
#     VACRVOPI > 25 ~ 3
# )
