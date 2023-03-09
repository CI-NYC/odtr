library(haven)
library(tidyverse)
library(lubridate)

sows_raw <- read_sas("data/src/sow.sas7bdat")
xbot <- readRDS("data/drv/xbot.rds")

ids <- paste0("0", xbot$baseline$PATID)

sows <- 
    mutate(sows_raw, across(SOANX:SOUSENOW, as.numeric)) |> 
    (function(x) {
        mutate(x, SOWS_total = rowSums(select(x, SOANX:SOUSENOW), na.rm = TRUE))
    })() |> 
    select(PATID, RANDDT, PROTSEG, VISNO, SOWASMDT, SOASMTM, SOWS_total) |> 
    mutate(SOASMTM = hm(SOASMTM))

redefined <- 
    group_by(sows, PATID, SOWASMDT) |> 
    filter(row_number(as.numeric(SOASMTM)) == 1) |> 
    ungroup() |> 
    mutate(weeks = interval(RANDDT, SOWASMDT) %/% weeks(1)) |> 
    group_by(PATID, weeks) |> 
    filter(row_number(weeks) == 1) |> 
    ungroup() |> 
    filter(weeks <= 24) |> 
    select(PATID, weeks, SOWS_total) |> 
    pivot_wider(names_from = weeks, values_from = SOWS_total, names_sort = TRUE) |>
    mutate(`0` = ifelse(is.na(`0`), `-1`, `0`)) |> 
    select(-`-1`) |> 
    right_join(tibble(PATID = ids))

og <- 
    filter(sows, VISNO %in% str_pad(0:24, width = 2, pad = "0")) |> 
    mutate(VISNO = as.numeric(VISNO)) |> 
    select(PATID, VISNO, SOWS_total) |> 
    pivot_wider(names_from = VISNO, values_from = SOWS_total, names_sort = TRUE) |> 
    right_join(tibble(PATID = ids))

out <- list(
    original = og, 
    redefined = redefined
)

saveRDS(out, "data/drv/sows.rds")
