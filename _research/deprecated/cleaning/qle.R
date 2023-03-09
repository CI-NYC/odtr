library(haven)
library(tidyverse)
library(lubridate)

qle_raw <- read_sas("data/src/qle.sas7bdat")
xbot <- readRDS("data/drv/xbot.rds")

ids <- paste0("0", xbot$baseline$PATID)

qle <- 
    select(qle_raw, PATID, RANDDT, QLEASMDT, VISNO, QLPAIN) |> 
    right_join(tibble(PATID = ids)) |> 
    mutate(QLPAIN = na_if(QLPAIN, ""))

og <- 
    filter(qle, VISNO %in% str_pad(0:24, width = 2, pad = "0")) |> 
    mutate(VISNO = as.numeric(VISNO)) |> 
    select(PATID, VISNO, QLPAIN) |> 
    pivot_wider(names_from = VISNO, values_from = QLPAIN, names_sort = TRUE) |> 
    right_join(tibble(PATID = ids))

redefined <- 
    group_by(qle, PATID, QLEASMDT) |> 
    filter(row_number(as.numeric(QLEASMDT)) == 1) |> 
    ungroup() |> 
    mutate(weeks = interval(RANDDT, QLEASMDT) %/% weeks(1)) |> 
    group_by(PATID, weeks) |> 
    filter(row_number(weeks) == 1) |> 
    ungroup() |> 
    filter(weeks <= 24) |> 
    select(PATID, weeks, QLPAIN) |> 
    pivot_wider(names_from = weeks, values_from = QLPAIN, names_sort = TRUE) |>
    mutate(`0` = ifelse(is.na(`0`), `-1`, `0`)) |> 
    select(-`-1`) |> 
    right_join(tibble(PATID = ids))

out <- list(
    original = og, 
    redefined = redefined
)

saveRDS(out, "data/drv/qlpain.rds")
