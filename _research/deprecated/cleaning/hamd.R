library(haven)
library(tidyverse)
library(lubridate)

ham_raw <- read_sas("data/src/ham.sas7bdat")
xbot <- readRDS("data/drv/xbot.rds")

ids <- xbot$baseline$PATID

redefined <- 
    select(ham_raw, PATID, RANDDT, VISNO, HAMASMDT, HAMSCORE) |> 
    mutate(weeks = interval(RANDDT, HAMASMDT) %/% weeks(1)) |> 
    filter(weeks <= 24) |> 
    right_join(tibble(PATID = paste0("0", ids))) |> 
    group_by(PATID, weeks) |> 
    filter(row_number(HAMASMDT) == 1) |> 
    ungroup() |> 
    select(PATID, weeks, HAMSCORE) |> 
    pivot_wider(names_from = weeks, values_from = HAMSCORE, names_sort = TRUE) |> 
    mutate(`0` = ifelse(is.na(`0`), `-1`, `0`)) |> 
    select(-`-1`)

# using_past_3 <- 
#     select(tmp, -PATID) |> 
#     slide(~ as.numeric(unlist(.x))) |> 
#     map(\(x) slide_dbl(x, mean, na.rm = TRUE, .before = 3)) |> 
#     map(\(x) ifelse(is.nan(x), NA_real_, x)) |> 
#     map_dfr(function (x) {
#         as.data.frame(as.list(x), col.names = as.character(0:24), 
#                       check.names = FALSE)
#     }) |> 
#     as_tibble() |> 
#     mutate(PATID = tmp$PATID, .before = 1)
    
og <- 
    select(ham_raw, PATID, VISNO, HAMSCORE) |> 
    filter(VISNO %in% str_pad(0:24, width = 2, pad = "0")) |> 
    mutate(VISNO = as.numeric(VISNO)) |> 
    pivot_wider(names_from = VISNO, values_from = HAMSCORE, names_sort = TRUE) |> 
    right_join(tibble(PATID = paste0("0", ids)))

out <- list(
    original = og, 
    redefined = redefined
)

saveRDS(out, "data/drv/hamd.rds")
