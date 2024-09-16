library(here)
library(lmtp)
library(glue)
library(tidyverse)
library(kableExtra)

analysis_root <- "_research/paper_analysis"

med <- "bup"
res <- readRDS(here(analysis_root, glue("data/drv/odtr_dSL_{med}_cf.rds")))

combined <- readRDS(here(analysis_root, "data/src/clean_weeks_with_relapse_wide_080922.rds"))
odtr <- res$theta_odtr
constant <- res$theta_constant
dynamic <- res$theta_rudolph

# strategy counts ---------------------------------------------------------

d_dynamic <- res$rudolph
d <- res$odtr_cv$A_opt

bup <- filter(combined, medicine == "bup")

tau <- 5
A <- glue("wk{2:tau}.dose_increase_this_week")
Y <- glue("wk{3:(tau + 1)}.relapse_this_week")

# Counts for following a strategy
strategy_n <- function(x) {
    colnames(x) <- c("Wk. 2", 3:tau)
    colSums(x, na.rm = TRUE) |>
        as_tibble(rownames = "time") |>
        pivot_wider(
            names_from = "time",
            values_from = "value"
        )
}

increased <- bup[, A] == 1

map_dfr(
    list(
        `Increase`= odtr$density_ratios != 0 & increased,
        `Constant`= odtr$density_ratios != 0 & !increased,
        `\\emph{Total}` = odtr$density_ratios != 0, 
        `\\emph{Total}` = constant$density_ratios != 0, 
        `\\hspace{1em}$\\d_0^a$` = bup[, Y] == 0 & 
            d[, A] == 1 & 
            constant$density_ratios != 0, 
        `\\hspace{1em}$\\d1^a$` =  bup[, Y] == 0 & 
            d_dynamic[, A] == 1 & 
            constant$density_ratios != 0,
        `Increase`= dynamic$density_ratios != 0 & increased,
        `Constant`= dynamic$density_ratios != 0 & !increased,
        `\\emph{Total}` = dynamic$density_ratios != 0
    ), strategy_n, .id = "subset"
) |> 
    kbl(format = "latex", booktabs = TRUE, escape = FALSE) |>
    pack_rows("Learned ODTR, $\\mathbf{D}_{n, \\text{opt}}$", 1, 3) |>
    pack_rows("Constant dosing rule", 4, 6) |> 
    pack_rows("\\citet{rudolph2022dose} dynamic dosing rule, $\\d_1$", 7, 9)

# results -----------------------------------------------------------------

risk_ratio <- function(x, y) {
    x$theta <- 1 - x$theta
    y$theta <- 1 - y$theta
    lmtp_contrast(x, ref = y, type = "rr")
}

risk_ratio(odtr, constant)
risk_ratio(odtr, dynamic)
risk_ratio(dynamic, constant)

results <- map_dfr(list(constant, dynamic, odtr), tidy) |> 
    mutate(estimate = 1 - estimate, 
           conf.low = estimate + qnorm(0.025)*std.error, 
           conf.high = estimate - qnorm(0.025)*std.error)

data.frame(Variable = c("constant", "dynamic", "odtr"), Visualization = "") |> 
    kbl(booktabs = TRUE, format = "latex") |> 
    column_spec(
        2,
        image = spec_pointrange(
            x = results$theta,
            xmin = results$conf.low,
            xmax = results$conf.high,
            #vline = 1, 
            dir = "figures", 
            file_type = "pdf", 
            ann = TRUE
        )
    )

odtr_lasso_coef <- function(fit_t) {
    coef(fit_t[[1]]$learners$regr.cv_glmnet$model$regr.cv_glmnet$model)
}

lapply(res$odtr$decision_fits, odtr_lasso_coef)

