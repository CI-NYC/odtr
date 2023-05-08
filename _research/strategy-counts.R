library(tidyverse)
library(lmtp)
library(glue)
library(kableExtra)

combined <- readRDS("data/drv/clean_weeks_with_relapse_wide_080922.rds")

d <- readRDS("data/drv/optimal-rule-050523.rds")
odtr <- readRDS("data/drv/survival-combined-odtr-042823.rds")
dynamic <- readRDS("data/drv/survival-combined-dynamic-030923.rds")
constant <- readRDS("data/drv/survival-combined-constant-030923.rds")

lmtp_contrast(constant, ref = odtr, type = "rr")
lmtp_contrast(constant, ref = dynamic, type = "rr")
lmtp_contrast(dynamic, ref = odtr, type = "rr")

bup <- filter(combined, medicine == "bup")

A <- glue("wk{2:11}.dose_increase_this_week")
Y <- glue("wk{3:12}.relapse_this_week")

# Was there previous time opioid use? 
condA <- bup[, glue("wk{2:11}.use_this_week")] == 1
colnames(condA) <- A

# Was dose under the allowable maximum? 
condB <- bup[, glue("wk{2:11}.dose_this_week")] < 32
colnames(condB) <- A

# Was dose under the dose threshold? 
condC <- bup[, glue("wk{2:11}.dose_this_week")] < 16
colnames(condC) <- A

d_dynamic <- bup[, A]
d_dynamic[, A] <- apply(condC | (condA & condB), 2, \(x) as.numeric(x), simplify = FALSE)

# Counts for following a strategy
strategy_n <- function(x) {
    colnames(x) <- c("Wk. 2", 3:11)
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
        `\\emph{Total}` = constant$density_ratios != 0, 
        `\\hspace{1em}$\\d_0^a$` = bup[, Y] == 0 & 
            d[, A] == 1 & 
            constant$density_ratios != 0, 
        `\\hspace{1em}$\\d1^a$` =  bup[, Y] == 0 & 
            d_dynamic[, A] == 1 & 
            constant$density_ratios != 0,
        `Increase`= odtr$density_ratios != 0 & increased,
        `Constant`= odtr$density_ratios != 0 & !increased,
        `\\emph{Total}` = odtr$density_ratios != 0, 
        `Increase`= dynamic$density_ratios != 0 & increased,
        `Constant`= dynamic$density_ratios != 0 & !increased,
        `\\emph{Total}` = dynamic$density_ratios != 0
    ), strategy_n, .id = "subset"
) |> 
    kbl(format = "latex", booktabs = TRUE, escape = FALSE) |>
    pack_rows("Constant dosing strategy", 1, 3) |> 
    pack_rows("Learning optimal strategy, $d_0$", 4, 6) |>
    pack_rows("Increase dose in response to use strategy, $d1$", 7, 9)

y_odtr <- 
    tidy(odtr) |> 
    mutate(estimate = 1 - estimate, 
           conf.low = estimate + qnorm(0.025)*std.error, 
           conf.high = estimate - qnorm(0.025)*std.error)

y_dyn <- 
    tidy(dynamic) |> 
    mutate(estimate = 1 - estimate, 
           conf.low = estimate + qnorm(0.025)*std.error, 
           conf.high = estimate - qnorm(0.025)*std.error)

y_cons <- 
    tidy(constant) |> 
    mutate(estimate = 1 - estimate, 
           conf.low = estimate + qnorm(0.025)*std.error, 
           conf.high = estimate - qnorm(0.025)*std.error)

y_d <- bind_rows(y_odtr, y_cons, y_dyn)

data.frame(Variable = c("constant", "dynamic", "odtr"), Visualization = "") |> 
    kbl(booktabs = TRUE, format = "latex") |> 
    column_spec(
        2,
        image = spec_pointrange(
            x = y_d$estimate,
            xmin = y_d$conf.low,
            xmax = y_d$conf.high,
            #vline = 1, 
            dir = "figures", 
            file_type = "pdf", 
            ann = TRUE
        )
    )

# \begin{tabular}[t]{l>{}l}
# \toprule
# Variable & Visualization\\
# \midrule
# constant & \includegraphics[width=0.67in, height=0.17in]{file:////Users/ntw2117/Documents/sdaepi/odtr/figures/pointrange_1149153bc8b24.pdf}\\
# dynamic & \includegraphics[width=0.67in, height=0.17in]{file:////Users/ntw2117/Documents/sdaepi/odtr/figures/pointrange_114912720dc10.pdf}\\
# odtr & \includegraphics[width=0.67in, height=0.17in]{file:////Users/ntw2117/Documents/sdaepi/odtr/figures/pointrange_1149142506386.pdf}\\
# \bottomrule
# \end{tabular}
