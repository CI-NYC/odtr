library(devtools)
library(lmtp)
library(glue)
suppressPackageStartupMessages(library(tidyverse))
library(kableExtra)

source_gist("https://gist.github.com/nt-williams/295b2ba8d4dc3acfb588865dd1dbdc41")

med <- "bup"
node <- 1:5
tau <- 5

combined <- readRDS("data/drv/clean_weeks_with_relapse_wide_080922.rds")

odtr <- glue("data/drv/optimal-rule-{med}-week{tau+1}-lmtp-{node}.rds") |> 
    map(readRDS)

constant <- glue("data/src/{node}_{med}_{tau}_constant.rds") |> 
    map(readRDS) |> 
    map("fit")

dynamic <- glue("data/src/{node}_{med}_{tau}_hybrid.rds") |> 
    map(readRDS) |> 
    map("fit")

# strategy counts ---------------------------------------------------------
bup <- filter(combined, medicine == "bup")

A <- glue("wk{2:tau}.dose_increase_this_week")
Y <- glue("wk{3:(tau + 1)}.relapse_this_week")

# Was there previous time opioid use? 
condA <- bup[, glue("wk{2:tau}.use_this_week")] == 1
colnames(condA) <- A

# Was dose under the allowable maximum? 
condB <- bup[, glue("wk{2:tau}.dose_this_week")] < 32
colnames(condB) <- A

# Was dose under the dose threshold? 
condC <- bup[, glue("wk{2:tau}.dose_this_week")] < 16
colnames(condC) <- A

d_hybrid <- bup[, A]
d_hybrid[, A] <- apply(condC | (condA & condB), 2, \(x) as.numeric(x), simplify = FALSE)

d <- as.data.frame(readRDS(glue("data/drv/optimal-rule-{med}-week{tau+1}_1.rds"))$A_opt)

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
        `Increase`= odtr[[1]]$density_ratios != 0 & increased,
        `Constant`= odtr[[1]]$density_ratios != 0 & !increased,
        `\\emph{Total}` = odtr[[1]]$density_ratios != 0, 
        `\\emph{Total}` = constant[[1]]$density_ratios != 0, 
        `\\hspace{1em}$\\d_0^a$` = bup[, Y] == 0 & 
            d[, A] == 1 & 
            constant[[1]]$density_ratios != 0, 
        `\\hspace{1em}$\\d1^a$` =  bup[, Y] == 0 & 
            d_hybrid[, A] == 1 & 
            constant[[1]]$density_ratios != 0,
        `Increase`= dynamic[[1]]$density_ratios != 0 & increased,
        `Constant`= dynamic[[1]]$density_ratios != 0 & !increased,
        `\\emph{Total}` = dynamic[[1]]$density_ratios != 0
    ), strategy_n, .id = "subset"
) |> 
    kbl(format = "latex", booktabs = TRUE, escape = FALSE) |>
    pack_rows("Learned ODTR, $\\mathbf{D}_{n, \\text{opt}}$", 1, 3) |>
    pack_rows("Constant dosing rule", 4, 6) |> 
    pack_rows("\\citet{rudolph2022dose} dynamic dosing rule, $\\d_1$", 7, 9)

# results -----------------------------------------------------------------

map2(odtr, constant, function(x, y) {
    x$theta <- 1 - x$theta
    y$theta <- 1 - y$theta
    lmtp_contrast(x, ref = y, type = "rr")
}) |> 
    rubins_rules(label = "constant v. odtr", log = T)

map2(odtr, dynamic, function(x, y) {
    x$theta <- 1 - x$theta
    y$theta <- 1 - y$theta
    lmtp_contrast(x, ref = y, type = "rr")
}) |> 
    rubins_rules(label = "odtr v. dynamic", log = T)

map2(dynamic, constant, function(x, y) {
    x$theta <- 1 - x$theta
    y$theta <- 1 - y$theta
    lmtp_contrast(x, ref = y, type = "rr")
}) |> 
    rubins_rules(label = "constant v. dynamic", log = T)

results <- bind_rows(rubins_rules(constant, label = "constant"),
                     rubins_rules(dynamic, label = "hybrid"),
                     rubins_rules(odtr, label = "odtr")) |> 
    mutate(theta = 1 - theta, 
           conf.low = theta + qnorm(0.025)*se, 
           conf.high = theta - qnorm(0.025)*se)

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
