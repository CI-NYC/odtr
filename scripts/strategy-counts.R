library(tidyverse)
library(lmtp)
library(glue)
library(kableExtra)

imputed <- readRDS("data/src/clean•patients•imputed•010422.rds")
combined <- readRDS("data/src/clean•weeks•with•relapse•wide•010422.rds")

dat <- imputed$data

bup <- filter(combined, medicine == "bup")

# obs <- tibble::tibble(
#     time = 4:12, 
#     strategy = "Observed", 
#     estimate = map_dbl(Y, \(x) mean(bup[[x]]))
# )

demog <- c("sex", "age", "xrace")

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

process_missing <- function(data) {
    as.data.frame(cbind(who = data$who, sl3::make_sl3_Task(data, names(data)[-1])$data))
}

dat <- 
    filter(dat, medicine == "bup") |> 
    select(who, all_of(demog), all_of(comorbidities)) |> 
    mutate(sex = if_else(sex == "female", 1, 0)) |> 
    process_missing()

W <- names(dat)[-1]
A <- glue("wk{3:11}.dose_increase_this_week")
L <- lapply(3:11, \(x) c(glue("wk{x-1}.dose_this_week"), glue("wk{x}.use_this_week")))
Y <- glue("wk{4:12}.relapse_this_week")

sl <- c("SL.glm", "SL.xgboost", "SL.earth", "SL.mean")

dat <- left_join(bup, dat)

baseline <- as.data.frame(model.matrix(reformulate(W), dat)[, -1])
W <- names(baseline)

dat <- cbind(baseline, dat[, c(A, unlist(L), Y)])

condA <- dat[, glue("wk{3:11}.use_this_week")] == 1
colnames(condA) <- A

condB <- dat[, glue("wk{2:10}.dose_this_week")] < 32
colnames(condB) <- A

dynamic <- dat
dynamic[, A] <- apply(condA & condB, 2, \(x) as.numeric(x), simplify = FALSE)

odtr <- dat
d <- readRDS("data/drv/optimal-rule.rds")
for (a in A) {
    odtr[[a]] <- d[, a]
}

constant <- dat
for (a in A) {
    constant[[a]] <- 0
}

lmtp <- function(data, tau, shifted) {
    stack <- "SL.mean"
    
    lmtp_tmle(
        data, 
        A, Y, W, L, 
        shifted = shifted, 
        outcome_type = "survival", 
        folds = 1, 
        learners_outcome = stack, 
        learners_trt = stack, 
        .learners_outcome_folds = 2, 
        .learners_trt_folds = 2
    )
}

fits <- list(
    dynamic = lmtp(dat, 11, dynamic), 
    odtr = lmtp(dat, 11, odtr), 
    constant = lmtp(dat, 11, constant)
)

# Counts for following a strategy
strategy_n <- function(x) {
    colnames(x) <- c("Wk. 3", 4:11)
    colSums(x, na.rm = TRUE) |>
        as_tibble(rownames = "time") |>
        pivot_wider(
            names_from = "time",
            values_from = "value"
        )
}

increased <- dat[, glue("wk{3:11}.dose_increase_this_week")] == 1

map_dfr(
    list(
        `\\emph{Total}` = fits$constant$density_ratios != 0, 
        `\\hspace{1em}$\\d_0^a$` = bup[, glue("wk{4:12}.relapse_this_week")] == 0 & 
            odtr[, glue("wk{3:11}.dose_increase_this_week")] == 1 & 
            fits$constant$density_ratios != 0, 
        `\\hspace{1em}$\\d1^a$` =  bup[, glue("wk{4:12}.relapse_this_week")] == 0 & 
            dynamic[, glue("wk{3:11}.dose_increase_this_week")] == 1 & 
            fits$constant$density_ratios != 0,
        `Increase`= fits$odtr$density_ratios != 0 & increased,
        `Constant`= fits$odtr$density_ratios != 0 & !increased,
        `\\emph{Total}` = fits$odtr$density_ratios != 0, 
        `Increase`= fits$dynamic$density_ratios != 0 & increased,
        `Constant`= fits$dynamic$density_ratios != 0 & !increased,
        `\\emph{Total}` = fits$dynamic$density_ratios != 0
    ), strategy_n, .id = "subset"
) |> 
    kbl(format = "latex", booktabs = TRUE, escape = FALSE) |>
    pack_rows("Constant dosing strategy", 1, 3) |> 
    pack_rows("Learning optimal strategy, $d_0$", 4, 6) |>
    pack_rows("Increase dose in response to use strategy, $d1$", 7, 9)

odtr <- readRDS("data/drv/survival-combined-odtr.rds")
simple <- readRDS("data/drv/survival-combined-dynamic.rds")
constant <- readRDS("data/drv/survival-combined-constant.rds")

y_odtr <- 
    tidy(odtr) |> 
    mutate(estimate = 1 - estimate, 
           conf.low = estimate + qnorm(0.025)*std.error, 
           conf.high = estimate - qnorm(0.025)*std.error)

y_simp <- 
    tidy(simple) |> 
    mutate(estimate = 1 - estimate, 
           conf.low = estimate + qnorm(0.025)*std.error, 
           conf.high = estimate - qnorm(0.025)*std.error)

y_cons <- 
    tidy(constant) |> 
    mutate(estimate = 1 - estimate, 
           conf.low = estimate + qnorm(0.025)*std.error, 
           conf.high = estimate - qnorm(0.025)*std.error)

y_d <- bind_rows(y_cons, y_simp, y_odtr)

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
# constant & \includegraphics[width=0.67in, height=0.17in]{file:////Volumes/CUMC - Seni/odtr/figures/pointrange_b9c649cb1cb3.pdf}\\
# dynamic & \includegraphics[width=0.67in, height=0.17in]{file:////Volumes/CUMC - Seni/odtr/figures/pointrange_b9c6340d77a8.pdf}\\
# odtr & \includegraphics[width=0.67in, height=0.17in]{file:////Volumes/CUMC - Seni/odtr/figures/pointrange_b9c61b128299.pdf}\\
# \bottomrule
# \end{tabular}
