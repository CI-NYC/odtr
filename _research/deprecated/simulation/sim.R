library(sl3)
library(ltmle)
library(lmtp)

devtools::load_all("pkg")

id <- Sys.getenv("SGE_TASK_ID")
if (id == "undefined" || id == "") id <- 1
 
seed <- round(runif(1, 1, 1e5))
seed <- 356
set.seed(seed)

ObsData <- sampleLuedtke2015(1e5)
sem <- Npsem$new(paste0("W", 1:4), A = c("A1", "A2"), Y = "Y")

Sl <- make_learner(
    Lrnr_sl,
    learners = make_learner_stack(
        Lrnr_glm, Lrnr_rpart, Lrnr_earth, Lrnr_xgboost
    ),
    # metalearner = make_learner("Lrnr_nnls", convex = TRUE),
    keep_extra = FALSE
)

# Qbar0(
#     ObsData[, c("A1", "A2")],
#     ObsData[, paste0("W", 1:4)]
# ) |> head()

d <- odtr(ObsData, sem, 1, Lrnr_mean$new(), Lrnr_xgboost$new(), "binomial")

table(
    as.character(d$A2),
    unlist(lapply(strsplit(ObsData$d, split = ""), \(x) x[2]))
)

mean(as.character(d$A2) == unlist(lapply(strsplit(ObsData$d, split = ""), \(x) x[2])))

table(
    paste0(d$A1, d$A2), 
    ObsData$d
)

mean(paste0(d$A1, d$A2) == ObsData$d)

Odtr_dat <- ObsData
Odtr_dat$A1 <- d$A1
Odtr_dat$A2 <- d$A2

Ey_d <- lmtp_tmle(
    ObsData, c("A1", "A2"), "Y", paste0("W", 1:4),
    shifted = Odtr_dat, folds = 1,
    learners_outcome = Lrnr_xgboost$new(), learners_trt = Lrnr_xgboost$new()
)

write.csv(
    data.frame(
        seed = seed,
        n = 1000,
        psi = Ey_d$treatment$estimate,
        std_error = Ey_d$treatment$std.dev, 
        conf_low = Ey_d$treatment$CI[1], 
        conf_high = Ey_d$treatment$CI[2]
    ),
    glue::glue("data/simulation/{id}.csv"),
    row.names = FALSE
)

quit("no")
