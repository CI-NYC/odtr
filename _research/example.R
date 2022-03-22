library(sl3)
library(lmtp)

devtools::load_all("pkg")

dat <- sampleLuedtke2015(1e3)
sem <- Npsem$new(paste0("W", 1:4), A = c("A1", "A2"), Y = "Y")

Sl <- make_learner(
    Lrnr_sl,
    learners = make_learner_stack(
        Lrnr_glm, Lrnr_mean, Lrnr_xgboost
    ),
    metalearner = make_learner("Lrnr_nnls", convex = TRUE),
    keep_extra = FALSE
)

optimal <- odtr(dat, sem, 1, Lrnr_glm$new(), Sl, "binomial")

Odtr_dat <- dat
Odtr_dat$A1 <- optimal$A1
Odtr_dat$A2 <- optimal$A2

Ey_d <- lmtp_tmle(
    dat, c("A1", "A2"), "Y", paste0("W", 1:4), 
    shifted = Odtr_dat, folds = 1,
    learners_outcome = Sl, 
    learners_trt = Sl
)

Ey_d
