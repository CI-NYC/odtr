source("_research/leudtke-simulation.R")
source("R/npsem.R")
source("R/folds.R")
source("R/crossfit.R")
source("R/g.R")
source("R/Q.R")

library(sl3)

dat <- sim.data.mtp(1e3)

np <- Npsem$new(paste0("W", 1:4), A = c("A1", "A2"), Y = "Y")
folds <- make_folds(dat, 10)

# Fit propensity at all time points
g0 <- crossFitg0(dat, np, Lrnr_glm$new(), folds)

Sl <- make_learner(
    Lrnr_sl,
    learners = make_learner_stack(
        Lrnr_glm, Lrnr_mean
    ),
    metalearner = make_learner("Lrnr_nnls", convex = TRUE),
    keep_extra = FALSE
)

# Fit the first (really last) outcome regression
Q0_2 <- crossFitQ0(dat, np$Y, np$A[2], np$history("Y"), Sl, folds, "binomial")

# Define the DR blip function, regress and predict on the set V
Qv_2 <- crossFitQv(dat, g0, Q0_2$Q0, 2, np, Sl, folds)

# Identify optimal treatment decision for time 2
OdtrA_2 <- ifelse(Qv_2[, 1] > 0, 1, 0)

# Need estimates of Y_d2. This requires saving the fits from Q0_2 and predicting using OdtrA_2? 
Yd_2 <- matrix(nrow = nrow(dat), ncol = 1)
for (v in seq_along(folds)) {
    valid <- origami::validation(dat, folds[[v]])
    valid[[np$A[2]]] <- origami::validation(OdtrA_2, folds[[v]])
    Yd_2[folds[[v]]$validation_set, 1] <- predictt(Q0_2$fits[[v]], valid)
} 

dat[["tmp_Yd_2"]] <- Yd_2[, 1]

# Fit the next outcome regression with the pseudo outcome
Q0_1 <- crossFitQ0(dat, "tmp_Yd_2", np$A[1], np$history("L", 2), Sl, folds, "continuous")

# Define DR blip, regress and predict
Qv_1 <- crossFitQv(dat, g0, Q0_1$Q0, 1, np, Sl, folds)

# Identify optimal treatment decision for time 1
OdtrA_1 <- ifelse(Qv_1[, 1] > 0, 1, 0)

# evaluate Y_d ------------------------------------------------------------

Odtr_dat <- dat
Odtr_dat$A1 <- OdtrA_1
Odtr_dat$A2 <- OdtrA_2

library(lmtp)

Ey_d <- lmtp_tmle(dat, c("A1", "A2"), "Y", paste0("W", 1:4), shifted = Odtr_dat)
