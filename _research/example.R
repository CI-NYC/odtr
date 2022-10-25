library(lmtp)

devtools::load_all("odtr")

dat <- sampleLuedtke2015(1e3)
sem <- Npsem$new(paste0("W", 1:4), A = c("A1", "A2"), Y = "Y")

optimal <- odtr(dat, sem, 1, "SL.glm", c("SL.glm", "SL.mean", "SL.xgboost"), "binomial")

Odtr_dat <- dat
Odtr_dat$A1 <- optimal$A1
Odtr_dat$A2 <- optimal$A2

Ey_d <- lmtp_tmle(
    dat, c("A1", "A2"), "Y", paste0("W", 1:4), 
    shifted = Odtr_dat, folds = 1,
    learners_outcome = c("SL.glm", "SL.mean", "SL.xgboost"), 
    learners_trt = c("SL.glm", "SL.mean", "SL.xgboost")
)

Ey_d
