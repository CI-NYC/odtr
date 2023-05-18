dat <- odtr_sim(1e4)
sem <- Npsem$new(paste0("W", 1:3), A = c("A1", "A2"), Y = "Y")
opt <- odtr(dat, sem, 1, c("glm", "nnet"), c("glm", "nnet"), "binomial")

d <- dat
tr <- strsplit(dat$d, "")
d$A1 <- as.numeric(purrr::map_chr(tr, 1))
d$A2 <- as.numeric(purrr::map_chr(tr, 2))

table(d[, c("A1", "A2")])
table(opt$A_opt)
mean(dat$Ymax)

library(devtools)
library(future)
library(mlr3extralearners)

load_all("odtr")

source("_research/sim/gendata.R")

truth <- {
    tmp <- gendata(1e7)
    mean(tmp$Ymax)
}

learners <- c("glm", "xgboost", "nnet")

L <- list(c("W1"), 
          c("L1"))
vars <- Npsem$new(L = L, A = c("A1", "A2"), Y = "Y")

tmp <- gendata(1000)
opt <- odtr(tmp, vars, 1, learners, learners, "binomial")

d <- tmp
tr <- strsplit(tmp$d, "")
d$A1 <- as.numeric(purrr::map_chr(tr, 1))
d$A2 <- as.numeric(purrr::map_chr(tr, 2))

table(d[, c("A1", "A2")])
table(opt$A_opt)
