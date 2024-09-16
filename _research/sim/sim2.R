library(odtr)
library(purrr)
library(mlr3extralearners)

source("_research/sim/gendata2.R")

id <- Sys.getenv("SGE_TASK_ID")
if (id == "undefined" || id == "") id <- 1

learners <- c("mean", "glm", "earth", "lightgbm")

A <- c("A1", "A2")
L <- list(c("W1", "W2"), c("W3"))
Y <- "Y"

returns <- map_dfr(c(500, 1000, 1e4), function(n) {
    seed <- floor(runif(1, 1000, 1e5))
    set.seed(seed)
    
    tmp <- gendata(n)
    
    opt <- odtr(data = tmp, 
                trt = A, 
                outcome = Y, 
                baseline = NULL, 
                time_varying = L, 
                learners_trt = learners, 
                learners_outcome = learners, 
                learners_rule = learners, 
                folds = 10, 
                outcome_type = "continuous")
    
    data.frame(id = id, 
               seed = seed, 
               n = n,
               psi = opt$psi, 
               conf.low = opt$conf.low, 
               conf.high = opt$conf.high)
})

saveRDS(returns, glue::glue("_research/sim/data/raw/sim_{id}.rds"))
