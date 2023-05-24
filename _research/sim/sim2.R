library(odtr)
library(purrr)

source("_research/sim/gendata2.R")

id <- Sys.getenv("SGE_TASK_ID")
if (id == "undefined" || id == "") id <- 1

learners <- c("SL.nnet", "SL.earth", "SL.gam", "SL.bayesglm", "SL.glm.interaction", "SL.glm")

L <- list(c("W1", "W2"), c("W3"))
vars <- Vars$new(L = L, A = c("A1", "A2"), Y = "Y")

returns <- map_dfr(c(500, 1000, 1e4), function(n) {
    seed <- floor(runif(1, 1000, 1e5))
    set.seed(seed)
    tmp <- gendata(n)
    opt <- optimal_rule(tmp, vars, 1, "SL.mean", learners, "continuous")
    data.frame(id = id, 
               seed = seed, 
               n = n,
               psi = opt$psi, 
               conf.low = opt$conf.low, 
               conf.high = opt$conf.high)
})

saveRDS(returns, glue::glue("_research/sim/data/raw/sim_{id}.rds"))
