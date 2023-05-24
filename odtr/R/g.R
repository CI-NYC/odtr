crossFitg0 <- function(data, Vars, learners, folds) {
    g0 <- matrix(nrow = nrow(data), ncol = length(Vars$A))
    bnd <- 5 / sqrt(nrow(data)) / log(nrow(data))

    for (t in 1:length(Vars$A)) {
        for (v in seq_along(folds)) {
            train <- origami::training(data, folds[[v]])
            valid <- origami::validation(data, folds[[v]])
            
            train <- train[at_risk(train, Vars, t), ]
            risk <- at_risk(valid, Vars, t)
            valid <- valid[risk, ]
        
            fit <- crossFit(train, list(valid), Vars$A[t], Vars$history("A", t), "binomial", learners)
            g0[folds[[v]]$validation_set[risk], t] <- pmax(pmin(fit[[1]], 1 - bnd), bnd)
        } 
    }
    g0
}
