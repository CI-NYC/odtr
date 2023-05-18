crossFitg0 <- function(data, Npsem, learners, folds) {
    g0 <- matrix(nrow = nrow(data), ncol = length(Npsem$A))
    bnd <- 5 / sqrt(nrow(data)) / log(nrow(data))

    for (t in 1:length(Npsem$A)) {
        for (v in seq_along(folds)) {
            train <- origami::training(data, folds[[v]])
            valid <- origami::validation(data, folds[[v]])
            
            # Need to handle risk here
            train <- train[at_risk(train, Npsem, t), ]
            risk <- at_risk(valid, Npsem, t)
            valid <- valid[risk, ]
        
            fit <- crossFit(train, list(valid), Npsem$A[t], Npsem$history("A", t), "binomial", learners)
            g0[folds[[v]]$validation_set[risk], t] <- pmax(pmin(fit[[1]], 1 - bnd), bnd)
        } 
    }
    g0
}
