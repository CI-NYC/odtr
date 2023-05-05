crossFitg0 <- function(data, Npsem, learners, folds) {
    g0 <- matrix(nrow = nrow(data), ncol = length(Npsem$A))
    for (t in 1:length(Npsem$A)) {
        for (v in seq_along(folds)) {
            train <- origami::training(data, folds[[v]])
            valid <- origami::validation(data, folds[[v]])
            
            # Need to handle risk here
            train <- train[at_risk(train, Npsem, t), ]
            risk <- at_risk(valid, Npsem, t)
            valid <- valid[risk, ]
            
            fit <- mlr3superlearner::mlr3superlearner(train[, c(Npsem$history("A", t), Npsem$A[t])], 
                                                      Npsem$A[t], 
                                                      learners, 
                                                      "binomial", 
                                                      newdata = list(valid))
            
            g0[folds[[v]]$validation_set[risk], t] <- fit$preds[[1]]
        } 
    }
    g0
}
