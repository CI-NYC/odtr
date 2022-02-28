crossFitg0 <- function(data, Npsem, learners, folds) {
    g0 <- matrix(nrow = nrow(data), ncol = length(Npsem$A))
    for (t in 1:length(Npsem$A)) {
        for (v in seq_along(folds)) {
            train <- origami::training(data, folds[[v]])
            valid <- origami::validation(data, folds[[v]])
            
            g0[folds[[v]]$validation_set, t] <- 
                crossFit(train, list(valid), Npsem$A[t], Npsem$history("A", t), "binomial", learners)[[1]]
        } 
    }
    g0
}
