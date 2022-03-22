crossFit <- function(train, valid, y, x, weights = NULL, 
                     type = c("binomial", "continuous"), learners, returnFit = FALSE) {
    fit <- regress(train, y, x, weights, match.arg(type), learners)
    preds <- lapply(valid, function(x) predictt(fit, x))
    
    if (isFALSE(returnFit)) {
        return(preds)
    }
    
    list(
        preds = preds, 
        fit = fit
    )
}

regress <- function(train, y, x, weights, type, learners) {
    if (!is.null(weights)) {
        train[, "tmp_odtr_weights"] <- weights
        weights <- "tmp_odtr_weights"
    }
    
    task <- sl3::sl3_Task$new(
        data = train,
        covariates = x,
        outcome = y,
        outcome_type = type, 
        weights = weights
    )

    fit <- learners$train(task)
    fit
}

predictt <- function(fit, data) {
    task <- sl3::sl3_Task$new(
        data = data,
        covariates = fit$training_task$nodes$covariates
    )
    
    fit$predict(task)
}
