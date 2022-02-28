crossFit <- function(train, valid, y, x, type = c("binomial", "continuous"), learners, returnFit = FALSE) {
    fit <- regress(train, y, x, match.arg(type), learners)
    preds <- lapply(valid, \(x) predictt(fit, x))
    
    if (isFALSE(returnFit)) {
        return(preds)
    }
    
    list(
        preds = preds, 
        fit = fit
    )
}

regress <- function(train, y, x, type, learners) {
    task <- sl3::sl3_Task$new(
        data = train,
        covariates = x,
        outcome = y,
        outcome_type = type
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
