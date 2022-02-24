crossFit <- function(train, valid, y, x, type = c("binomial", "continuous"), learners) {
    fit <- regress(train, y, x, match.arg(type), learners)
    predictt(fit, valid)
}

regress <- function(train, y, x, type, learners) {
    task <- sl3::sl3_Task$new(
        data = train,
        covariates = x,
        outcome = y,
        outcome_type = type
    )
    
    learners$train(task)
}

predictt <- function(fit, data) {
    task <- sl3::sl3_Task$new(
        data = data,
        covariates = fit$training_task$nodes$covariates
    )
    
    fit$predict(task)
}
