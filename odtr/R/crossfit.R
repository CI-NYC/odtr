crossFit <- function(train, valid, y, x,
                     type = c("binomial", "continuous"), 
                     learners, returnFit = FALSE) {
    fit <- regress(train, y, x, match.arg(type), learners)
    preds <- lapply(valid, function(x) predictt(fit, x))
    
    if (isFALSE(returnFit)) {
        return(preds)
    }
    
    list(
        preds = preds, 
        fit = fit
    )
}

regress <- function(train, y, x, type, learners) {
    # if (!is.null(weights)) {
    #     train[, "tmp_odtr_weights"] <- weights
    #     weights <- "tmp_odtr_weights"
    # }
    
    family <- ifelse(type == "binomial", binomial(), gaussian())
    
    fit <- SuperLearner::SuperLearner(
        train[[y]], train[, x], family = family[[1]], SL.library = learners,
        # method = "method.CC_LS",
        method = "method.NNLS",
        env = environment(SuperLearner::SuperLearner)
    )
    
    # task <- sl3::sl3_Task$new(
    #     data = train,
    #     covariates = x,
    #     outcome = y,
    #     outcome_type = type, 
    #     weights = weights
    # )
    # 
    # fit <- learners$train(task)
    fit
}

predictt <- function(fit, data) {
    predict(fit, data[, fit$varNames])$pred[, 1]
    # task <- sl3::sl3_Task$new(
    #     data = data,
    #     covariates = fit$training_task$nodes$covariates
    # )
    # 
    # fit$predict(task)
}
