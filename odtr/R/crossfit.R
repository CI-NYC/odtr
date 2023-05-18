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
    family <- ifelse(type == "binomial", binomial(), gaussian())
    SuperLearner::SuperLearner(
        train[[y]], train[, x, drop = F], family = family[[1]], SL.library = learners,
        method = "method.CC_LS",
        env = environment(SuperLearner::SuperLearner)
    )
}

predictt <- function(fit, data) {
    predict(fit, data[, fit$varNames, drop = F])$pred[, 1]
}
