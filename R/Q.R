# Assume only a 2-time point scenario...
# Assuming treatment is binary
crossFitQ0 <- function(data, y, a, vars, learners, folds, outcome_type = c("binomial", "continuous")) {
    Q0 <- matrix(nrow = nrow(data), ncol = 3)
    fits <- list()
    for (v in seq_along(folds)) {
        train <- origami::training(data, folds[[v]])
        valid0 <- valid1 <- valid <- 
            origami::validation(data, folds[[v]])
        
        valid0[[a]] <- 0
        valid1[[a]] <- 1

        preds <- 
            crossFit(train, list(valid, valid0, valid1), y, vars, match.arg(outcome_type), learners, TRUE)

        fits[[v]] <- preds$fit
        
        Q0[folds[[v]]$validation_set, 1] <- preds$preds[[1]]
        Q0[folds[[v]]$validation_set, 2] <- preds$preds[[2]]
        Q0[folds[[v]]$validation_set, 3] <- preds$preds[[3]]
    } 
    colnames(Q0) <- c("A=a", "A=0", "A=1")
    list(
        Q0 = Q0, 
        fits = fits
    )
}

crossFitQv <- function(data, g, Q0, t, Npsem, learners, folds) {
    Qv <- matrix(nrow = nrow(data), ncol = 1)
    for (v in seq_along(folds)) {
        train <- origami::training(data, folds[[v]])
        valid <- origami::validation(data, folds[[v]])
        g_train <- origami::training(g, folds[[v]])
        Q0_train <- origami::training(Q0, folds[[v]])

        A <- train[[Npsem$A[t]]]
        Y <- train[[Npsem$Y]]
        D1 <- (2 * A - 1) / g_train[, t] * (Y - Q0_train[, "A=a"]) + Q0_train[, "A=1"] - Q0_train[, "A=0"]
        train[, "tmp_psuedo_D_blip"] <- D1

        # Need to ask about what variables can be included here?
        Qv[folds[[v]]$validation_set, 1] <- 
            crossFit(train, list(valid), "tmp_psuedo_D_blip", c(Npsem$W, unlist(Npsem$L)), "continuous", learners)[[1]]
    }
    Qv
}

addYd_t <- function(data, npsem, OdtrA_t, Q0_tFits, folds, t) {
    Yd_t <- matrix(nrow = nrow(data), ncol = 1)
    for (v in seq_along(folds)) {
        valid <- origami::validation(data, folds[[v]])
        valid[[npsem$A[t]]] <- origami::validation(OdtrA_t, folds[[v]])
        Yd_t[folds[[v]]$validation_set, 1] <- predictt(Q0_tFits[[v]], valid)
    } 
    
    data[[g("tmp_Yd_{t}")]] <- Yd_t[, 1]
    data
}
