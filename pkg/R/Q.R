crossFitQ0 <- function(data, y, a, vars, t, Npsem, learners, folds, outcome_type = c("binomial", "continuous")) {
    Q0 <- matrix(nrow = nrow(data), ncol = 3)
    fits <- list()
    for (v in seq_along(folds)) {
        train <- origami::training(data, folds[[v]])
        valid0 <- valid1 <- valid <- 
            origami::validation(data, folds[[v]])
        
        risk_train <- at_risk(train, Npsem, t)
        risk_valid <- at_risk(valid, Npsem, t)
        
        valid0[[a]] <- 0
        valid1[[a]] <- 1

        preds <- 
            crossFit(train[risk_train, ], list(valid[risk_valid, ], valid0[risk_valid, ], valid1[risk_valid, ]), 
                     y, vars, NULL, match.arg(outcome_type), learners, TRUE)

        fits[[v]] <- preds$fit
        
        Q0[folds[[v]]$validation_set[risk_valid], 1] <- preds$preds[[1]]
        Q0[folds[[v]]$validation_set[risk_valid], 2] <- preds$preds[[2]]
        Q0[folds[[v]]$validation_set[risk_valid], 3] <- preds$preds[[3]]
        
        Q0[folds[[v]]$validation_set[!risk_valid], 1] <- 0
        Q0[folds[[v]]$validation_set[!risk_valid], 2] <- 0
        Q0[folds[[v]]$validation_set[!risk_valid], 3] <- 0
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

        risk_train <- at_risk(train, Npsem, t)
        risk_valid <- at_risk(valid, Npsem, t)

        A <- train[[Npsem$A[t]]][risk_train]
        Y <- train[[Npsem$Y]][risk_train]
        
        H <- g_train[risk_train, t]
        H <- A*(1 / H) + (1 - A)*(1 / (1 - H))
        
        D1 <- (H * (Y - Q0_train[risk_train, "A=a"])) + Q0_train[risk_train, "A=1"] - Q0_train[risk_train, "A=0"]
        
        # if (t == 2) {
        #     weights <- H
        # } else {
        #     weights <- rep(1, length(A))
        # }
        
        weights <- NULL
        
        train[risk_train, "tmp_psuedo_D_blip"] <- D1

        Qv[folds[[v]]$validation_set[risk_valid], 1] <- 
            crossFit(train[risk_train, ], list(valid[risk_valid, ]), "tmp_psuedo_D_blip", 
                     Npsem$history("A", t), weights, "continuous", learners)[[1]]

        #Qv[folds[[v]]$validation_set[!risk_valid], 1] <- 0
    }
    Qv
}

addYd_t <- function(data, Npsem, OdtrA_t, Q0_tFits, folds, t) {
    Yd_t <- matrix(nrow = nrow(data), ncol = 1)
    for (v in seq_along(folds)) {
        valid <- origami::validation(data, folds[[v]])
        risk <- at_risk(valid, Npsem, t)
        valid[[Npsem$A[t]]][risk] <- origami::validation(OdtrA_t, folds[[v]])[risk]
        
        Yd_t[folds[[v]]$validation_set[risk], 1] <- predictt(Q0_tFits[[v]], valid[risk, ])
        Yd_t[folds[[v]]$validation_set[!risk], 1] <- 0
    } 

    data[[g("tmp_Yd_{t}")]] <- Yd_t[, 1]
    data
}
