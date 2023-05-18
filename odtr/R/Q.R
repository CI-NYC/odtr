crossFitQ <- function(data, g, Npsem, learners, folds, 
                      outcome_type = c("binomial", "continuous"), maximize = TRUE) {
    tau <- length(Npsem$A)
    
    Q0 <- lapply(1:3, function(x) matrix(nrow = nrow(data), ncol = tau))
    m <- matrix(nrow = nrow(data), ncol = tau + 1)
    m[, tau + 1] <- data[[Npsem$Y]]
    
    A_opt <- matrix(nrow = nrow(data), ncol = tau)
    
    for (t in tau:1) {
        . <- data
        if (t == tau) {
            y <- Npsem$Y
            type <- match.arg(outcome_type)
        } else {
            y <- g("tmp_Yd_{t+1}")
            .[[y]] <- m[, t + 1]
            type <- "continuous"
        }
        vars <- Npsem$history("L", t + 1)
        a_r <- at_risk(., Npsem, t)
        
        valid1 <- valid0 <- .
        valid0[[Npsem$A[t]]] <- 0
        valid1[[Npsem$A[t]]] <- 1
browser()

        # .f <- as.formula(paste0(y, paste0("~ -1 + .^", ncol(.[a_r, vars]))))
        # .f <- as.formula(paste0(y, paste0("~ -1 + .^", 2)))
        .f <- as.formula(paste0(y, "~ -1 + ."))
        fit <- mlr3superlearner(.[a_r, ], .f, 
                                learners, 
                                type, 
                                10,
                                newdata = list(.[a_r, ], valid0[a_r, ], valid1[a_r, ]))

        Q0[[1]][a_r, t] <- fit$preds[[1]];  Q0[[1]][!a_r, t] <- 0
        Q0[[2]][a_r, t] <- fit$preds[[2]];  Q0[[2]][!a_r, t] <- 0
        Q0[[3]][a_r, t] <- fit$preds[[3]];  Q0[[3]][!a_r, t] <- 0
        
        # vars <- setdiff(setdiff(vars, Npsem$A[t]), Npsem$history("L", t))
        vars <- setdiff(vars, Npsem$A[t])
        . <- cbind(.[, vars, drop = F], 
                   tmp_pseudo_blip_D = transform(g, t, data[, Npsem$A, drop = F], A_opt, m, Q0[[1]], Q0[[2]], Q0[[3]]))

        mtilde <- mlr3superlearner(.[a_r, ], 
                                   tmp_pseudo_blip_D ~ -1 + ., 
                                   learners, 
                                   "continuous", 
                                   10,
                                   newdata = list(.[a_r, ]))$preds[[1]]

        if (maximize) {
            A_opt[a_r, t] <- ifelse(mtilde > 0, 1, 0) 
        } else {
            A_opt[a_r, t] <- ifelse(mtilde < 0, 1, 0)
        }
        
        . <- data
        .[[Npsem$A[t]]][a_r] <- A_opt[, t][a_r]
        m[a_r, t] <- predict(fit, .[a_r, ])
        m[!a_r, t] <- 0
    }
    
    list(A_opt = A_opt, m = m, Q_a = Q0[[1]])
}

transform <- function(g, t, A, Aopt, ytilde, QA, Q0, Q1) {
    tau <- ncol(QA)
    wts <- matrix(1, nrow = nrow(g), ncol = length(1:tau))
    
    wts[, t] <- A[, t]*(1 / g[, t]) + (1 - A[, t])*(1 / (1 - g[, t]))
    if (t < tau) {
        for (s in (t + 1):tau) {
            wts[, s] <- as.numeric(A[, s] == Aopt[, s]) / g[, s]
            wts[is.na(wts[, s]), s] <- 0
        }
    }

    wts <- apply(wts, 2, function(x) pmin(x, quantile(x, 0.99, na.rm = T)))
    r <- ratio_sdr(wts, t, tau)
    m <- ytilde[, (t + 1):(tau + 1), drop = FALSE] - QA[, t:tau, drop = FALSE]
    eval1 <- rowSums(r * m, na.rm = TRUE) + Q1[, t] - Q0[, t]
    # m <- rep(0, nrow(wts))
    # for (s in t:tau) {
    #     m <- m + (apply(wts[, t:s, drop = F], 1, prod) *
    #                   (ytilde[, s + 1, drop = FALSE] - QA[, s, drop = FALSE]))
    # }
    # m <- m + (Q1[, t, drop = FALSE] - Q0[, t, drop = FALSE])
    # eval2 <- as.vector(m)
    eval1
}

ratio_sdr <- function(g, t, tau) {
    if ((t + 1) > tau) return(g[, tau])
    t(apply(g[, t:tau, drop = FALSE], 1, cumprod))
}
