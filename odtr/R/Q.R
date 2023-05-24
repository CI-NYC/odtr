crossFitQ <- function(data, g, Vars, learners, folds, 
                      outcome_type = c("binomial", "continuous"), maximize = TRUE) {
    tau <- length(Vars$A)
    
    QA <- matrix(nrow = nrow(data), ncol = tau)
    m <- matrix(nrow = nrow(data), ncol = tau + 1)
    m[, tau + 1] <- data[[Vars$Y]]

    A_opt <- matrix(nrow = nrow(data), ncol = tau)

    for (v in seq_along(folds)) {
        Q0 <- lapply(1:3, function(x) matrix(nrow = nrow(data), ncol = tau))
        mt <- matrix(nrow = nrow(data), ncol = tau + 1)
        mt[, tau + 1] <- data[[Vars$Y]]
        A_optt <- matrix(nrow = nrow(data), ncol = tau)
        
        for (t in tau:1) {
            vars <- Vars$history("L", t + 1)
            . <- data
            train <- tr(., folds[[v]])
            valid <- vl(., folds[[v]])
            
            ti <- folds[[v]]$training
            vi <- folds[[v]]$validation
            
            if (t == tau) {
                y <- Vars$Y
                type <- match.arg(outcome_type)
            } else {
                y <- g("tmp_Yd_{t+1}")
                train[[y]] <- tr(mt, folds[[v]])[, t + 1]
                type <- "continuous"
            }
            
            art <- at_risk(train, Vars, t)
            arv <- at_risk(valid, Vars, t)
            
            train1 <- train0 <- train
            train0[[Vars$A[t]]] <- 0
            train1[[Vars$A[t]]] <- 1

            fit <- crossFit(train[art, ], 
                            list(train[art, ], train0[art, ], train1[art, ], valid[arv, ]), 
                            y, vars, type, learners, TRUE)
            
            Q0[[1]][ti, ][art, t] <- fit$preds[[1]]
            Q0[[1]][ti, ][!art, t] <- 0
            Q0[[2]][ti, ][art, t] <- fit$preds[[2]]
            Q0[[2]][ti, ][!art, t] <- 0
            Q0[[3]][ti, ][art, t] <- fit$preds[[3]]
            Q0[[3]][ti, ][!art, t] <- 0
            QA[vi, ][arv, t] <- fit$preds[[4]]
            QA[vi, ][!arv, t] <- 0

            vars <- setdiff(vars, Vars$A[t])
            train <- cbind(train[, vars, drop = F], 
                           tmp_pseudo_blip_D = 
                               transform(tr(g, folds[[v]]), 
                                         t, 
                                         train[, Vars$A, drop = F], 
                                         tr(A_optt, folds[[v]]), 
                                         tr(m, folds[[v]]), 
                                         tr(Q0[[1]], folds[[v]]), 
                                         tr(Q0[[2]], folds[[v]]), 
                                         tr(Q0[[3]], folds[[v]])))
            
            mtilde <- crossFit(train[art, ], 
                               list(valid[arv, ], train[art, ]), 
                               "tmp_pseudo_blip_D", 
                               vars, 
                               "continuous", 
                               learners)
            
            if (maximize) {
                opt_trt <- function(x) ifelse(x > 0, 1, 0) 
            } else {
                opt_trt <- function(x) ifelse(x < 0, 1, 0)
            }
            
            A_opt[vi, ][arv, t] <- opt_trt(mtilde[[1]])
            A_optt[ti, ][art, t] <- opt_trt(mtilde[[2]])
            
            . <- data
            .[[Vars$A[t]]][vi][arv] <- A_opt[vi, t][arv]
            .[[Vars$A[t]]][ti][art] <- A_optt[ti, t][art]
            
            mt[ti, ][art, t] <- predictt(fit$fit, .[ti, ][art, ])
            mt[ti, ][!art, t] <- 0
            
            m[vi, ][arv, t] <- predictt(fit$fit, .[vi, ][arv, ])
            m[vi, ][!arv, t] <- 0
        }
    }
    list(A_opt = A_opt, m = m, Q_a = QA)
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
