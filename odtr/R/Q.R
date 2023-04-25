crossFitQ <- function(data, g, Npsem, learners, folds, outcome_type = c("binomial", "continuous"), maximize = TRUE) {
    tau <- length(Npsem$A)
    
    Q0 <- lapply(1:3, function(x) matrix(nrow = nrow(data), ncol = tau))
    m <- matrix(nrow = nrow(data), ncol = tau + 1)
    m[, tau + 1] <- data[[Npsem$Y]]
    
    A_opt <- matrix(nrow = nrow(data), ncol = tau)
    
    for (t in tau:1) {
        . <- data
        if (t == tau) {
            y <- Npsem$Y
            vars <- Npsem$history("Y")
            type <- match.arg(outcome_type)
        } else {
            y <- g("tmp_Yd_{t+1}")
            .[[y]] <- m[, t + 1]
            vars <- Npsem$history("L", t + 1)
            type <- "continuous"
        }
        
        a_r <- at_risk(., Npsem, t)
        
        valid1 <- valid0 <- .
        valid0[[Npsem$A[t]]] <- 0
        valid1[[Npsem$A[t]]] <- 1
        
        fit <- crossFit(
            .[a_r, ], 
            list(.[a_r, ], valid0[a_r, ], valid1[a_r, ]), 
            y, vars, type, learners, TRUE
        )
        
        Q0[[1]][a_r, t] <- fit$preds[[1]];  Q0[[1]][!a_r, t] <- 0
        Q0[[2]][a_r, t] <- fit$preds[[2]];  Q0[[2]][!a_r, t] <- 0
        Q0[[3]][a_r, t] <- fit$preds[[3]];  Q0[[3]][!a_r, t] <- 0
        
        . <- cbind(.[, c(Npsem$W, Npsem$L[[t]])], 
                   tmp_pseudo_blip_D = transform(g, t, data[, Npsem$A], A_opt, m, Q0[[1]], Q0[[2]], Q0[[3]]))
        
        mtilde <- crossFit(.[a_r, ],
                           list(.[a_r, ]),
                           "tmp_pseudo_blip_D",
                           c(Npsem$W, Npsem$L[[t]]),
                           "continuous",
                           learners)[[1]]
        
        if (maximize) {
            A_opt[, t] <- ifelse(mtilde > 0, 1, 0) 
        } else {
            A_opt[, t] <- ifelse(mtilde < 0, 1, 0)
        }
        
        . <- data
        .[[Npsem$A[t]]][a_r] <- A_opt[, t][a_r]
        m[a_r, t] <- predictt(fit$fit, .[a_r, ])
        m[!a_r, t] <- 0
    }
    
    A_opt
}

transform <- function(g, t, A, Aopt, ytilde, QA, Q0, Q1) {
    tau <- ncol(QA)
    wts <- matrix(1, nrow = nrow(g), ncol = length(1:tau))
    
    wts[, t] <- A[, t]*(1 / g[, t]) + (1 - A[, t])*(1 / (1 - g[, t]))
    if (t < tau) {
        for (s in (t + 1):tau) {
            wts[, s] <- as.numeric(A[, s] == Aopt[, s]) / g[, s]
        }
    }

    wts <- apply(wts, 1, prod)
    m <- (ytilde[, (t + 1):(tau + 1), drop = FALSE] - QA[, t:tau, drop = FALSE]) + 
        Q1[, t:tau, drop = FALSE] - Q0[, t:tau, drop = FALSE]
    
    rowSums(wts * m)
}
