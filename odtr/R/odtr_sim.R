#' @export
odtr_sim <- function(n) {
    W1 <- rbinom(n, 1, 0.75)
    W2 <- rbinom(n, 1, 0.5)
    W3 <- rbinom(n, 1, 0.35)
    A1 <- rbinom(n, 1, 1 / 2)
    A2 <- rbinom(n, 1, 1 / 2)
    
    A <- data.frame(A1 = A1, A2 = A2)
    W <- data.frame(
        W1 = W1,
        W2 = W2,
        W3 = W3
    )
    
    u <- runif(n)
    Y <- as.numeric(u < Qbar0(A, W))
    
    Y.00 <-
        Qbar0(
            data.frame(A1 = rep(0, n), A2 = rep(0, n)),
            data.frame(
                W1 = W1,
                W2 = W2,
                W3 = W3
            )
        )
    
    Y.01 <-
        Qbar0(
            data.frame(A1 = rep(0, n), A2 = rep(1, n)),
            data.frame(
                W1 = W1,
                W2 = W2,
                W3 = W3
            )
        )
    
    Y.10 <-
        Qbar0(
            data.frame(A1 = rep(1, n), A2 = rep(0, n)),
            data.frame(
                W1 = W1,
                W2 = W2,
                W3 = W3
            )
        )
    
    Y.11 <-
        Qbar0(
            data.frame(A1 = rep(1, n), A2 = rep(1, n)),
            data.frame(
                W1 = W1,
                W2 = W2,
                W3 = W3
            )
        )

    cf_probs <- data.frame(Y.00, Y.01, Y.10, Y.11)
    
    data.frame(
        W1 = as.numeric(W1),
        W2 = as.numeric(W2),
        W3 = as.numeric(W3),
        A1 = as.numeric(A1),
        A2 = as.numeric(A2),
        d = c("00", "01", "10", "11")[apply(cf_probs, 1, which.max)],
        Y = as.numeric(Y), 
        Ymax = apply(cf_probs, 1, max)
    )
}

Qbar0 <- function(A, W) {
    plogis(0.4 - 0.1*A$A1 + 0.2*A$A2 - 0.08*A$A2*W$W1 - 0.3*A$A1*W$W1 + 0.75*A$A1*W$W2)
}
