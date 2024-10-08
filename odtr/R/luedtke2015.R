#' @export
sampleLuedtke2015 <- function(n) {
    W1 <- 2 * runif(n) - 1
    W2 <- 2 * runif(n) - 1
    A1 <- rbinom(n, 1, 1 / 2)
    u3 <- runif(n)
    W3 <- (2 * u3 - 1) * (1.25 * A1 + 0.25)
    W3.0 <- (2 * u3 - 1) * (0 + 0.25)
    W3.1 <- (2 * u3 - 1) * (1.25 + 0.25)
    u4 <- runif(n)
    W4 <- (2 * u4 - 1) * (1.25 * A1 + 0.25)
    W4.0 <- (2 * u4 - 1) * (0 + 0.25)
    W4.1 <- (2 * u4 - 1) * (1.25 + 0.25)
    A2 <- rbinom(n, 1, 1 / 2)
    
    A <- data.frame(A1 = A1, A2 = A2)
    W <- data.frame(
        W1 = W1,
        W2 = W2,
        W3 = W3,
        W4 = W4
    )
    
    u <- runif(n)
    Y <- as.numeric(u < Qbar0(A, W))
    
    Y.00 <-
        Qbar0(
            data.frame(A1 = rep(0, n), A2 = rep(0, n)),
            data.frame(
                W1 = W1,
                W2 = W2,
                W3 = W3.0,
                W4 = W4.0
            )
        )
    
    Y.01 <-
        Qbar0(
            data.frame(A1 = rep(0, n), A2 = rep(1, n)),
            data.frame(
                W1 = W1,
                W2 = W2,
                W3 = W3.0,
                W4 = W4.0
            )
        )
    
    Y.10 <-
        Qbar0(
            data.frame(A1 = rep(1, n), A2 = rep(0, n)),
            data.frame(
                W1 = W1,
                W2 = W2,
                W3 = W3.1,
                W4 = W4.1
            )
        )
    
    Y.11 <-
        Qbar0(
            data.frame(A1 = rep(1, n), A2 = rep(1, n)),
            data.frame(
                W1 = W1,
                W2 = W2,
                W3 = W3.1,
                W4 = W4.1
            )
        )
    
    cf_probs <- data.frame(Y.00, Y.01, Y.10, Y.11)
    
    data.frame(
        W1 = W1,
        W2 = W2,
        W3 = W3,
        W4 = W4,
        A1 = A1,
        A2 = A2,
        d = c("00", "01", "10", "11")[apply(cf_probs, 1, which.max)],
        Y = Y
    )
}

Qbar0 <- function(A, W) {
    0.4 + (-0.4 * A$A1 - A$A2 * W$W3 - 1.5 * A$A1 * sign(W$W1) + 0.08 * A$A1 * A$A2 + A$A2 * (W$W3)^2 - 1.5 * A$A1 * W$W1 - 0.5 * A$A1 * (W$W2) ^ 2 - 0.1 * A$A2) * 0.069
}
