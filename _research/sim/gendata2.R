Q2 <- function(A1, A2, W1, W2, W3) {
    0.4 + -0.4*A1 - A2*W3 - 4*A1*W1 + 0.08*A1*A2 + A2*W3 - 4*A1*W1 - 2*A1*W2 - 0.1*A2 + 1.5*W1
}

Q1 <- function(A1, A2, W1, W2) {
    0.4 + -0.4*A1 - 4*A1*W1 + 0.08*A1*A2 - 4*A1*W1 - 2*A1*W2 - 0.1*A2 + 1.5*W1
}

gendata <- function(n) {
    W1 <- 2*runif(n) - 1
    W2 <- 2*runif(n) - 1
    A1 <- rbinom(n, 1, plogis(0.5 - 1.3*W1 + 0.4*W2))
    u3 <- runif(n)
    W3 <- (2*u3 - 1)*(1.25*A1 + 0.25)
    W3.0 <- (2*u3 - 1)*(0 + 0.25)
    W3.1 <- (2*u3 - 1)*(1.25 + 0.25)
    A2 <- rbinom(n, 1, plogis(0.5 + 0.4*A1 - 1.5*W3))
    
    u <- rnorm(n)
    Y <- Q2(A1, A2, W1, W2, W3) + u

    Y.00 <- Q2(rep(0, n), rep(0, n), W1, W2, W3.0)
    Y.01 <- Q2(rep(0, n), rep(1, n), W1, W2, W3.0)
    Y.10 <- Q2(rep(1, n), rep(0, n), W1, W2, W3.1)
    Y.11 <- Q2(rep(1, n), rep(1, n), W1, W2, W3.1)
    
    data.frame(W1 = W1, W2 = W2, W3 = W3, A1 = A1, A2 = A2, Y = Y, 
               Y.00 = Y.00, Y.01 = Y.01, Y.10 = Y.10, Y.11 = Y.11)
}

EYd0 <- function(n) {
    dat <- gendata(n)
    Q.00 <- with(dat, Q2(0, 0, W1, W2, W3))
    Q.01 <- with(dat, Q2(0, 1, W1, W2, W3))
    Q.10 <- with(dat, Q2(1, 0, W1, W2, W3))
    Q.11 <- with(dat, Q2(1, 1, W1, W2, W3))

    Q.0 <- with(dat, Q1(0, A2 = as.numeric(Q.01 > Q.00), W1, W2))
    Q.1 <- with(dat, Q1(1, A2 = as.numeric(Q.11 > Q.10), W1, W2))

    first.treat <- as.numeric(Q.1 > Q.0)
    second.treat <- rep(NA, length(first.treat))
    second.treat[first.treat == 0] <- as.numeric(Q.01 > Q.00)[first.treat == 0]
    second.treat[first.treat == 1] <- as.numeric(Q.11 > Q.10)[first.treat == 1]

    Y.d <- vector("numeric", n)
    Y.d[(first.treat == 0) & (second.treat == 0)] <- dat$Y.00[(first.treat == 0) & (second.treat == 0)]
    Y.d[(first.treat == 0) & (second.treat == 1)] <- dat$Y.01[(first.treat == 0) & (second.treat == 1)]
    Y.d[(first.treat == 1) & (second.treat == 0)] <- dat$Y.10[(first.treat == 1) & (second.treat == 0)]
    Y.d[(first.treat == 1) & (second.treat == 1)] <- dat$Y.11[(first.treat == 1) & (second.treat == 1)]

    mean(Y.d)
}
