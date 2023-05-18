gendata <- function(n) {
    W1 <- rbinom(n, 1, 0.5)
    A1 <- rbinom(n, 1, 0.5)
    # A1 <- rbinom(n, 1, plogis(-log(2) + log(1.7)*W1 - log(2)*W2))
    L1 <- rbinom(n, 1, plogis(log(5) - log(4)*A1 + log(1.75)*W1))
    A2 <- rbinom(n, 1, 0.5)
    # A2 <- rbinom(n, 1, plogis(-log(6) + log(2)*A1 + log(1.3)*W1*L1 - log(1.4)*W2 + log(3)*L2))
    
    Qbar0 <- function(W1, A1, L1, A2) {
        plogis(0.65 - 0.3*A1 + 0.2*A2 + 0.2*W1 + 0.6*L1 - 0.1*A2*W1 - 0.05*A1*W1 + 0.75*A1*L1 - 0.45*A2*L1)
    }
    
    u <- runif(n)
    Y <- as.numeric(u < Qbar0(W1, A1, L1, A2))
    
    Y.00 <- Qbar0(W1, 0, L1, 0)
    Y.01 <- Qbar0(W1, 0, L1, 1)
    Y.10 <- Qbar0(W1, 1, L1, 0)
    Y.11 <- Qbar0(W1, 1, L1, 1)
    
    probs <- data.frame(Y.00, Y.01, Y.10, Y.11)
    
    data.frame(W1 = as.numeric(W1), 
               A1 = as.numeric(A1), 
               L1 = as.numeric(L1), 
               A2 = as.numeric(A2), 
               Y = as.numeric(Y), 
               Y.00 = Y.00, 
               Y.01 = Y.01, 
               Y.10 = Y.10, 
               Y.11 = Y.11,
               d = c("00", "01", "10", "11")[apply(probs, 1, which.max)],
               Ymax = apply(probs, 1, max))
}
