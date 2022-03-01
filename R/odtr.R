odtr <- function(data, Npsem, V, g_learner, Q_learner, type = c("binomial", "continuous")) {
    checkmate::assertDataFrame(data[, Npsem$all_vars()], any.missing = FALSE)
    checkmate::assertR6(Npsem, "Npsem")
    checkmate::assertR6(g_learner, "Lrnr_base")
    checkmate::assertR6(Q_learner, "Lrnr_base")
    checkmate::assertNumber(V, lower = 1, upper = nrow(data) - 1)
    
    folds <- make_folds(data, V)
    
    g0 <- crossFitg0(dat, np, g_learner, folds)
    
    ans <- list()
    for (t in length(Npsem$A):1) {
        if (t == length(Npsem$A)) {
            y <- Npsem$Y
            vars <- Npsem$history("Y")
            type <- match.arg(type)
        } else {
            y <- g("tmp_Yd_{t+1}")
            vars <- Npsem$history("L", t + 1)
            type <- "continuous"
        }
        
        Q0 <- crossFitQ0(data, y, np$A[t], vars, Q_learner, folds, type)
        Qv <- crossFitQv(data, g0, Q0$Q0, t, Npsem, Q_learner, folds)
        ans[[t]] <- ifelse(Qv[, 1] > 0, 1, 0)
        
        if (t != 1) {
            data <- addYd_t(data, Npsem, ans[[t]], Q0$fits, folds, t)
        }
    }
    
    setNames(ans, Npsem$A)
}
