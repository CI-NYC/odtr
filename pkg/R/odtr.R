odtr <- function(data, Npsem, V, g_learner, Q_learner, type = c("binomial", "continuous")) {
    checkmate::assertDataFrame(data[, Npsem$all_vars()])
    checkmate::assertR6(Npsem, "Npsem")
    checkmate::assertR6(g_learner, "Lrnr_base")
    checkmate::assertR6(Q_learner, "Lrnr_base")
    checkmate::assertNumber(V, lower = 1, upper = nrow(data) - 1)

    tmp <- data.table::copy(data)
    if (!is.null(Npsem$risk)) {
        for (y in c(Npsem$risk, Npsem$Y)) {
            data.table::set(tmp, j = y, value = convert_to_surv(tmp[[y]]))
        }
    }
    
    folds <- make_folds(tmp, V)
    g0 <- crossFitg0(tmp, Npsem, g_learner, folds)
    
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

        Q0 <- crossFitQ0(tmp, y, Npsem$A[t], vars, t, Npsem, Q_learner, folds, type)
        Qv <- crossFitQv(tmp, g0, Q0$Q0, t, Npsem, Q_learner, folds)

        ans[[t]] <- ifelse(Qv[, 1] > 0, 1, 0)
        
        if (t != 1) {
            tmp <- addYd_t(tmp, Npsem, ans[[t]], Q0$fits, folds, t)
        }
    }
    
    setNames(ans, Npsem$A)
}
