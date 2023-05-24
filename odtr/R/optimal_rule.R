#' Optimal Dynamic Treatment Rules with Time-Varying Data
#'
#' @param data \[\code{data.frame}\]\cr
#'  A \code{data.frame} containing all necessary variables
#'  for the estimation problem. 
#' @param Vars \[\code{R6(Vars)}\]\cr
#'  An \code{Vars} object 
#' @param folds \[\code{integer}\]\cr
#'  Number of folds for cross-fitting
#' @param g_learner \[\code{character}\]\cr 
#'  A \code{SuperLearner} library
#' @param Q_learner \[\code{character}\]\cr 
#'  A \code{SuperLearner} library
#' @param type \[\code{character(1)}\]\cr
#'  Outcome variable type (i.e., continuous, binomial).
#' @param maximize \[\code{logical(1)}\]\cr
#'  Should the optimal treatment maximize (or minimize) the outcome?
#'
#' @return A \code{list} containing optimal treatment assignment and the estimated mean value of
#'  the outcome under the optimal treatment rule. 
#' @export
#'
#' @examples
optimal_rule <- function(data, Vars, folds, g_learner = "glm", Q_learner = "glm", 
                         type = c("binomial", "continuous"), maximize = TRUE) {
    checkmate::assertDataFrame(data[, Vars$all_vars()])
    checkmate::assertR6(Vars, "Vars")
    checkmate::assertNumber(folds, lower = 1, upper = nrow(data) - 1)

    tmp <- data.table::copy(data)
    if (!is.null(Vars$risk)) {
        for (y in c(Vars$risk, Vars$Y)) {
            data.table::set(tmp, j = y, value = convert_to_surv(tmp[[y]]))
        }
    }
    
    folds <- make_folds(tmp, folds)
    g0 <- crossFitg0(tmp, Vars, g_learner, folds)
    vals <- crossFitQ(tmp, g0, Vars, Q_learner, folds, match.arg(type), maximize)
    
    colnames(vals$A_opt) <- Vars$A
    inflnce <- eif(data[, Vars$A, drop = F], vals$A_opt, g0, vals$m, vals$Q_a)
    se <- sqrt(var(inflnce) / nrow(data))
    
    returns <- list(psi = mean(inflnce), 
                    std.error = se,
                    conf.low = mean(inflnce) - qnorm(0.975)*se, 
                    conf.high = mean(inflnce) + qnorm(0.975)*se, 
                    A_opt = as.data.frame(vals$A_opt), 
                    eif = inflnce,
                    m = vals$m, 
                    g = g0)
    class(returns) <- "odtr"
    returns
}

eif <- function(A, A_opt, g, mtilde, Q_a, Q_0, Q_1) {
    tau <- ncol(Q_a)
    r <- matrix(1, nrow = nrow(g), ncol = length(1:tau))
    
    t <- 1
    for (s in t:tau) {
        r[, s] <- as.numeric(A[, s] == A_opt[, s]) / g[, s]
        r[is.na(r[, s]), s] <- 0
    }
    
    r <- apply(r, 2, function(x) pmin(x, quantile(x, 0.99, na.rm = T)))
    m <- mtilde[, (t + 1):(tau + 1), drop = FALSE] - Q_a[, t:tau, drop = FALSE]
    rowSums(r * m, na.rm = TRUE) + mtilde[, t]
}
