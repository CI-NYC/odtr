#' Optimal Dynamic Treatment Rules with Time-Varying Data
#'
#' @param data \[\code{data.frame}\]\cr
#'  A \code{data.frame} containing all necessary variables
#'  for the estimation problem. 
#' @param trt 
#' @param outcome 
#' @param baseline 
#' @param time_varying
#' @param learners_trt \[\code{character}\]\cr 
#'  A \code{mlr3superlearner} library
#' @param learners_outcome \[\code{character}\]\cr 
#'  A \code{mlr3superlearner} library
#' @param learners_rule \[\code{character}\]\cr 
#'  A \code{mlr3superlearner} library
#' @param folds \[\code{integer}\]\cr
#'  Number of folds for cross-fitting
#' @param outcome_type \[\code{character(1)}\]\cr
#'  Outcome variable type (i.e., continuous, binomial).
#' @param maximize \[\code{logical(1)}\]\cr
#'  Should the optimal treatment maximize (or minimize) the outcome?
#'
#' @return A \code{list} containing optimal treatment assignment and the estimated mean value of the outcome under the optimal treatment rule. 
#' @export
#'
#' @examples
odtr <- function(data, trt, outcome, baseline, time_varying, 
                 learners_trt = "glm", learners_outcome = "glm", learners_rule = "glm",
                 folds, outcome_type = c("binomial", "continuous"), maximize = TRUE) {
    task <- odtr_Vars$new(baseline, time_varying, trt, outcome)
    
    checkmate::assertDataFrame(data[, task$all_vars()])
    checkmate::assertNumber(folds, lower = 1, upper = nrow(data) - 1)

    tmp <- data.table::copy(data)
    if (!is.null(task$risk)) {
        for (y in c(task$risk, task$Y)) {
            data.table::set(tmp, j = y, value = convert_to_surv(tmp[[y]]))
        }
    }
    
    folds <- make_folds(tmp, folds, tmp[[task$Y]])
    g0 <- crossFitg0(tmp, task, learners_trt, folds)
    vals <- crossFitQ(tmp, g0, task, learners_outcome, learners_rule, 
                      folds, match.arg(outcome_type), maximize)

    colnames(vals$A_opt) <- task$A
    inflnce <- eif(data[, task$A, drop = F], vals$A_opt, g0, vals$m, vals$Q_a)
    se <- sqrt(var(inflnce) / nrow(data))

    returns <- list(psi = mean(inflnce), 
                    std.error = se,
                    conf.low = mean(inflnce) - qnorm(0.975)*se, 
                    conf.high = mean(inflnce) + qnorm(0.975)*se, 
                    A_opt = as.data.frame(vals$A_opt), 
                    eif = inflnce,
                    m = vals$m, 
                    g = g0, 
                    decision_fits = vals$fits)
    class(returns) <- "odtr"
    returns
}

eif <- function(A, A_opt, g, mtilde, Q_a, Q_0, Q_1) {
    tau <- ncol(Q_a)
    r <- matrix(1, nrow = nrow(g), ncol = length(1:tau))
    
    for (s in 1:tau) {
        r[, s] <- as.numeric(A[, s] == A_opt[, s]) / g[, s]
        r[is.na(r[, s]), s] <- 0
    }
    
    r <- apply(r, 2, function(x) pmin(x, quantile(x, 0.99, na.rm = T)))
    m <- mtilde[, 2:(tau + 1), drop = FALSE] - Q_a[, 1:tau, drop = FALSE]
    rowSums(compute_weights(r, 1, tau) * m, na.rm = TRUE) + mtilde[, 1]
}
