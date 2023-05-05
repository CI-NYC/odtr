#' Optimal Dynamic Treatment Rules for Multiple Timepoints
#'
#' @param data \[\code{data.frame}\]\cr
#'  A \code{data.frame} containing all necessary variables
#'  for the estimation problem. 
#' @param Npsem \[\code{R6(Npsem)}\]\cr
#'  An \code{Npsem} object 
#' @param V \[\code{integer}\]\cr
#'  Number of folds for cross-fitting
#' @param g_learner \[\code{character}\]\cr A vector of \code{mlr3superlearner} algorithms 
#'  for estimation of the propensity score. Default is \code{"glm"}, a main effects GLM.
#' @param Q_learner \[\code{character}\]\cr A vector of \code{mlr3superlearner} algorithms 
#'  for estimation of the outcome regression. Default is \code{"glm"}, a main effects GLM.
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
#' dat <- odtr_sim(1e3)
#' sem <- Npsem$new(paste0("W", 1:4), A = c("A1", "A2"), Y = "Y")
#' optimal <- odtr(dat, sem, 1, "glm", "glm", "binomial")
odtr <- function(data, Npsem, V, g_learner = "glm", Q_learner = "glm", 
                 type = c("binomial", "continuous"), maximize = TRUE) {
    require("mlr3superlearner")
    checkmate::assertDataFrame(data[, Npsem$all_vars()])
    checkmate::assertR6(Npsem, "Npsem")
    checkmate::assertNumber(V, lower = 1, upper = nrow(data) - 1)

    tmp <- data.table::copy(data)
    if (!is.null(Npsem$risk)) {
        for (y in c(Npsem$risk, Npsem$Y)) {
            data.table::set(tmp, j = y, value = convert_to_surv(tmp[[y]]))
        }
    }
    
    folds <- make_folds(tmp, V)
    g0 <- crossFitg0(tmp, Npsem, g_learner, folds)
    vals <- crossFitQ(tmp, g0, Npsem, Q_learner, folds, "binomial", maximize)
    
    inflnce <- eif(data[, Npsem$A], vals$A_opt, g0, vals$m, vals$Q_a)
    
    colnames(vals$A_opt) <- Npsem$A
    se <- sqrt(var(inflnce) / nrow(data))
    list(psi = mean(inflnce), 
         std.error = se,
         conf.low = mean(inflnce) - qnorm(0.975)*se, 
         conf.high = mean(inflnce) + qnorm(0.975)*se, 
         A_opt = data.table::as.data.table(vals$A_opt))
}

eif <- function(A, A_opt, g, mtilde, Q_a, Q_0, Q_1) {
    tau <- ncol(Q_a)
    r <- matrix(1, nrow = nrow(g), ncol = length(1:tau))
    
    t <- 1
    for (s in t:tau) {
        r[, s] <- as.numeric(A[, s] == A_opt[, s]) / g[, s]
        r[is.na(r[, s]), s] <- 0
    }
    
    m <- mtilde[, (t + 1):(tau + 1), drop = FALSE] - Q_a[, t:tau, drop = FALSE]
    rowSums(r * m, na.rm = TRUE) + mtilde[, t]
}
