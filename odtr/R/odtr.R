#' Optimal Dynamic Treatment Rules for Multiple Timepoints
#'
#' @param data \[\code{data.frame}\]\cr
#'  A \code{data.frame} containing all necessary variables
#'  for the estimation problem. 
#' @param Npsem \[\code{R6(Npsem)}\]\cr
#'  An \code{Npsem} object 
#' @param V \[\code{integer}\]\cr
#'  Number of folds for cross-fitting
#' @param g_learner \[\code{character}\]\cr A vector of \code{SuperLearner} algorithms 
#'  for estimation of the propensity score. Default is \code{"SL.glm"}, a main effects GLM.
#' @param Q_learner \[\code{character}\]\cr A vector of \code{SuperLearner} algorithms 
#'  for estimation of the outcome regression. Default is \code{"SL.glm"}, a main effects GLM.
#' @param type \[\code{character(1)}\]\cr
#'  Outcome variable type (i.e., continuous, binomial).
#' @param maximize \[\code{logical(1)}\]\cr
#'  Should the optimal treatment maximize (or minimize) the outcome?
#'
#' @return A \code{data.frame} of the optimal treatment decisions
#' @export
#'
#' @examples
#' dat <- sampleLuedtke2015(1e3)
#' sem <- Npsem$new(paste0("W", 1:4), A = c("A1", "A2"), Y = "Y")
#' optimal <- odtr(dat, sem, 1, "SL.glm", "Sl.glm", "binomial")
odtr <- function(data, Npsem, V, g_learner = "SL.glm", Q_learner = "SL.glm", 
                 type = c("binomial", "continuous"), maximize = TRUE) {
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
    A_opt <- crossFitQ(tmp, g0, Npsem, Q_learner, folds, "binomial")
    
    colnames(A_opt) <- Npsem$A
    as.data.frame(A_opt)
}
