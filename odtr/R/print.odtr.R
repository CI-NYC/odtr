#' @exportS3Method
print.odtr <- function(obj, ...) {
    print(glue::glue("E[Y_odtr]: {format(obj$psi, digits = 4)} ({format(obj$conf.low, digits = 4)}, {format(obj$conf.high, digits = 4)})"))
    cat("\n")
    print(head(obj$A_opt))
}
