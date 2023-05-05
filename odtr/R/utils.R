g <- glue::glue

at_risk <- function(data, Npsem, tau) {
    if (is.null(Npsem$risk)) {
        return(rep(TRUE, nrow(data)))
    }
    
    if (tau == 1) {
        return(rep(TRUE, nrow(data)))
    }
    
    data[[Npsem$risk[tau - 1]]] == 1 & !is.na(data[[Npsem$risk[tau - 1]]])
}

convert_to_surv <- function(x) {
    data.table::fcase(x == 0, 1,
                      x == 1, 0)
}

bound <- function(x, p = 1e-05) {
    pmax(pmin(x, 1 - p), p)
}

trim <- function(x, trim) {
    pmin(pmax(x, quantile(x, 1 - trim)), quantile(x, trim))
}

