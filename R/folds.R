make_folds <- function(data, V, strata) {
    folds <- origami::make_folds(data, V = V, strata_ids = strata)
    if (V == 1) {
        folds[[1]]$training_set <- folds[[1]]$validation_set
    }
    folds
}
