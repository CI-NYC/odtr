Npsem <- R6::R6Class(
    "Npsem",
    public = list(
        W = NULL,
        L = NULL,
        A1 = NULL,
        A2 = NULL,
        Y = NULL,
        initialize = function(W, L, A1, A2, Y) {
            checkmate::assertCharacter(A1)
            
            if (!missing(A2)) {
                checkmate::assertCharacter(A2)
                self$A2 <- A2
            }
            
            if (!missing(W)) {
                checkmate::assertCharacter(W)
                self$W <- W
            }
            
            checkmate::assertCharacter(Y)
            checkmate::assertList(L, types = "character")
            
            self$L <- L
            self$A1 <- A1
            self$Y <- Y
        },
        # Get all parent nodes for a variable
        history = function(var = c("L", "A", "Y"), t) {
            switch(
                match.arg(var),
                L = private$parents_L(t),
                A = private$parents_A(t),
                Y = private$parents_Y()
            )
        },
        # Return the names of all variables
        all_vars = function() {
            c(self$W, unlist(self$L), self$A1, self$A2, self$Y)
        }
    ),
    private = list(
        parents_L = function(t) {
            if (t == 1) {
                return(self$W)
            }
            c(private$parents_A(t - 1), self$A1[t - 1])
        },
        parents_A = function(t) {
            c(private$parents_L(t), unlist(self$L[[t]]))
        },
        parents_Y = function() {
            c(self$W, unlist(self$L), self$A1)
        }
    )
)
