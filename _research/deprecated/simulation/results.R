library(tidyverse)
library(progressr)

read_zip <- function(tar) {
    files <- unzip(tar, list = TRUE)$Name
    p <- progressor(along = 1:length(files))
    map(files, function(file) {
        p()
        con <- unz(tar, file)
        read.csv(con)
    })
}

handlers(global = TRUE)

res <- read_zip("data/simulation/sim5.zip") |> 
    bind_rows()

hist(res$psi)
mean(res$psi)

map2_lgl(res$conf_low, res$conf_high, \(x, y) between(0.485, x, y)) |> 
    mean()

map2_lgl(0.485 - res$std_error*1.96, 
         0.485 + res$std_error*1.96, 
         \(x, y) between(0.485, x, y)) |> mean()


