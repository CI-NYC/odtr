suppressPackageStartupMessages(library(tidyverse))

source("_research/sim/gendata2.R")
devtools::source_gist("https://gist.github.com/nt-williams/3afb56f503c7f98077722baf9c7eb644")

res <- read_zip_rds("_research/sim/data/sim.zip") |> 
    bind_rows()

truth <- EYd0(1e7)

res <- group_by(res, n) |> 
    summarise(psi = mean(psi), 
              abs_bias = abs(mean(psi - truth)), 
              coverage = mean(map2_lgl(conf.low, conf.high, \(l, h) between(truth, l, h)))) |> 
    mutate(rootn_bias = abs_bias * sqrt(n))

saveRDS(res, "_research/sim/data/simulation_results.rds")
