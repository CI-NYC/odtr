library(tidyverse)
library(glue)

devtools::source_gist("https://gist.github.com/nt-williams/3afb56f503c7f98077722baf9c7eb644")

source("_research/sim/gendata2.R")

res <- map_dfr(1:2, function(x) {
    glue("_research/sim/data/sim_revision{x}_2.zip") |> 
    read_zip_rds() |> 
        bind_rows()
})


truth <- EYd0(1e7)

res <- group_by(res, n) |> 
    summarise(psi = mean(psi), 
              abs_bias = abs(mean(psi - truth)), 
              coverage = mean(map2_lgl(conf.low, conf.high, \(l, h) between(truth, l, h)))) |> 
    mutate(rootn_bias = abs_bias * sqrt(n))

saveRDS(res, "_research/sim/data/simulation_results_revision.rds")

kableExtra::kbl(res,
                format = "latex",
                booktabs = TRUE, 
                digits = 3)
