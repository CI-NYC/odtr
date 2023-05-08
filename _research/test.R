dat <- odtr_sim(1e4)
sem <- Npsem$new(paste0("W", 1:3), A = c("A1", "A2"), Y = "Y")
opt <- odtr(dat, sem, 1, c("glm", "nnet"), c("glm", "nnet"), "binomial")

d <- dat
tr <- strsplit(dat$d, "")
d$A1 <- as.numeric(purrr::map_chr(tr, 1))
d$A2 <- as.numeric(purrr::map_chr(tr, 2))

table(d[, c("A1", "A2")])
table(opt$A_opt)
mean(dat$Ymax)
