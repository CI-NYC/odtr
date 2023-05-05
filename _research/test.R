data <- sampleLuedtke2015(5e4)
sem <- Npsem$new(paste0("W", 1:4), A = c("A1", "A2"), Y = "Y")
optimal <- odtr(data, sem, 1, "SL.glm", c("SL.glm", "SL.lightgbm"), "binomial")

shifted <- data
for (a in c("A1", "A2")) {
    shifted[[a]] <- optimal[, a]
}

lmtp::lmtp_sdr(
    data, c("A1", "A2"), "Y", paste0("W", 1:4),
    shifted = shifted, folds = 1, learners_outcome = c("SL.glm", "SL.lightgbm")
)

d <- data
tr <- strsplit(data$d, "")
d$A1 <- as.numeric(purrr::map_chr(tr, 1))
d$A2 <- as.numeric(purrr::map_chr(tr, 2))

lmtp::lmtp_sdr(
    data, c("A1", "A2"), "Y", paste0("W", 1:4),
    shifted = d, folds = 1, learners_outcome = c("SL.glm", "SL.lightgbm")
)
