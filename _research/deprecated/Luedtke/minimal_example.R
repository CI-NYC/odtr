library(SuperLearner)

source("scripts/Luedtke/funs.R")
devtools::load_all("pkg")

V.list <- list(c('W1', 'W2'), c('W1', 'W2', 'W3', 'W4'))
L.list <- list(c('W1', 'W2'), c('W1', 'W2', 'W3', 'W4'))

QbV.SL.library <- c('SL.glm','SL.rpart')
cl.SL.library <- c('SL.glm.interaction','SL.mean')

seed <- 356
set.seed(seed)

ObsData <- sampleLuedtke2015(1e4)
L <- subset(ObsData, select = paste('W', 1:4, sep = ''))
A <- subset(ObsData, select = paste('A', 1:2, sep = ''))
Y <- ObsData$Y

Ey_d_luedtke <- EYd.mtp.tmle(L, A, Y, QbV.SL.library, NULL, 
                             Qbar0.list = list(Qbar0.mtp.ipcw1, Qbar0.mtp.ipcw2), 
                             V.list = V.list, L.list = L.list, Q.SL.library = QbV.SL.library, mtp.jointSL.list = NULL, 
                             g = 0.5, alpha.risk.est = 'emp.risk', bd = NULL, num.folds = 10, cv.ests = TRUE)

table(
    as.character(Ey_d_luedtke$d$A2), 
    unlist(lapply(strsplit(ObsData$d, split = ""), \(x) x[2]))
)

mean(as.character(Ey_d_luedtke$d$A2) == unlist(lapply(strsplit(ObsData$d, split = ""), \(x) x[2])))

table(
    paste0(Ey_d_luedtke$d$A1, Ey_d_luedtke$d$A2), 
    ObsData$d
)

mean(paste0(Ey_d_luedtke$d$A1, Ey_d_luedtke$d$A2) == ObsData$d)
