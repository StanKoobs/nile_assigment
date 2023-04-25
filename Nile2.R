##############################################################################
# Assignment 3 Advanced Econometrics III
##############################################################################
# Packages
library(ggplot2)
library(cowplot)
library(tseries)
library(resample)

OwnggTheme <- function() {
  theme(panel.background = element_rect(fill = "grey99"),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "gray94", fill = NA,
                                    size = 1.5),
        panel.grid.major = element_line(colour = "grey85",
                                        linetype = "dashed"),
        panel.grid.minor = element_line(colour = "grey85",
                                        linetype = "dashed"),
        axis.text.x = element_text(size= 12),
        axis.text.y = element_text(size= 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.ticks = element_line(size = 1.5),
        axis.ticks.length = unit(1.5, "mm"),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14)
  ) 
}

# Loading the data 
library(readxl)
df <- read_excel("C:/Users/72692svm/Downloads/Nile.xlsx")

y <- df$Nile
year <- df$...1
n <- length(y)

# The code below is copied from the previous assignment

##############################################################################
# Maximum likelihood estimation

logL <- function(sigma) {
  a1 <- 0
  P1 <- 10^7
  
  sigmaeps <- sigma[1]
  sigmaeta <- sigma[2]
  
  avec <- c(a1)
  Pvec <- c(P1)
  Ptvec <- c()
  Fvec <- c()
  vvec <- c()
  atvec <- c()
  Ktvec <- c()
  
  for (t in 1:length(y)) {
    # Update
    vvec[t] <- y[t] - avec[t]
    Fvec[t] <- Pvec[t] + sigmaeps
    Ktvec[t] <- Pvec[t] / Fvec[t]
    Ptvec[t] <- Pvec[t] * (1 - Ktvec[t])
    atvec[t] <- avec[t] + Ktvec[t] * vvec[t]
    
    # Predict
    avec[t + 1] <- avec[t] + Ktvec[t] * vvec[t]
    Pvec[t + 1] <- Pvec[t] * (1 - Ktvec[t]) + sigmaeta
  }
  
  logL <- -0.5 * sum(log(Fvec[-1])) - 0.5 * sum(vvec[-1]^2 / Fvec[-1])
  return(logL)
}

opt <- optim(c(10, 10), logL, method = "L-BFGS-B", lower = c(0.1, 0.1),
             control = list(fnscale = -1))
opt
# Same estimates as found in DK 

##############################################################################
# Kalman filter

a1 <- 0
P1 <- 10^7
sigmaeps <- 15099
sigmaeta <- 1469.1

avec <- c(a1)
Pvec <- c(P1)
Ptvec <- c()
Fvec <- c()
vvec <- c()
atvec <- c()
Ktvec <- c()

for (t in 1:length(y)) {
  # Update
  vvec[t] <- y[t] - avec[t]
  Fvec[t] <- Pvec[t] + sigmaeps
  Ktvec[t] <- Pvec[t] / Fvec[t]
  Ptvec[t] <- Pvec[t] * (1 - Ktvec[t])
  atvec[t] <- avec[t] + Ktvec[t] * vvec[t]
  
  # Predict
  avec[t + 1] <- avec[t] + Ktvec[t] * vvec[t]
  Pvec[t + 1] <- Pvec[t] * (1 - Ktvec[t]) + sigmaeta
}

##############################################################################
# State smoothing

rvec <- rep(0, n)
alphahat <- rep(0, n)
Lvec <- sigmaeps / Fvec
Nvec <- rep(0, n)
Vvec <- rep(0, n)

# Backward recursion
for (t in n:2) {
  rvec[t - 1] <- 1 / Fvec[t] * vvec[t] + Lvec[t] * rvec[t]
  alphahat[t] <- avec[t] + Pvec[t] * rvec[t - 1]
  
  Nvec[t - 1] <- 1 / Fvec[t] + Lvec[t]^2 * Nvec[t]
  Vvec[t] <- Pvec[t] - Pvec[t]^2 * Nvec[t - 1]
}
r0 <- 1 / Fvec[1] * vvec[1] + Lvec[1] * rvec[1]
alphahat[1] <- avec[1] + Pvec[1] * r0

N0 <- 1 / Fvec[1] + Lvec[1]^2 * Nvec[1]
Vvec[1] <- Pvec[1] - Pvec[1]^2 * N0

##############################################################################
# Question (c)
# Disturbance smoothing

# We can calculate smoothed observation disturbance directly
# However, we use the expression in terms of r_t and N_t

rvec <- rep(0, n)
Lvec <- sigmaeps / Fvec
Nvec <- rep(0, n)

# Backward recursion for only r and N
for (t in n:2) {
  rvec[t - 1] <- 1 / Fvec[t] * vvec[t] + Lvec[t] * rvec[t]
  Nvec[t - 1] <- 1 / Fvec[t] + Lvec[t]^2 * Nvec[t]
}
r0 <- 1 / Fvec[1] * vvec[1] + Lvec[1] * rvec[1]
N0 <- 1 / Fvec[1] + Lvec[1]^2 * Nvec[1]

# Calculate smoothed (observation / state) disturbances and corresponding
# variances
epshat <- sigmaeps * (1 / Fvec * vvec - Ktvec * rvec)
vareps <- sigmaeps - sigmaeps^2 * (1 / Fvec + Ktvec^2 * Nvec)
etahat <- sigmaeta * rvec
vareta <- sigmaeta - sigmaeta^2 * Nvec


##############################################################################
# Question (d)
# Simulation smoothing (from last part)

alphaplus <- rnorm(1, 0, 100)
yplus <- c()

epsplus <- rnorm(n, 0, sd = sqrt(sigmaeps))
etaplus <- rnorm(n, 0, sd = sqrt(sigmaeta))

for (t in 1:n) {
  alphaplus[t + 1] <- alphaplus[t] + etaplus[t]
  yplus[t] <- alphaplus[t] + epsplus[t]
}

# Kalman filter for the new data
avecp <- c(a1)
Pvecp <- c(P1)
Ptvecp <- c()
Fvecp <- c()
vvecp <- c()
atvecp <- c()
Ktvecp <- c()

for (t in 1:n) {
  # Update
  vvecp[t] <- yplus[t] - avecp[t]
  Fvecp[t] <- Pvecp[t] + sigmaeps
  Ktvecp[t] <- Pvecp[t] / Fvecp[t]
  Ptvecp[t] <- Pvecp[t] * (1 - Ktvecp[t])
  atvecp[t] <- avecp[t] + Ktvecp[t] * vvecp[t]
  
  # Predict
  avecp[t + 1] <- avecp[t] + Ktvecp[t] * vvecp[t]
  Pvecp[t + 1] <- Pvecp[t] * (1 - Ktvecp[t]) + sigmaeta
}

rvecp <- rep(0, n)
Lvecp <- sigmaeps / Fvecp
Nvecp <- rep(0, n)

# Backward recursion for only r and N
for (t in n:2) {
  rvecp[t - 1] <- 1 / Fvecp[t] * vvecp[t] + Lvecp[t] * rvecp[t]
  Nvecp[t - 1] <- 1 / Fvecp[t] + Lvecp[t]^2 * Nvecp[t]
}
r0 <- 1 / Fvecp[1] * vvecp[1] + Lvecp[1] * rvecp[1]
N0 <- 1 / Fvecp[1] + Lvecp[1]^2 * Nvecp[1]

# Calculate smoothed (observation / state) disturbances and corresponding
# variances
epshatp <- sigmaeps * (1 / Fvecp * vvecp - Ktvecp * rvecp)

epstilde <- epsplus - epshatp + epshat

alphatilde <- y - epstilde

ggplot() +
  geom_point(aes(x = year, y = alphaplus[1:n])) +
  geom_line(aes(x = year, y = alphahat))

ggplot() +
  geom_point(aes(x = year, y = alphatilde)) +
  geom_line(aes(x = year, y = alphahat))

##############################################################################
# We now also calculate them based on the approach from the lecture which is
# almost similar but does not use the epsilons. 

alphatildefunc <- function() {
  alphaplus <- y[1]
  
  epsplus <- rnorm(n, 0, sd = sqrt(sigmaeps))
  etaplus <- rnorm(n, 0, sd = sqrt(sigmaeta))
  
  for (t in 1:n) {
    alphaplus[t + 1] <- alphaplus[t] + etaplus[t]
    yplus[t] <- alphaplus[t] + epsplus[t]
  }
  
  # Kalman filter for the new data
  avecp <- c(a1)
  Pvecp <- c(P1)
  Ptvecp <- c()
  Fvecp <- c()
  vvecp <- c()
  atvecp <- c()
  Ktvecp <- c()
  
  for (t in 1:n) {
    # Update
    vvecp[t] <- yplus[t] - avecp[t]
    Fvecp[t] <- Pvecp[t] + sigmaeps
    Ktvecp[t] <- Pvecp[t] / Fvecp[t]
    Ptvecp[t] <- Pvecp[t] * (1 - Ktvecp[t])
    atvecp[t] <- avecp[t] + Ktvecp[t] * vvecp[t]
    
    # Predict
    avecp[t + 1] <- avecp[t] + Ktvecp[t] * vvecp[t]
    Pvecp[t + 1] <- Pvecp[t] * (1 - Ktvecp[t]) + sigmaeta
  }
  
  rvecp <- rep(0, n)
  alphahatp <- rep(0, n)
  Lvecp <- sigmaeps / Fvec
  Nvecp <- rep(0, n)
  Vvecp <- rep(0, n)
  
  # Backward recursion
  for (t in n:2) {
    rvecp[t - 1] <- 1 / Fvecp[t] * vvecp[t] + Lvecp[t] * rvecp[t]
    alphahatp[t] <- avecp[t] + Pvecp[t] * rvecp[t - 1]
    
    Nvecp[t - 1] <- 1 / Fvecp[t] + Lvecp[t]^2 * Nvecp[t]
    Vvecp[t] <- Pvecp[t] - Pvecp[t]^2 * Nvecp[t - 1]
  }
  r0p <- 1 / Fvecp[1] * vvecp[1] + Lvecp[1] * rvecp[1]
  alphahatp[1] <- avecp[1] + Pvecp[1] * r0p
  
  N0p <- 1 / Fvecp[1] + Lvecp[1]^2 * Nvecp[1]
  Vvecp[1] <- Pvecp[1] - Pvecp[1]^2 * N0p
 
  return(alphahat - alphahatp + alphaplus[-101])
}

# Number of iterations
M <- 10000

# Matrix with h_t estimates
hmat <- matrix(NA, nrow = M, ncol = n)

set.seed(1)
for (i in 1:M) {
  alphat <- alphatildefunc()
  hmat[i, ] <- sqrt(alphat / (2 * pi))
}

hmean <- colMeans(hmat)
htilde_naive <- sqrt(alphahat / (2 * pi))

ggplot() +
  geom_line(aes(x = year, y = htilde_naive), colour = "#F8766D", size = 1) +
  geom_line(aes(x = year, y = hmean), colour = "#529EFF", size = 1) +
  scale_x_continuous(n.breaks = 6) +
  OwnggTheme() +
  xlab("Year") +
  ylab("")

tildevars <- colVars(hmat)

ggplot() +
  geom_line(aes(x = year, y = tildevars), colour = "#529EFF", size = 1) +
  scale_x_continuous(n.breaks = 6) +
  OwnggTheme() +
  xlab("Year") +
  ylab("")

