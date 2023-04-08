##############################################################################
# Assignment 2 Advanced Econometrics III
##############################################################################
# Packages
library(ggplot2)
library(cowplot)
library(tseries)

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

# Synthetic data
sigmaeta <- 1
sigmaeps <- 1

n = 100
alpha1 <- 0
alpha <- c(alpha1)
y <- c()

for (t in 1:n) {
  alpha[t + 1] <- alpha[t] + rnorm(1, 0, sigmaeta)
  y[t] <- alpha[t] + rnorm(1, 0, sigmaeps)
}

# Loading the data 
library(readxl)
df <- read_excel("nile.xlsx")

y <- df$Nile
year <- df$...1

##############################################################################
# Question (a)

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

# We remove a_1 and plot for t = 2,...,n
p1 <- ggplot() +
  geom_point(aes(x = year, y = y), colour = "#529EFF", size = 1.5) +
  geom_line(aes(x = year[2:n], y = avec[2:n]), colour = "#F8766D", size = 1) +
  geom_line(aes(x = year[2:n], y = avec[2:n] + qnorm(0.95) * sqrt(Pvec[2:n])), 
            colour = "#F8766D", alpha = 0.7) +
  geom_line(aes(x = year[2:n], y = avec[2:n] - qnorm(0.95) * sqrt(Pvec[2:n])), 
            colour = "#F8766D", alpha = 0.7) +
  ylim(c(min(y), max(y))) + 
  scale_x_continuous(n.breaks = 6) +
  OwnggTheme() +
  xlab("Year") +
  ylab("")


p2 <- ggplot() + 
  geom_line(aes(x = year[2:n], y = Pvec[2:n])) +
  OwnggTheme() +
  scale_x_continuous(n.breaks = 6) +
  xlab("Year") +
  ylab("")

p3 <- ggplot() +
  geom_line(aes(x = year[2:n], y = vvec[2:n])) +
  geom_hline(yintercept = 0) +
  scale_x_continuous(n.breaks = 6) +
  OwnggTheme() +
  xlab("Year") +
  ylab("")

p4 <- ggplot() +
  geom_line(aes(x = year[2:n], y = Fvec[2:n])) +
  scale_x_continuous(n.breaks = 6) +
  OwnggTheme() +
  xlab("Year") +
  ylab("")

plot_grid(p1, p2, p3, p4)

##############################################################################
# Question (b)
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

p1 <- ggplot() +
  geom_point(aes(x = year, y = y), colour = "#529EFF", size = 1.5) +
  geom_line(aes(x = year, y = alphahat), colour = "#F8766D", size = 1) +
  geom_line(aes(x = year, y = alphahat + qnorm(0.95) * sqrt(Vvec)), 
            colour = "#F8766D", alpha = 0.7) +
  geom_line(aes(x = year, y = alphahat - qnorm(0.95) * sqrt(Vvec)), 
            colour = "#F8766D", alpha = 0.7) +
  ylim(c(min(y), max(y))) + 
  scale_x_continuous(n.breaks = 6) +
  OwnggTheme() +
  xlab("Year") +
  ylab("")

p2 <- ggplot() + 
  geom_line(aes(x = year, y = Vvec)) +
  OwnggTheme() +
  scale_x_continuous(n.breaks = 6) +
  xlab("Year") +
  ylab("")

p3 <- ggplot() + 
  geom_line(aes(x = year, y = c(r0, rvec[-n]))) +
  OwnggTheme() +
  scale_x_continuous(n.breaks = 6) +
  geom_hline(yintercept = 0) +
  xlab("Year") +
  ylab("")

p4 <- ggplot() + 
  geom_line(aes(x = year[-n], y = Nvec[-n])) +
  OwnggTheme() +
  scale_x_continuous(n.breaks = 6) +
  xlab("Year") +
  ylab("") 

plot_grid(p1, p2, p3, p4)

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

p1 <- ggplot() +
  geom_line(aes(x = year, y = epshat)) +
  OwnggTheme() +
  scale_x_continuous(n.breaks = 6) +
  geom_hline(yintercept = 0) +
  xlab("Year") +
  ylab("")

p2 <- ggplot() +
  geom_line(aes(x = year, y = sqrt(vareps))) +
  OwnggTheme() +
  scale_x_continuous(n.breaks = 6) +
  xlab("Year") +
  ylab("")

p3 <- ggplot() +
  geom_line(aes(x = year, y = etahat)) +
  OwnggTheme() +
  scale_x_continuous(n.breaks = 6) +
  geom_hline(yintercept = 0) +
  xlab("Year") +
  ylab("")

p4 <- ggplot() +
  geom_line(aes(x = year, y = sqrt(vareta))) +
  OwnggTheme() +
  scale_x_continuous(n.breaks = 6) +
  xlab("Year") +
  ylab("")

plot_grid(p1, p2, p3, p4)

##############################################################################
# Question (d)
# Simulation smoothing (not in the assignment itself but done for fun)

alphaplus <- c(1120)
yplus <- c()
set.seed(100)

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
# Question (e)
# Missing data

ymis <- y
ymis[c(21:40, 61:80)] <- NA

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

for (t in 1:n) {
  # Update
  if (!is.na(ymis[t])) {
    vvec[t] <- ymis[t] - avec[t]
    Fvec[t] <- Pvec[t] + sigmaeps
    Ktvec[t] <- Pvec[t] / Fvec[t]
    Ptvec[t] <- Pvec[t] * (1 - Ktvec[t])
    atvec[t] <- avec[t] + Ktvec[t] * vvec[t]
    
    # Predict
    avec[t + 1] <- avec[t] + Ktvec[t] * vvec[t]
    Pvec[t + 1] <- Pvec[t] * (1 - Ktvec[t]) + sigmaeta
  } else {
    avec[t + 1] <- avec[t]
    Pvec[t + 1] <- Pvec[t] + sigmaeta
  }
}

# Smoothing
rvec <- rep(0, n)
alphahat <- rep(0, n)
Lvec <- sigmaeps / Fvec
Nvec <- rep(0, n)
Vvec <- rep(0, n)

timemis <- c(21:40, 61:80)

# Backward recursion
for (t in n:2) {
  if (!(t %in% timemis)) {
    rvec[t - 1] <- 1 / Fvec[t] * vvec[t] + Lvec[t] * rvec[t]
    alphahat[t] <- avec[t] + Pvec[t] * rvec[t - 1]
    
    Nvec[t - 1] <- 1 / Fvec[t] + Lvec[t]^2 * Nvec[t]
    Vvec[t] <- Pvec[t] - Pvec[t]^2 * Nvec[t - 1]
  } else {
    rvec[t - 1] <- rvec[t]
    alphahat[t] <- avec[t] + Pvec[t] * rvec[t - 1]
    
    Nvec[t - 1] <- Nvec[t] 
    Vvec[t] <- Pvec[t] - Pvec[t]^2 * Nvec[t - 1]
  }
  
}
r0 <- 1 / Fvec[1] * vvec[1] + Lvec[1] * rvec[1]
alphahat[1] <- avec[1] + Pvec[1] * r0

N0 <- 1 / Fvec[1] + Lvec[1]^2 * Nvec[1]
Vvec[1] <- Pvec[1] - Pvec[1]^2 * N0

# Here we need the filtered estimates instead of the predicted estimates
# So we use avec[-1]
p1 <- ggplot() +
  geom_point(aes(x = year, y = ymis), colour = "#529EFF", size = 1.5) +
  geom_line(aes(x = year, y = ymis), colour = "#529EFF", size = 0.5, 
            alpha = .7) +
  geom_line(aes(x = year, y = avec[-1]), colour = "#F8766D", size = 1) +
  ylim(c(min(y), max(y))) + 
  scale_x_continuous(n.breaks = 6) +
  OwnggTheme() +
  xlab("Year") +
  ylab("")

p2 <- ggplot() +
  geom_line(aes(x = year[2:n], y = Pvec[2:n])) +
  OwnggTheme() +
  scale_x_continuous(n.breaks = 6) +
  xlab("Year") +
  ylab("")

p3 <- ggplot() +
  geom_point(aes(x = year, y = ymis), colour = "#529EFF", size = 1.5) +
  geom_line(aes(x = year, y = ymis), colour = "#529EFF", size = 0.5, 
            alpha = .7) +
  geom_line(aes(x = year, y = alphahat), colour = "#F8766D", size = 1) +
  ylim(c(min(y), max(y))) + 
  scale_x_continuous(n.breaks = 6) +
  OwnggTheme() +
  xlab("Year") +
  ylab("")

p4 <- ggplot() +
  geom_line(aes(x = year, y = Vvec)) +
  OwnggTheme() +
  scale_x_continuous(n.breaks = 6) +
  xlab("Year") +
  ylab("")

plot_grid(p1, p2, p3, p4)

##############################################################################
# Question (f)
# Forecasting new data

nnew <- 130
yearnew <- c(year, 1971:2000)

ymis <- c(y, rep(NA, 30))

avec <- c(a1)
Pvec <- c(P1)
Ptvec <- c()
Fvec <- c()
vvec <- c()
atvec <- c()
Ktvec <- c()

for (t in 1:nnew) {
  # Update
  if (!is.na(y[t])) {
    vvec[t] <- ymis[t] - avec[t]
    Fvec[t] <- Pvec[t] + sigmaeps
    Ktvec[t] <- Pvec[t] / Fvec[t]
    Ptvec[t] <- Pvec[t] * (1 - Ktvec[t])
    atvec[t] <- avec[t] + Ktvec[t] * vvec[t]
    
    # Predict
    avec[t + 1] <- avec[t] + Ktvec[t] * vvec[t]
    Pvec[t + 1] <- Pvec[t] * (1 - Ktvec[t]) + sigmaeta
  } else {
    Fvec[t] <- Pvec[t] + sigmaeps
    
    avec[t + 1] <- avec[t]
    Pvec[t + 1] <- Pvec[t] + sigmaeta
  }
}

p1 <- ggplot() +
  geom_point(aes(x = yearnew, y = ymis), colour = "#529EFF", size = 1.5) +
  geom_line(aes(x = yearnew, y = avec[-1]), colour = "#F8766D", size = 1) +
  geom_line(aes(x = 1971:2000, 
                y = avec[101:130] + qnorm(0.75) * sqrt(Pvec[101:130])), 
            colour = "#F8766D", alpha = .7) +
  geom_line(aes(x = 1971:2000, 
                y = avec[101:130] - qnorm(0.75) * sqrt(Pvec[101:130])), 
            colour = "#F8766D", alpha = .7) +
  ylim(c(min(y), max(y))) + 
  scale_x_continuous(n.breaks = 3) +
  OwnggTheme() +
  xlab("Year") +
  ylab("")

p2 <- ggplot() +
  geom_line(aes(x = yearnew[2:nnew], y = Pvec[2:nnew])) +
  OwnggTheme() +
  scale_x_continuous(n.breaks = 3) +
  xlab("Year") +
  ylab("")

p3 <- ggplot() +
  geom_line(aes(x = yearnew, y = avec[-1])) +
  ylim(c(min(avec[-1]), max(avec))) + 
  scale_x_continuous(n.breaks = 3) +
  OwnggTheme() +
  xlab("Year") +
  ylab("")

p4 <- ggplot() +
  geom_line(aes(x = yearnew[2:nnew], y = Fvec[2:nnew])) +
  OwnggTheme() +
  scale_x_continuous(n.breaks = 3) +
  xlab("Year") +
  ylab("")

plot_grid(p1, p2, p3, p4)

##############################################################################
# Question (g)
# Diagnostic checking standardised prediction errors 

# Using the code from part (a)
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

# Standardized residuals
evec <- vvec[-1] / sqrt(Fvec[-1])

p1 <- ggplot() +
  geom_line(aes(x = year[-1], y = evec)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  scale_x_continuous(n.breaks = 6) +
  OwnggTheme() +
  xlab("Year") +
  ylab("")

p2 <- ggplot(mapping = aes(x = evec)) +
  geom_histogram(aes(y = ..density..), fill = "white", colour = "black",
                 binwidth = 0.4) +
  geom_density() +
  OwnggTheme() +
  scale_x_continuous(n.breaks = 6) +
  xlab("") +
  ylab("")

p3 <- ggplot(mapping = (aes(sample = evec))) +
  stat_qq(geom = "line", colour = "grey45") +
  stat_qq_line() +
  OwnggTheme() +
  scale_x_continuous(n.breaks = 6) +
  xlab("") +
  ylab("") +
  geom_hline(yintercept = 0, colour = "grey60")


acfe <- acf(evec, pl = FALSE)
acfe <- as.numeric(acfe[1:10]$acf)

p4 <- ggplot(mapping = aes(x = 1:10, y = acfe)) +
  geom_col() +
  ylim(c(-1, 1)) +
  OwnggTheme() +
  xlab("") +
  ylab("") +
  geom_hline(yintercept = 0, colour = "grey60") +
  scale_x_continuous(n.breaks = 6)

plot_grid(p1, p2, p3, p4)

##############################################################################
# Question (h)
# Diagnostic checking auxiliary residuals

# Replication of code from part (c)

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
uvec <- 1 / Fvec * vvec - Ktvec * rvec
Dvec <- 1 / Fvec + Ktvec^2 * Nvec

ustar <- 1 / sqrt(Dvec) * uvec
# Remove last observation as we would otherwise divide by zero 
rstar <- 1 / sqrt(Nvec[-n]) * rvec[-n]

p1 <- ggplot() +
  geom_line(aes(x = year, y = ustar)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  scale_x_continuous(n.breaks = 6) +
  OwnggTheme() +
  xlab("Year") +
  ylab("")

p2 <- ggplot(mapping = aes(x = ustar)) +
  geom_histogram(aes(y = ..density..), fill = "white", colour = "black",
                 binwidth = 0.4) +
  geom_density() +
  OwnggTheme() +
  scale_x_continuous(n.breaks = 6) +
  xlab("") +
  ylab("")

p3 <- ggplot() +
  geom_line(aes(x = year[-n], y = rstar)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  scale_x_continuous(n.breaks = 6) +
  OwnggTheme() +
  xlab("Year") +
  ylab("")

p4 <- ggplot(mapping = aes(x = rstar)) +
  geom_histogram(aes(y = ..density..), fill = "white", colour = "black",
                 binwidth = 0.41) +
  geom_density() +
  OwnggTheme() +
  scale_x_continuous(n.breaks = 6) +
  xlab("") +
  ylab("") 

plot_grid(p1, p2, p3, p4)
