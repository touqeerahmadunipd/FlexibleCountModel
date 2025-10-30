###############################################################
## Motivating Example
## Simulate data from Negative Binomial distribution with
## different levels of zero inflation and different degrees of outliers
###############################################################

rm(list = ls())

# --- Load Required Packages ---
library(MASS)
library(GJRM)
library(pscl)
library(ggplot2)

# Install missing packages
if (!require("pscl")) install.packages("pscl")
if (!require("MASS")) install.packages("MASS")
if (!require("evd")) install.packages("evd")

library(pscl)
library(MASS)
library(evd)

# --- Set Working Directory ---
setwd("~/Documents/Touqeer-code/degpd-and-zidegpd-main")

# --- Source Custom R Functions ---
dest <- "./R/"
files <- list.files(dest, full.names = TRUE)
for (i in 1:length(files)) {
  source(files[i])
}

# --- Compile C++ Files ---
dest <- "./src/"
files <- list.files(dest, full.names = TRUE)
for (i in seq_along(files)) {
  Rcpp::sourceCpp(files[i])
}

###############################################################
## Function: simulate_one_rep
###############################################################
simulate_one_rep <- function(
    n = 500, 
    zero_inflation = 0.25, 
    outlier_ratio = 0.05, 
    outlier_level = "low", 
    beta = c(1, 0.5, -0.7), 
    mu = 2, 
    theta = 2
) {
  # Step 1: Set outlier values
  outlier_values <- switch(
    outlier_level,
    low = 30:49,
    medium = 50:69,
    high = 70:89
  )
  
  # Step 2: Sample size after excluding outliers
  num_outliers <- ceiling(outlier_ratio * n)
  n1 <- n - num_outliers
  
  # Step 3: Simulate covariates
  X1 <- sample(c(1, 2), n, replace = TRUE, prob = c(0.5, 0.5))
  X2 <- sample(c(1, 2, 3), n, replace = TRUE, prob = c(0.4, 0.35, 0.25))
  
  # Create model matrix
  model_mat <- model.matrix(~ X1 + X2)
  eta <- model_mat %*% beta
  
  # Step 4: Simulate NB data for n1 observations
  mu_nb <- exp(eta[1:n1])
  y_nb <- rnegbin(n1, mu = mu_nb, theta = theta)
  
  # Apply zero inflation
  zero_indices <- sample(1:n1, size = floor(zero_inflation * n1), replace = FALSE)
  y_nb[zero_indices] <- 0
  
  # Step 5: Simulate outliers
  outliers <- sample(outlier_values, size = num_outliers, replace = TRUE)
  
  # Step 6: Combine Y
  Y <- c(y_nb, outliers)
  
  # Step 7: Return data
  return(data.frame(Y = Y, X1 = factor(X1), X2 = factor(X2)))
}

###############################################################
## Simulation Checks
###############################################################
sim_data <- simulate_one_rep(outlier_level = "low")
plot(table(sim_data$Y))

sim_data <- simulate_one_rep(outlier_level = "medium")
plot(table(sim_data$Y))

sim_data <- simulate_one_rep(outlier_level = "high")
plot(table(sim_data$Y))

par(mfrow = c(1, 1))

###############################################################
## Simulate Data for Model Fitting
###############################################################
set.seed(123)
sim_data <- simulate_one_rep(n = 500, zero_inflation = 0.50, outlier_level = "high")

###############################################################
## ============ POISSON MODEL ============
###############################################################
poisson_model <- glm(Y ~ X1 + X2, family = poisson(link = "log"), data = sim_data)
lambda <- fitted(poisson_model)
y_pois <- sim_data$Y
n <- length(y_pois)

# Randomized quantile residuals
p.phat.l <- ppois(y_pois - 1, lambda = lambda)
p.phat.u <- ppois(y_pois, lambda = lambda)
u <- runif(n, p.phat.l, p.phat.u)
r_pois <- qnorm(u)
r_pois <- r_pois[!is.infinite(r_pois)]

qqnorm(r_pois, pch = 20, cex = 0.8, col = "gray", main = "Q-Q Plot (Poisson)")
qqline(r_pois, col = 2, lwd = 1.5)

res_df <- data.frame(residuals = r_pois)
ggplot(res_df, aes(sample = residuals)) +
  stat_qq(shape = 20, size = 1.5, color = "gray") +
  stat_qq_line(color = "red", linewidth = 0.5) +
  labs(title = "Q-Q Plot (Poisson)",
       x = "Theoretical Quantiles",
       y = "Empirical Quantiles") +
  theme_minimal()

# AIC/BIC
aic_pois <- AIC(poisson_model)
bic_pois <- BIC(poisson_model)

###############################################################
## ============ NEGATIVE BINOMIAL MODEL ============
###############################################################
sim_data <- simulate_one_rep(n = 500, zero_inflation = 0.50, outlier_level = "high")
nb_model <- glm.nb(Y ~ X1 + X2, data = sim_data)
mu_nb <- fitted(nb_model)
theta_nb <- nb_model$theta
y_nb <- sim_data$Y

p.phat.l <- pnbinom(y_nb - 1, size = theta_nb, mu = mu_nb)
p.phat.u <- pnbinom(y_nb, size = theta_nb, mu = mu_nb)
u <- runif(n, p.phat.l, p.phat.u)
r_nb <- qnorm(u)
r_nb <- r_nb[!is.infinite(r_nb)]

# Q-Q Plot
res_df <- data.frame(residuals = r_nb)
p <- ggplot(res_df, aes(sample = residuals)) +
  stat_qq(shape = 20, size = 1.5, color = "gray") +
  stat_qq_line(color = "red", linewidth = 0.5) +
  labs(title = "Q-Q Plot (Negative Binomial)",
       x = "Theoretical Quantiles",
       y = "Empirical Quantiles") +
  theme_minimal()
p

# AIC/BIC
aic_nb <- AIC(nb_model)
bic_nb <- BIC(nb_model)

###############################################################
## ============ ZERO-INFLATED NEGATIVE BINOMIAL ============
###############################################################
sim_data <- simulate_one_rep(n = 500, zero_inflation = 0.50, outlier_level = "high")

zinb_model <- zeroinfl(Y ~ X1 + X2 | 1, dist = "negbin", data = sim_data)
mu_zinb <- predict(zinb_model, type = "count")
pi_zinb <- predict(zinb_model, type = "zero")
theta_zinb <- zinb_model$theta
y_zinb <- sim_data$Y

# Define ZINB CDF function
pzinb <- function(y, mu, theta, pi) {
  cdf <- numeric(length(y))
  for (i in 1:length(y)) {
    if (y[i] == 0) {
      cdf[i] <- pi[i] + (1 - pi[i]) * dnbinom(0, size = theta, mu = mu[i])
    } else {
      cdf[i] <- pi[i] + (1 - pi[i]) * pnbinom(y[i], size = theta, mu = mu[i])
    }
  }
  return(cdf)
}

p.phat.l <- pzinb(y_zinb - 1, mu_zinb, theta_zinb, pi_zinb)
p.phat.u <- pzinb(y_zinb, mu_zinb, theta_zinb, pi_zinb)
p.phat.l[y_zinb == 0] <- 0
u <- runif(n, p.phat.l, p.phat.u)
r_zinb <- qnorm(u)
r_zinb <- r_zinb[!is.infinite(r_zinb)]

res_df <- data.frame(residuals = r_zinb)
p <- ggplot(res_df, aes(sample = residuals)) +
  stat_qq(shape = 20, size = 1.5, color = "gray") +
  stat_qq_line(color = "red", linewidth = 0.5) +
  labs(title = "Q-Q Plot (ZINB model with high level outliers)",
       x = "Theoretical Quantiles",
       y = "Empirical Quantiles") +
  theme_minimal()
p
#ggsave("qq_plot_zinb-high.pdf", plot = p, width = 6, height = 5, dpi = 300)

# AIC/BIC
aic_zinb <- AIC(zinb_model)
bic_zinb <- BIC(zinb_model)

###############################################################
## ============ DEGPD MODEL ============
###############################################################
inits1 <- list(logscale = 0.5, logshape = -0.2, logkappa = -0.5)

fmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1)
degpd1 <- devgam(fmla_degpd1, data = sim_data, family = "degpd",
                 degpd.args = list(m = 1), trace = 2)

sigma <- exp(degpd1$logscale$fitted)
shape <- exp(degpd1$logshape$fitted)
kappa <- exp(degpd1$logkappa$fitted)
y_degpd <- sim_data$Y

p.phat.l <- pdiscegpd(y_degpd - 1, kappa = kappa, sigma = sigma, xi = shape, type = 1)
p.phat.u <- pdiscegpd(y_degpd, kappa = kappa, sigma = sigma, xi = shape, type = 1)
u <- runif(n, p.phat.l, p.phat.u)
r_degpd <- qnorm(u)
r_degpd <- r_degpd[!is.infinite(r_degpd)]

res_df <- data.frame(residuals = r_degpd)
p <- ggplot(res_df, aes(sample = residuals)) +
  stat_qq(shape = 20, size = 1.5, color = "gray") +
  stat_qq_line(color = "red", linewidth = 0.5) +
  labs(title = "Q-Q Plot (DEGPD)",
       x = "Theoretical Quantiles",
       y = "Empirical Quantiles") +
  theme_minimal()
p
#ggsave("qq_plot_degpd-low.pdf", plot = p, width = 6, height = 5, dpi = 300)

aic_degpd <- AIC(degpd1)
bic_degpd <- BIC(degpd1)

###############################################################
## ============ ZIDEGPD MODEL ============
###############################################################
inits1 <- list(lsigma = log(1), lxi = log(0.1), lkappa = log(0.5), logitpi = qlogis(0.1))

fmla_zidegpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1, logitpi = ~1)
zidegpd1 <- devgam(fmla_zidegpd1, data = sim_data, family = "zidegpd",
                   degpd.args = list(m = 1), trace = 2,inits1)

sigma <- exp(zidegpd1$logscale$fitted)
shape <- exp(zidegpd1$logshape$fitted)
kappa <- exp(zidegpd1$logkappa$fitted)
pi <- exp(zidegpd1$logitpi$fitted) / (1 + exp(zidegpd1$logitpi$fitted))
y_degpd <- sim_data$Y

p.phat.l <- pzidiscegpd(y_degpd - 1, kappa = kappa, sigma = sigma, xi = shape, pi = pi, type = 1)
p.phat.u <- pzidiscegpd(y_degpd, kappa = kappa, sigma = sigma, xi = shape, pi = pi, type = 1)
u <- runif(n, p.phat.l, p.phat.u)
r_degpd <- qnorm(u)
r_degpd <- r_degpd[!is.infinite(r_degpd)]

res_df <- data.frame(residuals = r_degpd)
p <- ggplot(res_df, aes(sample = residuals)) +
  stat_qq(shape = 20, size = 1.5, color = "gray") +
  stat_qq_line(color = "red", linewidth = 0.5) +
  labs(title = "Q-Q Plot (ZIDEGPD)",
       x = "Theoretical Quantiles",
       y = "Empirical Quantiles") +
  theme_minimal()
p
#ggsave("qq_plot_zidegpd-low.pdf", plot = p, width = 6, height = 5, dpi = 300)

aic_degpd <- AIC(degpd1)
bic_degpd <- BIC(degpd1)
