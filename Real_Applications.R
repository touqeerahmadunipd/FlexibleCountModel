rm(list = ls())

# ===========================
# Load Required Packages
# ===========================
library(glmmTMB)
library(MASS)
library(GJRM)
library(pscl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggdist)
library(ggbeeswarm)
library(xtable)

# Install missing packages if necessary
if (!require("pscl")) install.packages("pscl")
if (!require("MASS")) install.packages("MASS")
if (!require("evd")) install.packages("evd")

# ===========================
# Set Working Directory
# ===========================
setwd("~/Documents/Touqeer-code/degpd-and-zidegpd-main")

# ===========================
# Source Custom R Functions
# ===========================
r_files <- list.files("./R/", full.names = TRUE)
lapply(r_files, source)

# ===========================
# Compile C++ Files
# ===========================
cpp_files <- list.files("./src/", full.names = TRUE)
lapply(cpp_files, Rcpp::sourceCpp)

# ===========================
# Load Data
# ===========================
setwd("~/Documents/Touqeer-Docs")
data <- read.csv("globalterrorismdb_0522dist.csv")
head(data)

# ===========================
# Subset for Afghanistan
# ===========================
afghanistan_data <- subset(data, country_txt == "Afghanistan")
head(afghanistan_data)
plot(table(afghanistan_data$nkill))

dat <- data.frame(
  nkill = afghanistan_data$nkillter,
  xcoord = afghanistan_data$longitude,
  ycoord = afghanistan_data$latitude
)
dat <- na.omit(dat)
dat$xcoord_ycoord <- dat$xcoord * dat$ycoord

plot(table(dat$nkill), main = "Distribution of nkill", ylab = "Frequency")

# ===========================
# 1. Poisson Model
# ===========================
poisson_model <- glm(
  nkill ~ xcoord * ycoord,
  family = poisson(link = "log"),
  data = dat
)
summary(poisson_model)

# ===========================
# 2. Negative Binomial Model
# ===========================
nb_model <- glmmTMB(
  nkill ~ xcoord * ycoord,
  data = dat,
  family = nbinom2()
)
summary(nb_model)

# ===========================
# 3. Zero-Inflated Negative Binomial
# ===========================
zinb_model <- glmmTMB(
  nkill ~ xcoord * ycoord,
  ziformula = ~1,
  family = nbinom2(),
  data = dat
)
summary(zinb_model)

# ===========================
# 4. DEGPD Model
# ===========================
fmla_degpd <- list(
  lsigma = nkill ~ xcoord + ycoord + xcoord_ycoord,
  lxi = ~1,
  lkappa = ~1
)

inits_degpd <- c(
  sqrt(sd(dat$nkill) / 2),
  sqrt(sd(dat$nkill) / 20),
  0.2
)

degpd_model <- devgam(
  fmla_degpd,
  data = dat,
  family = "degpd",
  degpd.args = list(m = 1)
)
summary(degpd_model)

# ===========================
# 5. Zero-Inflated DEGPD
# ===========================
zifmla_degpd <- list(
  lsigma = nkill ~ xcoord + ycoord + xcoord_ycoord,
  lxi = ~1,
  lkappa = ~1,
  logitpi = ~1
)

inits_zidegpd <- c(sd(dat$nkill) / 2, 0.5, sd(dat$nkill) / 20, mean(dat$nkill == 0))

zidegpd_model <- devgam(
  zifmla_degpd,
  data = dat,
  family = "zidegpd",
  degpd.args = list(m = 1),
  inits = inits_zidegpd
)
summary(zidegpd_model)

# ===========================
# 6. Compare Models
# ===========================
model_comp <- data.frame(
  Model = c("Poisson", "NegBin", "ZINB", "DEGPD", "ZIDEGPD"),
  AIC = c(AIC(poisson_model), AIC(nb_model), AIC(zinb_model), AIC(degpd_model), AIC(zidegpd_model)),
  BIC = c(BIC(poisson_model), BIC(nb_model), BIC(zinb_model), BIC(degpd_model), BIC(zidegpd_model))
)
print(model_comp)

latex_table <- xtable(
  model_comp,
  caption = "Model comparison based on AIC and BIC",
  label = "tab:model_comp"
)
print(latex_table, include.rownames = FALSE, booktabs = TRUE, comment = FALSE)

# ===========================
# Randomized Quantile Residuals and Q-Q Plots
# ===========================

plot_residuals <- function(residuals, title, file_name) {
  res_df <- data.frame(residuals = residuals)
  p <- ggplot(res_df, aes(sample = residuals)) +
    stat_qq(shape = 20, size = 1.5, color = "gray") +
    stat_qq_line(color = "red", linewidth = 0.5) +
    labs(title = title, x = "Theoretical Quantiles", y = "Empirical Quantiles") +
    theme_minimal()
  print(p)
  #ggsave(file_name, plot = p, width = 5, height = 5, dpi = 300)
}

# Poisson residuals
lambda <- fitted(poisson_model)
y_pois <- dat$nkill
n <- length(y_pois)
p.phat.l <- ppois(y_pois - 1, lambda)
p.phat.u <- ppois(y_pois, lambda)
u <- runif(n, p.phat.l, p.phat.u)
r_pois <- qnorm(u)
r_pois <- r_pois[!is.infinite(r_pois)]
plot_residuals(r_pois, "Q-Q Plot (Poisson)", "qq_plot_poisson_real.pdf")

# Negative Binomial residuals
nb_model <- glm.nb(nkill ~ xcoord + ycoord + xcoord_ycoord, data = dat)
mu_nb <- fitted(nb_model)
theta_nb <- nb_model$theta
y_nb <- dat$nkill
p.phat.l <- pnbinom(y_nb - 1, size = theta_nb, mu = mu_nb)
p.phat.u <- pnbinom(y_nb, size = theta_nb, mu = mu_nb)
u <- runif(n, p.phat.l, p.phat.u)
r_nb <- qnorm(u)
r_nb <- r_nb[!is.infinite(r_nb)]
plot_residuals(r_nb, "Q-Q Plot (Negative Binomial)", "qq_plot_nb_real.pdf")

# ZINB residuals
zinb_model <- zeroinfl(nkill ~ xcoord + ycoord + xcoord_ycoord | 1, dist = "negbin", data = dat)
mu_zinb <- predict(zinb_model, type = "count")
pi_zinb <- predict(zinb_model, type = "zero")
theta_zinb <- zinb_model$theta
y_zinb <- dat$nkill

pzinb <- function(y, mu, theta, pi) {
  sapply(1:length(y), function(i) {
    if (y[i] == 0) pi[i] + (1 - pi[i]) * dnbinom(0, size = theta, mu = mu[i])
    else pi[i] + (1 - pi[i]) * pnbinom(y[i], size = theta, mu = mu[i])
  })
}

p.phat.l <- pzinb(y_zinb - 1, mu_zinb, theta_zinb, pi_zinb)
p.phat.u <- pzinb(y_zinb, mu_zinb, theta_zinb, pi_zinb)
p.phat.l[y_zinb == 0] <- 0
u <- runif(n, p.phat.l, p.phat.u)
r_zinb <- qnorm(u)
r_zinb <- r_zinb[!is.infinite(r_zinb)]
plot_residuals(r_zinb, "Q-Q Plot (ZINB)", "qq_plot_zinb_real.pdf")

# DEGPD residuals
sigma <- exp(degpd_model$logscale$fitted)
shape <- exp(degpd_model$logshape$fitted)
kappa <- exp(degpd_model$logkappa$fitted)
y_degpd <- dat$nkill
p.phat.l <- pdiscegpd(y_degpd - 1, kappa = kappa, sigma = sigma, xi = shape, type = 1)
p.phat.u <- pdiscegpd(y_degpd, kappa = kappa, sigma = sigma, xi = shape, type = 1)
u <- runif(n, p.phat.l, p.phat.u)
r_degpd <- qnorm(u)
r_degpd <- r_degpd[!is.infinite(r_degpd)]
plot_residuals(r_degpd, "Q-Q Plot (DEGPD)", "qq_plot_degpd_real.pdf")

# ZIDEGPD residuals
sigma <- exp(zidegpd_model$logscale$fitted)
shape <- exp(zidegpd_model$logshape$fitted)
kappa <- exp(zidegpd_model$logkappa$fitted)
pi <- exp(zidegpd_model$logitpi$fitted) / (1 + exp(zidegpd_model$logitpi$fitted))
p.phat.l <- pzidiscegpd(y_degpd - 1, kappa = kappa, sigma = sigma, xi = shape, pi = pi, type = 1)
p.phat.u <- pzidiscegpd(y_degpd, kappa = kappa, sigma = sigma, xi = shape, pi = pi, type = 1)
u <- runif(n, p.phat.l, p.phat.u)
r_zidegpd <- qnorm(u)
r_zidegpd <- r_zidegpd[!is.infinite(r_zidegpd)]
plot_residuals(r_zidegpd, "Q-Q Plot (ZIDEGPD)", "qq_plot_zidegpd_real.pdf")
