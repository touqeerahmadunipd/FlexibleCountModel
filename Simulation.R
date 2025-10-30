##Simulate date from Negative binomial with different level of zero inflation and different degrees of outlierds
rm(list=ls())
library(MASS)
library(GJRM)
library(pscl)
# Load required packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggdist)
library(ggbeeswarm)
# Load all required packages
if (!require("pscl")) install.packages("pscl")
if (!require("MASS")) install.packages("MASS")
#if (!require("extRemes")) install.packages("extRemes")
if (!require("evd")) install.packages("evd")
library(pscl)
library(MASS)
#library(extRemes)
setwd("~/Documents/Touqeer-code/degpd-and-zidegpd-main")

dest <- "./R/"      # this function all R function 
files = list.files(dest, full.names = T)
for (i in 1:length(files)) {
  source(files[i])
}


#system("gfortran-15 --version")

# Compile your C++ files
dest <- "./src/"
files <- list.files(dest, full.names = TRUE)
for (i in seq_along(files)) {
  Rcpp::sourceCpp(files[i])
}




#
simulate_one_rep <- function(n = 500, zero_inflation = 0.25, outlier_ratio = 0.05, 
                             outlier_level = "low", beta = c(1, 0.5, -0.7), kappa = 1.5, xi = 0.2) {
  # Step 1: Set outlier values
  outlier_values <- switch(outlier_level,
                           low = 30:49,
                           medium = 50:69,
                           high = 70:89)
  
  # Step 2: Sample size after excluding outliers
  num_outliers <- ceiling(outlier_ratio * n)
  n1 <- n - num_outliers
  
  # Step 3: Simulate covariates
  X1 <- sample(c(1, 2), n, replace = TRUE, prob = c(0.5, 0.5))
  X2 <- sample(c(1, 2, 3), n, replace = TRUE, prob = c(0.4, 0.35, 0.25))
  #X1_factor <- factor(X1)
  #X2_factor <- factor(X2)
  
  # Create model matrix
  model_mat <- model.matrix(~ X1 + X2)
  eta <- model_mat %*% beta
  
  # Step 4: Simulate NB data for n1 observations
  sigma <- exp(eta[1:n1])  # only first n1 rows
  y_degpd <- rdiscegpd(n1, kappa=kappa, sigma = sigma, xi=xi, type = 1)
  
  # Apply zero-inflation
  zero_indices <- sample(1:n1, size = floor(zero_inflation * n1), replace = FALSE)
  y_degpd[zero_indices] <- 0
  
  # Step 5: Simulate outliers
  outliers <- sample(outlier_values, size = num_outliers, replace = TRUE)
  
  # Step 6: Combine Y
  Y <- c(y_degpd, outliers)
  
  # Step 7: Finalize covariates
  X1_final <- X1
  X2_final <- X2
  
  # Return data
  return(data.frame(Y = Y, X1 = factor(X1_final), X2 = factor(X2_final)))
}

#Simulation check with low outliers
sim_data <- simulate_one_rep(outlier_level = "low")
plot(table(sim_data$Y))

#Simulation check with medium outliers
sim_data <- simulate_one_rep(outlier_level = "medium")
plot(table(sim_data$Y))

#Simulation check with high outliers
sim_data <- simulate_one_rep(outlier_level = "high")
plot(table(sim_data$Y))



###. == Low outliers with 0.01 ratio and 25% zeros

# Initialize storage matrices
results_matrix <- matrix(NA, nrow = 1000, ncol = 10)
colnames(results_matrix) <- c("Poisson_AIC", "Poisson_BIC",
                              "NB_AIC", "NB_BIC",
                              "ZINB_AIC", "ZINB_BIC",
                              "DEGPD_AIC", "DEGPD_BIC",
                              "ZIDEGPD_AIC", "ZIDEGPD_BIC")

set.seed(123)

for (i in 1:1000) {
  # Simulate data once
  sim_data <- simulate_one_rep(outlier_level = "low", zero_inflation = 0.25, outlier_ratio = 0.01)
  
  # 1. Poisson Model
  poisson_model <- try(glm(Y ~ X1 + X2, family = poisson(link = "log"), data = sim_data), silent = TRUE)
  if (!inherits(poisson_model, "try-error")) {
    results_matrix[i, 1] <- AIC(poisson_model)
    results_matrix[i, 2] <- BIC(poisson_model)
  }
  
  # 2. Negative Binomial Model
  nb_model <- try(glm.nb(Y ~ X1 + X2, data = sim_data), silent = TRUE)
  if (!inherits(nb_model, "try-error")) {
    results_matrix[i, 3] <- AIC(nb_model)
    results_matrix[i, 4] <- BIC(nb_model)
  }
  
  # 3. ZINB Model
  zinb_model <- try(zeroinfl(Y ~ X1 + X2 | X1 + X2, data = sim_data, dist = "negbin"), silent = TRUE)
  if (!inherits(zinb_model, "try-error")) {
    results_matrix[i, 5] <- AIC(zinb_model)
    results_matrix[i, 6] <- BIC(zinb_model)
  }
  
  # 4. DEGPD Model
  fmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1)
  # fit<-fit.model(sim_data$Y, model = 1, family = c("discegpd"), init = c(1,1,0.2))
  # inits_degpd=fit$fit$mle
  
  inits_degpd <- c(sqrt(sd(sim_data$Y)/2), sqrt(sd(sim_data$Y)/20), 0.3)
  degpd1 <- try(devgam(fmla_degpd1, data = sim_data, family = "degpd",
                       degpd.args = list(m = 1), trace = 0), silent = TRUE)
  if (!inherits(degpd1, "try-error")) {
    results_matrix[i, 7] <- AIC(degpd1)
    results_matrix[i, 8] <- BIC(degpd1)
  }
  
  # 5. ZIDEGPD Model
  zifmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1, logitpi = ~1)
  inits_zidegpd <- c(sd(sim_data$Y)/2, sd(sim_data$Y)/2, 0.2, mean(sim_data$Y == 0))
  #fit<-fit.model(sim_data$Y, model = 1, family = c("zidiscegpd"), init = c(0.2,1,1,0.2))
  #inits_degpd=fit$fit$mle
  zidegpd1 <- try(devgam(zifmla_degpd1, data = sim_data, family = "zidegpd",
                         degpd.args = list(m = 1), trace = 0, inits = inits_zidegpd), silent = TRUE)
  if (!inherits(zidegpd1, "try-error")) {
    results_matrix[i, 9] <- AIC(zidegpd1)
    results_matrix[i, 10] <- BIC(zidegpd1)
  }
  
  if (i %% 100 == 0) cat("Completed:", i, "\n")
}




# Compute means for each model
mean_aic_bic <- data.frame(
  Model = c("Poisson", "Negative Binomial", "ZINB", "DEGPD", "ZIDEGPD"),
  Mean_AIC = c(
    mean(results_matrix[, 1], na.rm = TRUE),
    mean(results_matrix[, 3], na.rm = TRUE),
    mean(results_matrix[, 5], na.rm = TRUE),
    mean(results_matrix[, 7], na.rm = TRUE),
    mean(results_matrix[, 9], na.rm = TRUE)
  ),
  Mean_BIC = c(
    mean(results_matrix[, 2], na.rm = TRUE),
    mean(results_matrix[, 4], na.rm = TRUE),
    mean(results_matrix[, 6], na.rm = TRUE),
    mean(results_matrix[, 8], na.rm = TRUE),
    mean(results_matrix[, 10], na.rm = TRUE)
  )
)

# Print summary table
print(mean_aic_bic, row.names = FALSE)

# Combine AIC/BIC matrices into a long-format dataframe
aic_bic_df <- data.frame(
  Model = rep(c( "NB", "ZINB", "DEGPD", "ZIDEGPD"), each = nrow(results_matrix)),
  AIC = c(#results_matrix[, 1],
    results_matrix[, 3],
    results_matrix[, 5],
    results_matrix[, 7],
    results_matrix[, 9]),
  BIC = c(#results_matrix[, 2],
    results_matrix[, 4],
    results_matrix[, 6],
    results_matrix[, 8],
    results_matrix[, 10])
)

# Convert to long format
aic_bic_long <- pivot_longer(aic_bic_df, cols = c(AIC, BIC), names_to = "Criterion", values_to = "Value")

# Reorder Model to show DEGPD and ZIDEGPD on the left
aic_bic_long$Model <- factor(aic_bic_long$Model,
                             levels = c("DEGPD", "ZIDEGPD", "NB", "ZINB")
)

# Compute mean of the best model (lowest minimum) per criterion
min_per_model <- aic_bic_long %>%
  group_by(Model, Criterion) %>%
  summarise(Min_Value = min(Value, na.rm = TRUE), .groups = "drop")

lowest_model_per_criterion <- min_per_model %>%
  group_by(Criterion) %>%
  filter(Min_Value == min(Min_Value)) %>%
  select(Criterion, Model)

mean_value_of_best_model <- aic_bic_long %>%
  inner_join(lowest_model_per_criterion, by = c("Criterion", "Model")) %>%
  group_by(Criterion) %>%
  summarise(Line_Value = mean(Value, na.rm = TRUE), .groups = "drop")

# ======= 1. Flask-Shaped Boxplot (varwidth) =======
ggplot(aic_bic_long, aes(x = Model, y = Value, fill = Model)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, varwidth = TRUE) +
  #geom_hline(data = mean_value_of_best_model,
  #aes(yintercept = Line_Value),
  #color = "red", linetype = "dashed", linewidth = 0.5,
  # inherit.aes = FALSE) +
  facet_wrap(~ Criterion, scales = "free_y") +
  labs(title = "Flask-Shaped Boxplots (Varwidth)\nDashed Red Line = Mean of Best Model",
       x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# ======= 2. Violin + Boxplot =======
p<-ggplot(aic_bic_long, aes(x = Model, y = Value, fill = Model)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.4, outlier.shape = TRUE, fill = "white", color = "black", linewidth = 0.2, alpha=0.6) +
  #geom_hline(data = mean_value_of_best_model,
  #           aes(yintercept = Line_Value),
  #          color = "red", linetype = "dashed", linewidth = 0.3,inherit.aes = FALSE) +
  facet_wrap(~ Criterion, scales = "free_y") +
  labs(title = "AIC and BIC plots of the fitted models at 1% low outliers",
       x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

p




# Save to PDF
ggsave("AIC_BIC_0.01_low_0.25_zero_degpd.pdf", plot = p, width = 6, height = 5)








###. == medium outliers with 0.01 ratio and 25% zeros

# Initialize storage matrices
results_matrix <- matrix(NA, nrow = 1000, ncol = 10)
colnames(results_matrix) <- c("Poisson_AIC", "Poisson_BIC",
                              "NB_AIC", "NB_BIC",
                              "ZINB_AIC", "ZINB_BIC",
                              "DEGPD_AIC", "DEGPD_BIC",
                              "ZIDEGPD_AIC", "ZIDEGPD_BIC")

set.seed(123)

for (i in 1:1000) {
  # Simulate data once
  sim_data <- simulate_one_rep(outlier_level = "medium", zero_inflation = 0.25, outlier_ratio = 0.01)
  
  # 1. Poisson Model
  poisson_model <- try(glm(Y ~ X1 + X2, family = poisson(link = "log"), data = sim_data), silent = TRUE)
  if (!inherits(poisson_model, "try-error")) {
    results_matrix[i, 1] <- AIC(poisson_model)
    results_matrix[i, 2] <- BIC(poisson_model)
  }
  
  # 2. Negative Binomial Model
  nb_model <- try(glm.nb(Y ~ X1 + X2, data = sim_data), silent = TRUE)
  if (!inherits(nb_model, "try-error")) {
    results_matrix[i, 3] <- AIC(nb_model)
    results_matrix[i, 4] <- BIC(nb_model)
  }
  
  # 3. ZINB Model
  zinb_model <- try(zeroinfl(Y ~ X1 + X2 | X1 + X2, data = sim_data, dist = "negbin"), silent = TRUE)
  if (!inherits(zinb_model, "try-error")) {
    results_matrix[i, 5] <- AIC(zinb_model)
    results_matrix[i, 6] <- BIC(zinb_model)
  }
  
  # 4. DEGPD Model
  fmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1)
  inits_degpd <- c(sqrt(sd(sim_data$Y)/2), sqrt(sd(sim_data$Y)/20), 0.3)
  degpd1 <- try(devgam(fmla_degpd1, data = sim_data, family = "degpd",
                       degpd.args = list(m = 1), trace = 0), silent = TRUE)
  if (!inherits(degpd1, "try-error")) {
    results_matrix[i, 7] <- AIC(degpd1)
    results_matrix[i, 8] <- BIC(degpd1)
  }
  
  # 5. ZIDEGPD Model
  zifmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1, logitpi = ~1)
  inits_zidegpd <- c(sd(sim_data$Y)/2, sd(sim_data$Y)/2, 0.2, mean(sim_data$Y == 0))
  zidegpd1 <- try(devgam(zifmla_degpd1, data = sim_data, family = "zidegpd",
                         degpd.args = list(m = 1), trace = 0, inits = inits_zidegpd), silent = TRUE)
  if (!inherits(zidegpd1, "try-error")) {
    results_matrix[i, 9] <- AIC(zidegpd1)
    results_matrix[i, 10] <- BIC(zidegpd1)
  }
  
  if (i %% 100 == 0) cat("Completed:", i, "\n")
}


# Compute means for each model
mean_aic_bic <- data.frame(
  Model = c("Poisson", "Negative Binomial", "ZINB", "DEGPD", "ZIDEGPD"),
  Mean_AIC = c(
    mean(results_matrix[, 1], na.rm = TRUE),
    mean(results_matrix[, 3], na.rm = TRUE),
    mean(results_matrix[, 5], na.rm = TRUE),
    mean(results_matrix[, 7], na.rm = TRUE),
    mean(results_matrix[, 9], na.rm = TRUE)
  ),
  Mean_BIC = c(
    mean(results_matrix[, 2], na.rm = TRUE),
    mean(results_matrix[, 4], na.rm = TRUE),
    mean(results_matrix[, 6], na.rm = TRUE),
    mean(results_matrix[, 8], na.rm = TRUE),
    mean(results_matrix[, 10], na.rm = TRUE)
  )
)

# Print summary table
print(mean_aic_bic, row.names = FALSE)

# Combine AIC/BIC matrices into a long-format dataframe
aic_bic_df <- data.frame(
  Model = rep(c( "NB", "ZINB", "DEGPD", "ZIDEGPD"), each = nrow(results_matrix)),
  AIC = c(#results_matrix[, 1],
    results_matrix[, 3],
    results_matrix[, 5],
    results_matrix[, 7],
    results_matrix[, 9]),
  BIC = c(#results_matrix[, 2],
    results_matrix[, 4],
    results_matrix[, 6],
    results_matrix[, 8],
    results_matrix[, 10])
)

# Convert to long format
aic_bic_long <- pivot_longer(aic_bic_df, cols = c(AIC, BIC), names_to = "Criterion", values_to = "Value")

# Reorder Model to show DEGPD and ZIDEGPD on the left
aic_bic_long$Model <- factor(aic_bic_long$Model,
                             levels = c("DEGPD", "ZIDEGPD", "NB", "ZINB")
)

# Compute mean of the best model (lowest minimum) per criterion
min_per_model <- aic_bic_long %>%
  group_by(Model, Criterion) %>%
  summarise(Min_Value = min(Value, na.rm = TRUE), .groups = "drop")

lowest_model_per_criterion <- min_per_model %>%
  group_by(Criterion) %>%
  filter(Min_Value == min(Min_Value)) %>%
  select(Criterion, Model)

mean_value_of_best_model <- aic_bic_long %>%
  inner_join(lowest_model_per_criterion, by = c("Criterion", "Model")) %>%
  group_by(Criterion) %>%
  summarise(Line_Value = mean(Value, na.rm = TRUE), .groups = "drop")

# ======= 1. Flask-Shaped Boxplot (varwidth) =======
ggplot(aic_bic_long, aes(x = Model, y = Value, fill = Model)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, varwidth = TRUE) +
  #geom_hline(data = mean_value_of_best_model,
  #aes(yintercept = Line_Value),
  #color = "red", linetype = "dashed", linewidth = 0.5,
  # inherit.aes = FALSE) +
  facet_wrap(~ Criterion, scales = "free_y") +
  labs(title = "Flask-Shaped Boxplots (Varwidth)\nDashed Red Line = Mean of Best Model",
       x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# ======= 2. Violin + Boxplot =======
p<-ggplot(aic_bic_long, aes(x = Model, y = Value, fill = Model)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.4, outlier.shape = TRUE, fill = "white", color = "black", linewidth = 0.2, alpha=0.6) +
  #geom_hline(data = mean_value_of_best_model,
  #           aes(yintercept = Line_Value),
  #          color = "red", linetype = "dashed", linewidth = 0.3,inherit.aes = FALSE) +
  facet_wrap(~ Criterion, scales = "free_y") +
  labs(title = "AIC and BIC plots of the fitted models at 1% medium outliers",
       x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

p

# Save to PDF
ggsave("AIC_BIC_0.01_medium_0.25_zero_degpd.pdf", plot = p, width = 6, height = 5)




###. == high outliers with 0.01 ratio and 25% zeros. ===

# Initialize storage matrices
results_matrix <- matrix(NA, nrow = 1000, ncol = 10)
colnames(results_matrix) <- c("Poisson_AIC", "Poisson_BIC",
                              "NB_AIC", "NB_BIC",
                              "ZINB_AIC", "ZINB_BIC",
                              "DEGPD_AIC", "DEGPD_BIC",
                              "ZIDEGPD_AIC", "ZIDEGPD_BIC")

set.seed(123)

for (i in 1:1000) {
  # Simulate data once
  sim_data <- simulate_one_rep(outlier_level = "high", zero_inflation = 0.25, outlier_ratio = 0.01)
  
  # 1. Poisson Model
  poisson_model <- try(glm(Y ~ X1 + X2, family = poisson(link = "log"), data = sim_data), silent = TRUE)
  if (!inherits(poisson_model, "try-error")) {
    results_matrix[i, 1] <- AIC(poisson_model)
    results_matrix[i, 2] <- BIC(poisson_model)
  }
  
  # 2. Negative Binomial Model
  nb_model <- try(glm.nb(Y ~ X1 + X2, data = sim_data), silent = TRUE)
  if (!inherits(nb_model, "try-error")) {
    results_matrix[i, 3] <- AIC(nb_model)
    results_matrix[i, 4] <- BIC(nb_model)
  }
  
  # 3. ZINB Model
  zinb_model <- try(zeroinfl(Y ~ X1 + X2 | X1 + X2, data = sim_data, dist = "negbin"), silent = TRUE)
  if (!inherits(zinb_model, "try-error")) {
    results_matrix[i, 5] <- AIC(zinb_model)
    results_matrix[i, 6] <- BIC(zinb_model)
  }
  
  # 4. DEGPD Model
  fmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1)
  inits_degpd <- c(sqrt(sd(sim_data$Y)/2), sqrt(sd(sim_data$Y)/20), 0.3)
  degpd1 <- try(devgam(fmla_degpd1, data = sim_data, family = "degpd",
                       degpd.args = list(m = 1), trace = 0), silent = TRUE)
  if (!inherits(degpd1, "try-error")) {
    results_matrix[i, 7] <- AIC(degpd1)
    results_matrix[i, 8] <- BIC(degpd1)
  }
  
  # 5. ZIDEGPD Model
  zifmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1, logitpi = ~1)
  inits_zidegpd <- c(sd(sim_data$Y)/2, sd(sim_data$Y)/2, 0.2, mean(sim_data$Y == 0))
  zidegpd1 <- try(devgam(zifmla_degpd1, data = sim_data, family = "zidegpd",
                         degpd.args = list(m = 1), trace = 0, inits = inits_zidegpd), silent = TRUE)
  if (!inherits(zidegpd1, "try-error")) {
    results_matrix[i, 9] <- AIC(zidegpd1)
    results_matrix[i, 10] <- BIC(zidegpd1)
  }
  
  if (i %% 100 == 0) cat("Completed:", i, "\n")
}



# Compute means for each model
mean_aic_bic <- data.frame(
  Model = c("Poisson", "Negative Binomial", "ZINB", "DEGPD", "ZIDEGPD"),
  Mean_AIC = c(
    mean(results_matrix[, 1], na.rm = TRUE),
    mean(results_matrix[, 3], na.rm = TRUE),
    mean(results_matrix[, 5], na.rm = TRUE),
    mean(results_matrix[, 7], na.rm = TRUE),
    mean(results_matrix[, 9], na.rm = TRUE)
  ),
  Mean_BIC = c(
    mean(results_matrix[, 2], na.rm = TRUE),
    mean(results_matrix[, 4], na.rm = TRUE),
    mean(results_matrix[, 6], na.rm = TRUE),
    mean(results_matrix[, 8], na.rm = TRUE),
    mean(results_matrix[, 10], na.rm = TRUE)
  )
)

# Print summary table
print(mean_aic_bic, row.names = FALSE)

# Combine AIC/BIC matrices into a long-format dataframe
aic_bic_df <- data.frame(
  Model = rep(c( "NB", "ZINB", "DEGPD", "ZIDEGPD"), each = nrow(results_matrix)),
  AIC = c(#results_matrix[, 1],
    results_matrix[, 3],
    results_matrix[, 5],
    results_matrix[, 7],
    results_matrix[, 9]),
  BIC = c(#results_matrix[, 2],
    results_matrix[, 4],
    results_matrix[, 6],
    results_matrix[, 8],
    results_matrix[, 10])
)

# Convert to long format
aic_bic_long <- pivot_longer(aic_bic_df, cols = c(AIC, BIC), names_to = "Criterion", values_to = "Value")

# Reorder Model to show DEGPD and ZIDEGPD on the left
aic_bic_long$Model <- factor(aic_bic_long$Model,
                             levels = c("DEGPD", "ZIDEGPD", "NB", "ZINB")
)

# Compute mean of the best model (lowest minimum) per criterion
min_per_model <- aic_bic_long %>%
  group_by(Model, Criterion) %>%
  summarise(Min_Value = min(Value, na.rm = TRUE), .groups = "drop")

lowest_model_per_criterion <- min_per_model %>%
  group_by(Criterion) %>%
  filter(Min_Value == min(Min_Value)) %>%
  select(Criterion, Model)

mean_value_of_best_model <- aic_bic_long %>%
  inner_join(lowest_model_per_criterion, by = c("Criterion", "Model")) %>%
  group_by(Criterion) %>%
  summarise(Line_Value = mean(Value, na.rm = TRUE), .groups = "drop")

# ======= 1. Flask-Shaped Boxplot (varwidth) =======
ggplot(aic_bic_long, aes(x = Model, y = Value, fill = Model)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, varwidth = TRUE) +
  #geom_hline(data = mean_value_of_best_model,
  #aes(yintercept = Line_Value),
  #color = "red", linetype = "dashed", linewidth = 0.5,
  # inherit.aes = FALSE) +
  facet_wrap(~ Criterion, scales = "free_y") +
  labs(title = "Flask-Shaped Boxplots (Varwidth)\nDashed Red Line = Mean of Best Model",
       x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# ======= 2. Violin + Boxplot =======
p<-ggplot(aic_bic_long, aes(x = Model, y = Value, fill = Model)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.4, outlier.shape = TRUE, fill = "white", color = "black", linewidth = 0.2, alpha=0.6) +
  #geom_hline(data = mean_value_of_best_model,
  #           aes(yintercept = Line_Value),
  #          color = "red", linetype = "dashed", linewidth = 0.3,inherit.aes = FALSE) +
  facet_wrap(~ Criterion, scales = "free_y") +
  labs(title = "AIC and BIC plots of the fitted models at 1% high outliers",
       x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

p




# Save to PDF
ggsave("AIC_BIC_0.01_high_0.25_zero_degpd.pdf", plot = p, width = 6, height = 5)






##==============================================================================================================================================






###. == Low outliers with 5% ratio and 25% zeros

# Initialize storage matrices
results_matrix <- matrix(NA, nrow = 1000, ncol = 10)
colnames(results_matrix) <- c("Poisson_AIC", "Poisson_BIC",
                              "NB_AIC", "NB_BIC",
                              "ZINB_AIC", "ZINB_BIC",
                              "DEGPD_AIC", "DEGPD_BIC",
                              "ZIDEGPD_AIC", "ZIDEGPD_BIC")

set.seed(123)

for (i in 1:1000) {
  # Simulate data once
  sim_data <- simulate_one_rep(outlier_level = "low", zero_inflation = 0.25, outlier_ratio = 0.05)
  
  # 1. Poisson Model
  poisson_model <- try(glm(Y ~ X1 + X2, family = poisson(link = "log"), data = sim_data), silent = TRUE)
  if (!inherits(poisson_model, "try-error")) {
    results_matrix[i, 1] <- AIC(poisson_model)
    results_matrix[i, 2] <- BIC(poisson_model)
  }
  
  # 2. Negative Binomial Model
  nb_model <- try(glm.nb(Y ~ X1 + X2, data = sim_data), silent = TRUE)
  if (!inherits(nb_model, "try-error")) {
    results_matrix[i, 3] <- AIC(nb_model)
    results_matrix[i, 4] <- BIC(nb_model)
  }
  
  # 3. ZINB Model
  zinb_model <- try(zeroinfl(Y ~ X1 + X2 | X1 + X2, data = sim_data, dist = "negbin"), silent = TRUE)
  if (!inherits(zinb_model, "try-error")) {
    results_matrix[i, 5] <- AIC(zinb_model)
    results_matrix[i, 6] <- BIC(zinb_model)
  }
  
  # 4. DEGPD Model
  fmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1)
  inits_degpd <- c(sqrt(sd(sim_data$Y)/2), sqrt(sd(sim_data$Y)/20), 0.3)
  degpd1 <- try(devgam(fmla_degpd1, data = sim_data, family = "degpd",
                       degpd.args = list(m = 1), trace = 0), silent = TRUE)
  if (!inherits(degpd1, "try-error")) {
    results_matrix[i, 7] <- AIC(degpd1)
    results_matrix[i, 8] <- BIC(degpd1)
  }
  
  # 5. ZIDEGPD Model
  zifmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1, logitpi = ~1)
  inits_zidegpd <- c(sd(sim_data$Y)/2, sd(sim_data$Y)/2, 0.2, mean(sim_data$Y == 0))
  zidegpd1 <- try(devgam(zifmla_degpd1, data = sim_data, family = "zidegpd",
                         degpd.args = list(m = 1), trace = 0, inits = inits_zidegpd), silent = TRUE)
  if (!inherits(zidegpd1, "try-error")) {
    results_matrix[i, 9] <- AIC(zidegpd1)
    results_matrix[i, 10] <- BIC(zidegpd1)
  }
  
  if (i %% 100 == 0) cat("Completed:", i, "\n")
}







# Load required packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggdist)
library(ggbeeswarm)

# Compute means for each model
mean_aic_bic <- data.frame(
  Model = c("Poisson", "Negative Binomial", "ZINB", "DEGPD", "ZIDEGPD"),
  Mean_AIC = c(
    mean(results_matrix[, 1], na.rm = TRUE),
    mean(results_matrix[, 3], na.rm = TRUE),
    mean(results_matrix[, 5], na.rm = TRUE),
    mean(results_matrix[, 7], na.rm = TRUE),
    mean(results_matrix[, 9], na.rm = TRUE)
  ),
  Mean_BIC = c(
    mean(results_matrix[, 2], na.rm = TRUE),
    mean(results_matrix[, 4], na.rm = TRUE),
    mean(results_matrix[, 6], na.rm = TRUE),
    mean(results_matrix[, 8], na.rm = TRUE),
    mean(results_matrix[, 10], na.rm = TRUE)
  )
)

# Print summary table
print(mean_aic_bic, row.names = FALSE)

# Combine AIC/BIC matrices into a long-format dataframe
aic_bic_df <- data.frame(
  Model = rep(c( "NB", "ZINB", "DEGPD", "ZIDEGPD"), each = nrow(results_matrix)),
  AIC = c(#results_matrix[, 1],
    results_matrix[, 3],
    results_matrix[, 5],
    results_matrix[, 7],
    results_matrix[, 9]),
  BIC = c(#results_matrix[, 2],
    results_matrix[, 4],
    results_matrix[, 6],
    results_matrix[, 8],
    results_matrix[, 10])
)

# Convert to long format
aic_bic_long <- pivot_longer(aic_bic_df, cols = c(AIC, BIC), names_to = "Criterion", values_to = "Value")

# Reorder Model to show DEGPD and ZIDEGPD on the left
aic_bic_long$Model <- factor(aic_bic_long$Model,
                             levels = c("DEGPD", "ZIDEGPD", "NB", "ZINB")
)

# Compute mean of the best model (lowest minimum) per criterion
min_per_model <- aic_bic_long %>%
  group_by(Model, Criterion) %>%
  summarise(Min_Value = min(Value, na.rm = TRUE), .groups = "drop")

lowest_model_per_criterion <- min_per_model %>%
  group_by(Criterion) %>%
  filter(Min_Value == min(Min_Value)) %>%
  select(Criterion, Model)

mean_value_of_best_model <- aic_bic_long %>%
  inner_join(lowest_model_per_criterion, by = c("Criterion", "Model")) %>%
  group_by(Criterion) %>%
  summarise(Line_Value = mean(Value, na.rm = TRUE), .groups = "drop")

# ======= 1. Flask-Shaped Boxplot (varwidth) =======
ggplot(aic_bic_long, aes(x = Model, y = Value, fill = Model)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, varwidth = TRUE) +
  #geom_hline(data = mean_value_of_best_model,
  #aes(yintercept = Line_Value),
  #color = "red", linetype = "dashed", linewidth = 0.5,
  # inherit.aes = FALSE) +
  facet_wrap(~ Criterion, scales = "free_y") +
  labs(title = "Flask-Shaped Boxplots (Varwidth)\nDashed Red Line = Mean of Best Model",
       x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# ======= 2. Violin + Boxplot =======
p<-ggplot(aic_bic_long, aes(x = Model, y = Value, fill = Model)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.4, outlier.shape = TRUE, fill = "white", color = "black", linewidth = 0.2, alpha=0.6) +
  #geom_hline(data = mean_value_of_best_model,
  #           aes(yintercept = Line_Value),
  #          color = "red", linetype = "dashed", linewidth = 0.3,inherit.aes = FALSE) +
  facet_wrap(~ Criterion, scales = "free_y") +
  labs(title = "AIC and BIC plots of the fitted models at 5% low outliers",
       x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

p




# Save to PDF
ggsave("AIC_BIC_0.05_low_0.25_zero_degpd.pdf", plot = p, width = 6, height = 5)








###. == medium outliers with 0.05 ratio and 25% zeros

# Initialize storage matrices
results_matrix <- matrix(NA, nrow = 1000, ncol = 10)
colnames(results_matrix) <- c("Poisson_AIC", "Poisson_BIC",
                              "NB_AIC", "NB_BIC",
                              "ZINB_AIC", "ZINB_BIC",
                              "DEGPD_AIC", "DEGPD_BIC",
                              "ZIDEGPD_AIC", "ZIDEGPD_BIC")

set.seed(123)

for (i in 1:1000) {
  # Simulate data once
  sim_data <- simulate_one_rep(outlier_level = "medium", zero_inflation = 0.25, outlier_ratio = 0.05)
  
  # 1. Poisson Model
  poisson_model <- try(glm(Y ~ X1 + X2, family = poisson(link = "log"), data = sim_data), silent = TRUE)
  if (!inherits(poisson_model, "try-error")) {
    results_matrix[i, 1] <- AIC(poisson_model)
    results_matrix[i, 2] <- BIC(poisson_model)
  }
  
  # 2. Negative Binomial Model
  nb_model <- try(glm.nb(Y ~ X1 + X2, data = sim_data), silent = TRUE)
  if (!inherits(nb_model, "try-error")) {
    results_matrix[i, 3] <- AIC(nb_model)
    results_matrix[i, 4] <- BIC(nb_model)
  }
  
  # 3. ZINB Model
  zinb_model <- try(zeroinfl(Y ~ X1 + X2 | X1 + X2, data = sim_data, dist = "negbin"), silent = TRUE)
  if (!inherits(zinb_model, "try-error")) {
    results_matrix[i, 5] <- AIC(zinb_model)
    results_matrix[i, 6] <- BIC(zinb_model)
  }
  
  # 4. DEGPD Model
  fmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1)
  inits_degpd <- c(sqrt(sd(sim_data$Y)/2), sqrt(sd(sim_data$Y)/20), 0.3)
  degpd1 <- try(devgam(fmla_degpd1, data = sim_data, family = "degpd",
                       degpd.args = list(m = 1), trace = 0), silent = TRUE)
  if (!inherits(degpd1, "try-error")) {
    results_matrix[i, 7] <- AIC(degpd1)
    results_matrix[i, 8] <- BIC(degpd1)
  }
  
  # 5. ZIDEGPD Model
  zifmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1, logitpi = ~1)
  inits_zidegpd <- c(sd(sim_data$Y)/2, sd(sim_data$Y)/2, 0.2, mean(sim_data$Y == 0))
  zidegpd1 <- try(devgam(zifmla_degpd1, data = sim_data, family = "zidegpd",
                         degpd.args = list(m = 1), trace = 0, inits = inits_zidegpd), silent = TRUE)
  if (!inherits(zidegpd1, "try-error")) {
    results_matrix[i, 9] <- AIC(zidegpd1)
    results_matrix[i, 10] <- BIC(zidegpd1)
  }
  
  if (i %% 100 == 0) cat("Completed:", i, "\n")
}






# Compute means for each model
mean_aic_bic <- data.frame(
  Model = c("Poisson", "Negative Binomial", "ZINB", "DEGPD", "ZIDEGPD"),
  Mean_AIC = c(
    mean(results_matrix[, 1], na.rm = TRUE),
    mean(results_matrix[, 3], na.rm = TRUE),
    mean(results_matrix[, 5], na.rm = TRUE),
    mean(results_matrix[, 7], na.rm = TRUE),
    mean(results_matrix[, 9], na.rm = TRUE)
  ),
  Mean_BIC = c(
    mean(results_matrix[, 2], na.rm = TRUE),
    mean(results_matrix[, 4], na.rm = TRUE),
    mean(results_matrix[, 6], na.rm = TRUE),
    mean(results_matrix[, 8], na.rm = TRUE),
    mean(results_matrix[, 10], na.rm = TRUE)
  )
)

# Print summary table
print(mean_aic_bic, row.names = FALSE)

# Combine AIC/BIC matrices into a long-format dataframe
aic_bic_df <- data.frame(
  Model = rep(c( "NB", "ZINB", "DEGPD", "ZIDEGPD"), each = nrow(results_matrix)),
  AIC = c(#results_matrix[, 1],
    results_matrix[, 3],
    results_matrix[, 5],
    results_matrix[, 7],
    results_matrix[, 9]),
  BIC = c(#results_matrix[, 2],
    results_matrix[, 4],
    results_matrix[, 6],
    results_matrix[, 8],
    results_matrix[, 10])
)

# Convert to long format
aic_bic_long <- pivot_longer(aic_bic_df, cols = c(AIC, BIC), names_to = "Criterion", values_to = "Value")

# Reorder Model to show DEGPD and ZIDEGPD on the left
aic_bic_long$Model <- factor(aic_bic_long$Model,
                             levels = c("DEGPD", "ZIDEGPD", "NB", "ZINB")
)

# Compute mean of the best model (lowest minimum) per criterion
min_per_model <- aic_bic_long %>%
  group_by(Model, Criterion) %>%
  summarise(Min_Value = min(Value, na.rm = TRUE), .groups = "drop")

lowest_model_per_criterion <- min_per_model %>%
  group_by(Criterion) %>%
  filter(Min_Value == min(Min_Value)) %>%
  select(Criterion, Model)

mean_value_of_best_model <- aic_bic_long %>%
  inner_join(lowest_model_per_criterion, by = c("Criterion", "Model")) %>%
  group_by(Criterion) %>%
  summarise(Line_Value = mean(Value, na.rm = TRUE), .groups = "drop")

# ======= 1. Flask-Shaped Boxplot (varwidth) =======
ggplot(aic_bic_long, aes(x = Model, y = Value, fill = Model)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, varwidth = TRUE) +
  #geom_hline(data = mean_value_of_best_model,
  #aes(yintercept = Line_Value),
  #color = "red", linetype = "dashed", linewidth = 0.5,
  # inherit.aes = FALSE) +
  facet_wrap(~ Criterion, scales = "free_y") +
  labs(title = "Flask-Shaped Boxplots (Varwidth)\nDashed Red Line = Mean of Best Model",
       x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# ======= 2. Violin + Boxplot =======
p<-ggplot(aic_bic_long, aes(x = Model, y = Value, fill = Model)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.4, outlier.shape = TRUE, fill = "white", color = "black", linewidth = 0.2, alpha=0.6) +
  #geom_hline(data = mean_value_of_best_model,
  #           aes(yintercept = Line_Value),
  #          color = "red", linetype = "dashed", linewidth = 0.3,inherit.aes = FALSE) +
  facet_wrap(~ Criterion, scales = "free_y") +
  labs(title = "AIC and BIC plots of the fitted models at 5% medium outliers",
       x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

p




# Save to PDF
ggsave("AIC_BIC_0.05_medium_0.25_zero_degpd.pdf", plot = p, width = 6, height = 5)




###. == high outliers with 0.01 ratio and 25% zeros. ===

# Initialize storage matrices
results_matrix <- matrix(NA, nrow = 1000, ncol = 10)
colnames(results_matrix) <- c("Poisson_AIC", "Poisson_BIC",
                              "NB_AIC", "NB_BIC",
                              "ZINB_AIC", "ZINB_BIC",
                              "DEGPD_AIC", "DEGPD_BIC",
                              "ZIDEGPD_AIC", "ZIDEGPD_BIC")

set.seed(123)

for (i in 1:1000) {
  # Simulate data once
  sim_data <- simulate_one_rep(outlier_level = "high", zero_inflation = 0.25, outlier_ratio = 0.05)
  
  # 1. Poisson Model
  poisson_model <- try(glm(Y ~ X1 + X2, family = poisson(link = "log"), data = sim_data), silent = TRUE)
  if (!inherits(poisson_model, "try-error")) {
    results_matrix[i, 1] <- AIC(poisson_model)
    results_matrix[i, 2] <- BIC(poisson_model)
  }
  
  # 2. Negative Binomial Model
  nb_model <- try(glm.nb(Y ~ X1 + X2, data = sim_data), silent = TRUE)
  if (!inherits(nb_model, "try-error")) {
    results_matrix[i, 3] <- AIC(nb_model)
    results_matrix[i, 4] <- BIC(nb_model)
  }
  
  # 3. ZINB Model
  zinb_model <- try(zeroinfl(Y ~ X1 + X2 | X1 + X2, data = sim_data, dist = "negbin"), silent = TRUE)
  if (!inherits(zinb_model, "try-error")) {
    results_matrix[i, 5] <- AIC(zinb_model)
    results_matrix[i, 6] <- BIC(zinb_model)
  }
  
  # 4. DEGPD Model
  fmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1)
  inits_degpd <- c(sqrt(sd(sim_data$Y)/2), sqrt(sd(sim_data$Y)/20), 0.3)
  degpd1 <- try(devgam(fmla_degpd1, data = sim_data, family = "degpd",
                       degpd.args = list(m = 1), trace = 0), silent = TRUE)
  if (!inherits(degpd1, "try-error")) {
    results_matrix[i, 7] <- AIC(degpd1)
    results_matrix[i, 8] <- BIC(degpd1)
  }
  
  # 5. ZIDEGPD Model
  zifmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1, logitpi = ~1)
  inits_zidegpd <- c(sd(sim_data$Y)/2, sd(sim_data$Y)/2, 0.2, mean(sim_data$Y == 0))
  zidegpd1 <- try(devgam(zifmla_degpd1, data = sim_data, family = "zidegpd",
                         degpd.args = list(m = 1), trace = 0, inits = inits_zidegpd), silent = TRUE)
  if (!inherits(zidegpd1, "try-error")) {
    results_matrix[i, 9] <- AIC(zidegpd1)
    results_matrix[i, 10] <- BIC(zidegpd1)
  }
  
  if (i %% 100 == 0) cat("Completed:", i, "\n")
}



# Compute means for each model
mean_aic_bic <- data.frame(
  Model = c("Poisson", "Negative Binomial", "ZINB", "DEGPD", "ZIDEGPD"),
  Mean_AIC = c(
    mean(results_matrix[, 1], na.rm = TRUE),
    mean(results_matrix[, 3], na.rm = TRUE),
    mean(results_matrix[, 5], na.rm = TRUE),
    mean(results_matrix[, 7], na.rm = TRUE),
    mean(results_matrix[, 9], na.rm = TRUE)
  ),
  Mean_BIC = c(
    mean(results_matrix[, 2], na.rm = TRUE),
    mean(results_matrix[, 4], na.rm = TRUE),
    mean(results_matrix[, 6], na.rm = TRUE),
    mean(results_matrix[, 8], na.rm = TRUE),
    mean(results_matrix[, 10], na.rm = TRUE)
  )
)

# Print summary table
print(mean_aic_bic, row.names = FALSE)

# Combine AIC/BIC matrices into a long-format dataframe
aic_bic_df <- data.frame(
  Model = rep(c( "NB", "ZINB", "DEGPD", "ZIDEGPD"), each = nrow(results_matrix)),
  AIC = c(#results_matrix[, 1],
    results_matrix[, 3],
    results_matrix[, 5],
    results_matrix[, 7],
    results_matrix[, 9]),
  BIC = c(#results_matrix[, 2],
    results_matrix[, 4],
    results_matrix[, 6],
    results_matrix[, 8],
    results_matrix[, 10])
)

# Convert to long format
aic_bic_long <- pivot_longer(aic_bic_df, cols = c(AIC, BIC), names_to = "Criterion", values_to = "Value")

# Reorder Model to show DEGPD and ZIDEGPD on the left
aic_bic_long$Model <- factor(aic_bic_long$Model,
                             levels = c("DEGPD", "ZIDEGPD", "NB", "ZINB")
)

# Compute mean of the best model (lowest minimum) per criterion
min_per_model <- aic_bic_long %>%
  group_by(Model, Criterion) %>%
  summarise(Min_Value = min(Value, na.rm = TRUE), .groups = "drop")

lowest_model_per_criterion <- min_per_model %>%
  group_by(Criterion) %>%
  filter(Min_Value == min(Min_Value)) %>%
  select(Criterion, Model)

mean_value_of_best_model <- aic_bic_long %>%
  inner_join(lowest_model_per_criterion, by = c("Criterion", "Model")) %>%
  group_by(Criterion) %>%
  summarise(Line_Value = mean(Value, na.rm = TRUE), .groups = "drop")

# ======= 1. Flask-Shaped Boxplot (varwidth) =======
ggplot(aic_bic_long, aes(x = Model, y = Value, fill = Model)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, varwidth = TRUE) +
  #geom_hline(data = mean_value_of_best_model,
  #aes(yintercept = Line_Value),
  #color = "red", linetype = "dashed", linewidth = 0.5,
  # inherit.aes = FALSE) +
  facet_wrap(~ Criterion, scales = "free_y") +
  labs(title = "Flask-Shaped Boxplots (Varwidth)\nDashed Red Line = Mean of Best Model",
       x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# ======= 2. Violin + Boxplot =======
p<-ggplot(aic_bic_long, aes(x = Model, y = Value, fill = Model)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.4, outlier.shape = TRUE, fill = "white", color = "black", linewidth = 0.2, alpha=0.6) +
  #geom_hline(data = mean_value_of_best_model,
  #           aes(yintercept = Line_Value),
  #          color = "red", linetype = "dashed", linewidth = 0.3,inherit.aes = FALSE) +
  facet_wrap(~ Criterion, scales = "free_y") +
  labs(title = "AIC and BIC plots of the fitted models at 5% high outliers",
       x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

p




# Save to PDF
ggsave("AIC_BIC_0.05_high_0.25_zero_degpd.pdf", plot = p, width = 6, height = 5)










##===================================================================================================================================================================



###. == Low outliers with 10% ratio and 25% zeros

# Initialize storage matrices
results_matrix <- matrix(NA, nrow = 1000, ncol = 10)
colnames(results_matrix) <- c("Poisson_AIC", "Poisson_BIC",
                              "NB_AIC", "NB_BIC",
                              "ZINB_AIC", "ZINB_BIC",
                              "DEGPD_AIC", "DEGPD_BIC",
                              "ZIDEGPD_AIC", "ZIDEGPD_BIC")

set.seed(123)

for (i in 1:1000) {
  # Simulate data once
  sim_data <- simulate_one_rep(outlier_level = "low", zero_inflation = 0.25, outlier_ratio = 0.10)
  
  # 1. Poisson Model
  poisson_model <- try(glm(Y ~ X1 + X2, family = poisson(link = "log"), data = sim_data), silent = TRUE)
  if (!inherits(poisson_model, "try-error")) {
    results_matrix[i, 1] <- AIC(poisson_model)
    results_matrix[i, 2] <- BIC(poisson_model)
  }
  
  # 2. Negative Binomial Model
  nb_model <- try(glm.nb(Y ~ X1 + X2, data = sim_data), silent = TRUE)
  if (!inherits(nb_model, "try-error")) {
    results_matrix[i, 3] <- AIC(nb_model)
    results_matrix[i, 4] <- BIC(nb_model)
  }
  
  # 3. ZINB Model
  zinb_model <- try(zeroinfl(Y ~ X1 + X2 | X1 + X2, data = sim_data, dist = "negbin"), silent = TRUE)
  if (!inherits(zinb_model, "try-error")) {
    results_matrix[i, 5] <- AIC(zinb_model)
    results_matrix[i, 6] <- BIC(zinb_model)
  }
  
  # 4. DEGPD Model
  fmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1)
  inits_degpd <- c(sqrt(sd(sim_data$Y)/2), sqrt(sd(sim_data$Y)/20), 0.3)
  degpd1 <- try(devgam(fmla_degpd1, data = sim_data, family = "degpd",
                       degpd.args = list(m = 1), trace = 0), silent = TRUE)
  if (!inherits(degpd1, "try-error")) {
    results_matrix[i, 7] <- AIC(degpd1)
    results_matrix[i, 8] <- BIC(degpd1)
  }
  
  # 5. ZIDEGPD Model
  zifmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1, logitpi = ~1)
  inits_zidegpd <- c(sd(sim_data$Y)/2, sd(sim_data$Y)/2, 0.2, mean(sim_data$Y == 0))
  zidegpd1 <- try(devgam(zifmla_degpd1, data = sim_data, family = "zidegpd",
                         degpd.args = list(m = 1), trace = 0, inits = inits_zidegpd), silent = TRUE)
  if (!inherits(zidegpd1, "try-error")) {
    results_matrix[i, 9] <- AIC(zidegpd1)
    results_matrix[i, 10] <- BIC(zidegpd1)
  }
  
  if (i %% 100 == 0) cat("Completed:", i, "\n")
}




# Compute means for each model
mean_aic_bic <- data.frame(
  Model = c("Poisson", "Negative Binomial", "ZINB", "DEGPD", "ZIDEGPD"),
  Mean_AIC = c(
    mean(results_matrix[, 1], na.rm = TRUE),
    mean(results_matrix[, 3], na.rm = TRUE),
    mean(results_matrix[, 5], na.rm = TRUE),
    mean(results_matrix[, 7], na.rm = TRUE),
    mean(results_matrix[, 9], na.rm = TRUE)
  ),
  Mean_BIC = c(
    mean(results_matrix[, 2], na.rm = TRUE),
    mean(results_matrix[, 4], na.rm = TRUE),
    mean(results_matrix[, 6], na.rm = TRUE),
    mean(results_matrix[, 8], na.rm = TRUE),
    mean(results_matrix[, 10], na.rm = TRUE)
  )
)

# Print summary table
print(mean_aic_bic, row.names = FALSE)

# Combine AIC/BIC matrices into a long-format dataframe
aic_bic_df <- data.frame(
  Model = rep(c( "NB", "ZINB", "DEGPD", "ZIDEGPD"), each = nrow(results_matrix)),
  AIC = c(#results_matrix[, 1],
    results_matrix[, 3],
    results_matrix[, 5],
    results_matrix[, 7],
    results_matrix[, 9]),
  BIC = c(#results_matrix[, 2],
    results_matrix[, 4],
    results_matrix[, 6],
    results_matrix[, 8],
    results_matrix[, 10])
)

# Convert to long format
aic_bic_long <- pivot_longer(aic_bic_df, cols = c(AIC, BIC), names_to = "Criterion", values_to = "Value")

# Reorder Model to show DEGPD and ZIDEGPD on the left
aic_bic_long$Model <- factor(aic_bic_long$Model,
                             levels = c("DEGPD", "ZIDEGPD", "NB", "ZINB")
)

# Compute mean of the best model (lowest minimum) per criterion
min_per_model <- aic_bic_long %>%
  group_by(Model, Criterion) %>%
  summarise(Min_Value = min(Value, na.rm = TRUE), .groups = "drop")

lowest_model_per_criterion <- min_per_model %>%
  group_by(Criterion) %>%
  filter(Min_Value == min(Min_Value)) %>%
  select(Criterion, Model)

mean_value_of_best_model <- aic_bic_long %>%
  inner_join(lowest_model_per_criterion, by = c("Criterion", "Model")) %>%
  group_by(Criterion) %>%
  summarise(Line_Value = mean(Value, na.rm = TRUE), .groups = "drop")

# ======= 1. Flask-Shaped Boxplot (varwidth) =======
ggplot(aic_bic_long, aes(x = Model, y = Value, fill = Model)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, varwidth = TRUE) +
  #geom_hline(data = mean_value_of_best_model,
  #aes(yintercept = Line_Value),
  #color = "red", linetype = "dashed", linewidth = 0.5,
  # inherit.aes = FALSE) +
  facet_wrap(~ Criterion, scales = "free_y") +
  labs(title = "Flask-Shaped Boxplots (Varwidth)\nDashed Red Line = Mean of Best Model",
       x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# ======= 2. Violin + Boxplot =======
p<-ggplot(aic_bic_long, aes(x = Model, y = Value, fill = Model)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.4, outlier.shape = TRUE, fill = "white", color = "black", linewidth = 0.2, alpha=0.6) +
  #geom_hline(data = mean_value_of_best_model,
  #           aes(yintercept = Line_Value),
  #          color = "red", linetype = "dashed", linewidth = 0.3,inherit.aes = FALSE) +
  facet_wrap(~ Criterion, scales = "free_y") +
  labs(title = "AIC and BIC plots of the fitted models at 10% low outliers",
       x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

p




# Save to PDF
ggsave("AIC_BIC_0.10_low_0.25_zero_degpd.pdf", plot = p, width = 6, height = 5)








###. == medium outliers with 0.10 ratio and 25% zeros

# Initialize storage matrices
results_matrix <- matrix(NA, nrow = 1000, ncol = 10)
colnames(results_matrix) <- c("Poisson_AIC", "Poisson_BIC",
                              "NB_AIC", "NB_BIC",
                              "ZINB_AIC", "ZINB_BIC",
                              "DEGPD_AIC", "DEGPD_BIC",
                              "ZIDEGPD_AIC", "ZIDEGPD_BIC")

set.seed(123)

for (i in 1:1000) {
  # Simulate data once
  sim_data <- simulate_one_rep(outlier_level = "medium", zero_inflation = 0.25, outlier_ratio = 0.10)
  
  # 1. Poisson Model
  poisson_model <- try(glm(Y ~ X1 + X2, family = poisson(link = "log"), data = sim_data), silent = TRUE)
  if (!inherits(poisson_model, "try-error")) {
    results_matrix[i, 1] <- AIC(poisson_model)
    results_matrix[i, 2] <- BIC(poisson_model)
  }
  
  # 2. Negative Binomial Model
  nb_model <- try(glm.nb(Y ~ X1 + X2, data = sim_data), silent = TRUE)
  if (!inherits(nb_model, "try-error")) {
    results_matrix[i, 3] <- AIC(nb_model)
    results_matrix[i, 4] <- BIC(nb_model)
  }
  
  # 3. ZINB Model
  zinb_model <- try(zeroinfl(Y ~ X1 + X2 | X1 + X2, data = sim_data, dist = "negbin"), silent = TRUE)
  if (!inherits(zinb_model, "try-error")) {
    results_matrix[i, 5] <- AIC(zinb_model)
    results_matrix[i, 6] <- BIC(zinb_model)
  }
  
  # 4. DEGPD Model
  fmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1)
  inits_degpd <- c(sqrt(sd(sim_data$Y)/2), sqrt(sd(sim_data$Y)/20), 0.3)
  degpd1 <- try(devgam(fmla_degpd1, data = sim_data, family = "degpd",
                       degpd.args = list(m = 1), trace = 0), silent = TRUE)
  if (!inherits(degpd1, "try-error")) {
    results_matrix[i, 7] <- AIC(degpd1)
    results_matrix[i, 8] <- BIC(degpd1)
  }
  
  # 5. ZIDEGPD Model
  zifmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1, logitpi = ~1)
  inits_zidegpd <- c(sd(sim_data$Y)/2, sd(sim_data$Y)/2, 0.2, mean(sim_data$Y == 0))
  zidegpd1 <- try(devgam(zifmla_degpd1, data = sim_data, family = "zidegpd",
                         degpd.args = list(m = 1), trace = 0, inits = inits_zidegpd), silent = TRUE)
  if (!inherits(zidegpd1, "try-error")) {
    results_matrix[i, 9] <- AIC(zidegpd1)
    results_matrix[i, 10] <- BIC(zidegpd1)
  }
  
  if (i %% 100 == 0) cat("Completed:", i, "\n")
}






# Compute means for each model
mean_aic_bic <- data.frame(
  Model = c("Poisson", "Negative Binomial", "ZINB", "DEGPD", "ZIDEGPD"),
  Mean_AIC = c(
    mean(results_matrix[, 1], na.rm = TRUE),
    mean(results_matrix[, 3], na.rm = TRUE),
    mean(results_matrix[, 5], na.rm = TRUE),
    mean(results_matrix[, 7], na.rm = TRUE),
    mean(results_matrix[, 9], na.rm = TRUE)
  ),
  Mean_BIC = c(
    mean(results_matrix[, 2], na.rm = TRUE),
    mean(results_matrix[, 4], na.rm = TRUE),
    mean(results_matrix[, 6], na.rm = TRUE),
    mean(results_matrix[, 8], na.rm = TRUE),
    mean(results_matrix[, 10], na.rm = TRUE)
  )
)

# Print summary table
print(mean_aic_bic, row.names = FALSE)

# Combine AIC/BIC matrices into a long-format dataframe
aic_bic_df <- data.frame(
  Model = rep(c( "NB", "ZINB", "DEGPD", "ZIDEGPD"), each = nrow(results_matrix)),
  AIC = c(#results_matrix[, 1],
    results_matrix[, 3],
    results_matrix[, 5],
    results_matrix[, 7],
    results_matrix[, 9]),
  BIC = c(#results_matrix[, 2],
    results_matrix[, 4],
    results_matrix[, 6],
    results_matrix[, 8],
    results_matrix[, 10])
)

# Convert to long format
aic_bic_long <- pivot_longer(aic_bic_df, cols = c(AIC, BIC), names_to = "Criterion", values_to = "Value")

# Reorder Model to show DEGPD and ZIDEGPD on the left
aic_bic_long$Model <- factor(aic_bic_long$Model,
                             levels = c("DEGPD", "ZIDEGPD", "NB", "ZINB")
)

# Compute mean of the best model (lowest minimum) per criterion
min_per_model <- aic_bic_long %>%
  group_by(Model, Criterion) %>%
  summarise(Min_Value = min(Value, na.rm = TRUE), .groups = "drop")

lowest_model_per_criterion <- min_per_model %>%
  group_by(Criterion) %>%
  filter(Min_Value == min(Min_Value)) %>%
  select(Criterion, Model)

mean_value_of_best_model <- aic_bic_long %>%
  inner_join(lowest_model_per_criterion, by = c("Criterion", "Model")) %>%
  group_by(Criterion) %>%
  summarise(Line_Value = mean(Value, na.rm = TRUE), .groups = "drop")

# ======= 1. Flask-Shaped Boxplot (varwidth) =======
ggplot(aic_bic_long, aes(x = Model, y = Value, fill = Model)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, varwidth = TRUE) +
  #geom_hline(data = mean_value_of_best_model,
  #aes(yintercept = Line_Value),
  #color = "red", linetype = "dashed", linewidth = 0.5,
  # inherit.aes = FALSE) +
  facet_wrap(~ Criterion, scales = "free_y") +
  labs(title = "Flask-Shaped Boxplots (Varwidth)\nDashed Red Line = Mean of Best Model",
       x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# ======= 2. Violin + Boxplot =======
p<-ggplot(aic_bic_long, aes(x = Model, y = Value, fill = Model)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.4, outlier.shape = TRUE, fill = "white", color = "black", linewidth = 0.2, alpha=0.6) +
  #geom_hline(data = mean_value_of_best_model,
  #           aes(yintercept = Line_Value),
  #          color = "red", linetype = "dashed", linewidth = 0.3,inherit.aes = FALSE) +
  facet_wrap(~ Criterion, scales = "free_y") +
  labs(title = "AIC and BIC plots of the fitted models at 10% medium outliers",
       x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

p




# Save to PDF
ggsave("AIC_BIC_0.10_medium_0.25_zero_degpd.pdf", plot = p, width = 6, height = 5)




###. == high outliers with 0.10 ratio and 25% zeros. ===

# Initialize storage matrices
results_matrix <- matrix(NA, nrow = 1000, ncol = 10)
colnames(results_matrix) <- c("Poisson_AIC", "Poisson_BIC",
                              "NB_AIC", "NB_BIC",
                              "ZINB_AIC", "ZINB_BIC",
                              "DEGPD_AIC", "DEGPD_BIC",
                              "ZIDEGPD_AIC", "ZIDEGPD_BIC")

set.seed(123)

for (i in 1:1000) {
  # Simulate data once
  sim_data <- simulate_one_rep(outlier_level = "high", zero_inflation = 0.25, outlier_ratio = 0.10)
  
  # 1. Poisson Model
  poisson_model <- try(glm(Y ~ X1 + X2, family = poisson(link = "log"), data = sim_data), silent = TRUE)
  if (!inherits(poisson_model, "try-error")) {
    results_matrix[i, 1] <- AIC(poisson_model)
    results_matrix[i, 2] <- BIC(poisson_model)
  }
  
  # 2. Negative Binomial Model
  nb_model <- try(glm.nb(Y ~ X1 + X2, data = sim_data), silent = TRUE)
  if (!inherits(nb_model, "try-error")) {
    results_matrix[i, 3] <- AIC(nb_model)
    results_matrix[i, 4] <- BIC(nb_model)
  }
  
  # 3. ZINB Model
  zinb_model <- try(zeroinfl(Y ~ X1 + X2 | X1 + X2, data = sim_data, dist = "negbin"), silent = TRUE)
  if (!inherits(zinb_model, "try-error")) {
    results_matrix[i, 5] <- AIC(zinb_model)
    results_matrix[i, 6] <- BIC(zinb_model)
  }
  
  # 4. DEGPD Model
  fmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1)
  inits_degpd <- c(sqrt(sd(sim_data$Y)/2), sqrt(sd(sim_data$Y)/20), 0.3)
  degpd1 <- try(devgam(fmla_degpd1, data = sim_data, family = "degpd",
                       degpd.args = list(m = 1), trace = 0), silent = TRUE)
  if (!inherits(degpd1, "try-error")) {
    results_matrix[i, 7] <- AIC(degpd1)
    results_matrix[i, 8] <- BIC(degpd1)
  }
  
  # 5. ZIDEGPD Model
  zifmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1, logitpi = ~1)
  inits_zidegpd <- c(sd(sim_data$Y)/2, sd(sim_data$Y)/2, 0.2, mean(sim_data$Y == 0))
  zidegpd1 <- try(devgam(zifmla_degpd1, data = sim_data, family = "zidegpd",
                         degpd.args = list(m = 1), trace = 0, inits = inits_zidegpd), silent = TRUE)
  if (!inherits(zidegpd1, "try-error")) {
    results_matrix[i, 9] <- AIC(zidegpd1)
    results_matrix[i, 10] <- BIC(zidegpd1)
  }
  
  if (i %% 100 == 0) cat("Completed:", i, "\n")
}



# Compute means for each model
mean_aic_bic <- data.frame(
  Model = c("Poisson", "Negative Binomial", "ZINB", "DEGPD", "ZIDEGPD"),
  Mean_AIC = c(
    mean(results_matrix[, 1], na.rm = TRUE),
    mean(results_matrix[, 3], na.rm = TRUE),
    mean(results_matrix[, 5], na.rm = TRUE),
    mean(results_matrix[, 7], na.rm = TRUE),
    mean(results_matrix[, 9], na.rm = TRUE)
  ),
  Mean_BIC = c(
    mean(results_matrix[, 2], na.rm = TRUE),
    mean(results_matrix[, 4], na.rm = TRUE),
    mean(results_matrix[, 6], na.rm = TRUE),
    mean(results_matrix[, 8], na.rm = TRUE),
    mean(results_matrix[, 10], na.rm = TRUE)
  )
)

# Print summary table
print(mean_aic_bic, row.names = FALSE)

# Combine AIC/BIC matrices into a long-format dataframe
aic_bic_df <- data.frame(
  Model = rep(c( "NB", "ZINB", "DEGPD", "ZIDEGPD"), each = nrow(results_matrix)),
  AIC = c(#results_matrix[, 1],
    results_matrix[, 3],
    results_matrix[, 5],
    results_matrix[, 7],
    results_matrix[, 9]),
  BIC = c(#results_matrix[, 2],
    results_matrix[, 4],
    results_matrix[, 6],
    results_matrix[, 8],
    results_matrix[, 10])
)

# Convert to long format
aic_bic_long <- pivot_longer(aic_bic_df, cols = c(AIC, BIC), names_to = "Criterion", values_to = "Value")

# Reorder Model to show DEGPD and ZIDEGPD on the left
aic_bic_long$Model <- factor(aic_bic_long$Model,
                             levels = c("DEGPD", "ZIDEGPD", "NB", "ZINB")
)

# Compute mean of the best model (lowest minimum) per criterion
min_per_model <- aic_bic_long %>%
  group_by(Model, Criterion) %>%
  summarise(Min_Value = min(Value, na.rm = TRUE), .groups = "drop")

lowest_model_per_criterion <- min_per_model %>%
  group_by(Criterion) %>%
  filter(Min_Value == min(Min_Value)) %>%
  select(Criterion, Model)

mean_value_of_best_model <- aic_bic_long %>%
  inner_join(lowest_model_per_criterion, by = c("Criterion", "Model")) %>%
  group_by(Criterion) %>%
  summarise(Line_Value = mean(Value, na.rm = TRUE), .groups = "drop")

# ======= 1. Flask-Shaped Boxplot (varwidth) =======
ggplot(aic_bic_long, aes(x = Model, y = Value, fill = Model)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, varwidth = TRUE) +
  #geom_hline(data = mean_value_of_best_model,
  #aes(yintercept = Line_Value),
  #color = "red", linetype = "dashed", linewidth = 0.5,
  # inherit.aes = FALSE) +
  facet_wrap(~ Criterion, scales = "free_y") +
  labs(title = "Flask-Shaped Boxplots (Varwidth)\nDashed Red Line = Mean of Best Model",
       x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# ======= 2. Violin + Boxplot =======
p<-ggplot(aic_bic_long, aes(x = Model, y = Value, fill = Model)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.4, outlier.shape = TRUE, fill = "white", color = "black", linewidth = 0.2, alpha=0.6) +
  #geom_hline(data = mean_value_of_best_model,
  #           aes(yintercept = Line_Value),
  #          color = "red", linetype = "dashed", linewidth = 0.3,inherit.aes = FALSE) +
  facet_wrap(~ Criterion, scales = "free_y") +
  labs(title = "AIC and BIC plots of the fitted models at 10% high outliers",
       x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

p




# Save to PDF
ggsave("AIC_BIC_0.10_high_0.25_zero_degpd.pdf", plot = p, width = 6, height = 5)















##. 50 % zeros========================================================================================================================





###. == Low outliers with 0.01 ratio and 25% zeros

# Initialize storage matrices
results_matrix <- matrix(NA, nrow = 1000, ncol = 10)
colnames(results_matrix) <- c("Poisson_AIC", "Poisson_BIC",
                              "NB_AIC", "NB_BIC",
                              "ZINB_AIC", "ZINB_BIC",
                              "DEGPD_AIC", "DEGPD_BIC",
                              "ZIDEGPD_AIC", "ZIDEGPD_BIC")

set.seed(123)

for (i in 1:1000) {
  # Simulate data once
  sim_data <- simulate_one_rep(outlier_level = "low", zero_inflation = 0.50, outlier_ratio = 0.01)
  
  # 1. Poisson Model
  poisson_model <- try(glm(Y ~ X1 + X2, family = poisson(link = "log"), data = sim_data), silent = TRUE)
  if (!inherits(poisson_model, "try-error")) {
    results_matrix[i, 1] <- AIC(poisson_model)
    results_matrix[i, 2] <- BIC(poisson_model)
  }
  
  # 2. Negative Binomial Model
  nb_model <- try(glm.nb(Y ~ X1 + X2, data = sim_data), silent = TRUE)
  if (!inherits(nb_model, "try-error")) {
    results_matrix[i, 3] <- AIC(nb_model)
    results_matrix[i, 4] <- BIC(nb_model)
  }
  
  # 3. ZINB Model
  zinb_model <- try(zeroinfl(Y ~ X1 + X2 | X1 + X2, data = sim_data, dist = "negbin"), silent = TRUE)
  if (!inherits(zinb_model, "try-error")) {
    results_matrix[i, 5] <- AIC(zinb_model)
    results_matrix[i, 6] <- BIC(zinb_model)
  }
  
  # 4. DEGPD Model
  fmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1)
  inits_degpd <- c(sqrt(sd(sim_data$Y)/2), sqrt(sd(sim_data$Y)/20), 0.3)
  degpd1 <- try(devgam(fmla_degpd1, data = sim_data, family = "degpd",
                       degpd.args = list(m = 1), trace = 0), silent = TRUE)
  if (!inherits(degpd1, "try-error")) {
    results_matrix[i, 7] <- AIC(degpd1)
    results_matrix[i, 8] <- BIC(degpd1)
  }
  
  # 5. ZIDEGPD Model
  zifmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1, logitpi = ~1)
  inits_zidegpd <- c(sd(sim_data$Y)/2, sd(sim_data$Y)/2, 0.2, mean(sim_data$Y == 0))
  zidegpd1 <- try(devgam(zifmla_degpd1, data = sim_data, family = "zidegpd",
                         degpd.args = list(m = 1), trace = 0, inits = inits_zidegpd), silent = TRUE)
  if (!inherits(zidegpd1, "try-error")) {
    results_matrix[i, 9] <- AIC(zidegpd1)
    results_matrix[i, 10] <- BIC(zidegpd1)
  }
  
  if (i %% 100 == 0) cat("Completed:", i, "\n")
}






# Compute means for each model
mean_aic_bic <- data.frame(
  Model = c("Poisson", "Negative Binomial", "ZINB", "DEGPD", "ZIDEGPD"),
  Mean_AIC = c(
    mean(results_matrix[, 1], na.rm = TRUE),
    mean(results_matrix[, 3], na.rm = TRUE),
    mean(results_matrix[, 5], na.rm = TRUE),
    mean(results_matrix[, 7], na.rm = TRUE),
    mean(results_matrix[, 9], na.rm = TRUE)
  ),
  Mean_BIC = c(
    mean(results_matrix[, 2], na.rm = TRUE),
    mean(results_matrix[, 4], na.rm = TRUE),
    mean(results_matrix[, 6], na.rm = TRUE),
    mean(results_matrix[, 8], na.rm = TRUE),
    mean(results_matrix[, 10], na.rm = TRUE)
  )
)

# Print summary table
print(mean_aic_bic, row.names = FALSE)

# Combine AIC/BIC matrices into a long-format dataframe
aic_bic_df <- data.frame(
  Model = rep(c( "NB", "ZINB", "DEGPD", "ZIDEGPD"), each = nrow(results_matrix)),
  AIC = c(#results_matrix[, 1],
    results_matrix[, 3],
    results_matrix[, 5],
    results_matrix[, 7],
    results_matrix[, 9]),
  BIC = c(#results_matrix[, 2],
    results_matrix[, 4],
    results_matrix[, 6],
    results_matrix[, 8],
    results_matrix[, 10])
)

# Convert to long format
aic_bic_long <- pivot_longer(aic_bic_df, cols = c(AIC, BIC), names_to = "Criterion", values_to = "Value")

# Reorder Model to show DEGPD and ZIDEGPD on the left
aic_bic_long$Model <- factor(aic_bic_long$Model,
                             levels = c("DEGPD", "ZIDEGPD", "NB", "ZINB")
)

# Compute mean of the best model (lowest minimum) per criterion
min_per_model <- aic_bic_long %>%
  group_by(Model, Criterion) %>%
  summarise(Min_Value = min(Value, na.rm = TRUE), .groups = "drop")

lowest_model_per_criterion <- min_per_model %>%
  group_by(Criterion) %>%
  filter(Min_Value == min(Min_Value)) %>%
  select(Criterion, Model)

mean_value_of_best_model <- aic_bic_long %>%
  inner_join(lowest_model_per_criterion, by = c("Criterion", "Model")) %>%
  group_by(Criterion) %>%
  summarise(Line_Value = mean(Value, na.rm = TRUE), .groups = "drop")

# ======= 1. Flask-Shaped Boxplot (varwidth) =======
ggplot(aic_bic_long, aes(x = Model, y = Value, fill = Model)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, varwidth = TRUE) +
  #geom_hline(data = mean_value_of_best_model,
  #aes(yintercept = Line_Value),
  #color = "red", linetype = "dashed", linewidth = 0.5,
  # inherit.aes = FALSE) +
  facet_wrap(~ Criterion, scales = "free_y") +
  labs(title = "Flask-Shaped Boxplots (Varwidth)\nDashed Red Line = Mean of Best Model",
       x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# ======= 2. Violin + Boxplot =======
p<-ggplot(aic_bic_long, aes(x = Model, y = Value, fill = Model)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.4, outlier.shape = TRUE, fill = "white", color = "black", linewidth = 0.2, alpha=0.6) +
  #geom_hline(data = mean_value_of_best_model,
  #           aes(yintercept = Line_Value),
  #          color = "red", linetype = "dashed", linewidth = 0.3,inherit.aes = FALSE) +
  facet_wrap(~ Criterion, scales = "free_y") +
  labs(title = "AIC and BIC plots of the fitted models at 1% low outliers",
       x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

p




# Save to PDF
ggsave("AIC_BIC_0.01_low_0.50_zero_degpd.pdf", plot = p, width = 6, height = 5)








###. == medium outliers with 0.01 ratio and 25% zeros

# Initialize storage matrices
results_matrix <- matrix(NA, nrow = 1000, ncol = 10)
colnames(results_matrix) <- c("Poisson_AIC", "Poisson_BIC",
                              "NB_AIC", "NB_BIC",
                              "ZINB_AIC", "ZINB_BIC",
                              "DEGPD_AIC", "DEGPD_BIC",
                              "ZIDEGPD_AIC", "ZIDEGPD_BIC")

set.seed(123)

for (i in 1:1000) {
  # Simulate data once
  sim_data <- simulate_one_rep(outlier_level = "medium", zero_inflation = 0.50, outlier_ratio = 0.01)
  
  # 1. Poisson Model
  poisson_model <- try(glm(Y ~ X1 + X2, family = poisson(link = "log"), data = sim_data), silent = TRUE)
  if (!inherits(poisson_model, "try-error")) {
    results_matrix[i, 1] <- AIC(poisson_model)
    results_matrix[i, 2] <- BIC(poisson_model)
  }
  
  # 2. Negative Binomial Model
  nb_model <- try(glm.nb(Y ~ X1 + X2, data = sim_data), silent = TRUE)
  if (!inherits(nb_model, "try-error")) {
    results_matrix[i, 3] <- AIC(nb_model)
    results_matrix[i, 4] <- BIC(nb_model)
  }
  
  # 3. ZINB Model
  zinb_model <- try(zeroinfl(Y ~ X1 + X2 | X1 + X2, data = sim_data, dist = "negbin"), silent = TRUE)
  if (!inherits(zinb_model, "try-error")) {
    results_matrix[i, 5] <- AIC(zinb_model)
    results_matrix[i, 6] <- BIC(zinb_model)
  }
  
  # 4. DEGPD Model
  fmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1)
  inits_degpd <- c(sqrt(sd(sim_data$Y)/2), sqrt(sd(sim_data$Y)/20), 0.3)
  degpd1 <- try(devgam(fmla_degpd1, data = sim_data, family = "degpd",
                       degpd.args = list(m = 1), trace = 0), silent = TRUE)
  if (!inherits(degpd1, "try-error")) {
    results_matrix[i, 7] <- AIC(degpd1)
    results_matrix[i, 8] <- BIC(degpd1)
  }
  
  # 5. ZIDEGPD Model
  zifmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1, logitpi = ~1)
  inits_zidegpd <- c(sd(sim_data$Y)/2, sd(sim_data$Y)/2, 0.2, mean(sim_data$Y == 0))
  zidegpd1 <- try(devgam(zifmla_degpd1, data = sim_data, family = "zidegpd",
                         degpd.args = list(m = 1), trace = 0, inits = inits_zidegpd), silent = TRUE)
  if (!inherits(zidegpd1, "try-error")) {
    results_matrix[i, 9] <- AIC(zidegpd1)
    results_matrix[i, 10] <- BIC(zidegpd1)
  }
  
  if (i %% 100 == 0) cat("Completed:", i, "\n")
}






# Compute means for each model
mean_aic_bic <- data.frame(
  Model = c("Poisson", "Negative Binomial", "ZINB", "DEGPD", "ZIDEGPD"),
  Mean_AIC = c(
    mean(results_matrix[, 1], na.rm = TRUE),
    mean(results_matrix[, 3], na.rm = TRUE),
    mean(results_matrix[, 5], na.rm = TRUE),
    mean(results_matrix[, 7], na.rm = TRUE),
    mean(results_matrix[, 9], na.rm = TRUE)
  ),
  Mean_BIC = c(
    mean(results_matrix[, 2], na.rm = TRUE),
    mean(results_matrix[, 4], na.rm = TRUE),
    mean(results_matrix[, 6], na.rm = TRUE),
    mean(results_matrix[, 8], na.rm = TRUE),
    mean(results_matrix[, 10], na.rm = TRUE)
  )
)

# Print summary table
print(mean_aic_bic, row.names = FALSE)

# Combine AIC/BIC matrices into a long-format dataframe
aic_bic_df <- data.frame(
  Model = rep(c( "NB", "ZINB", "DEGPD", "ZIDEGPD"), each = nrow(results_matrix)),
  AIC = c(#results_matrix[, 1],
    results_matrix[, 3],
    results_matrix[, 5],
    results_matrix[, 7],
    results_matrix[, 9]),
  BIC = c(#results_matrix[, 2],
    results_matrix[, 4],
    results_matrix[, 6],
    results_matrix[, 8],
    results_matrix[, 10])
)

# Convert to long format
aic_bic_long <- pivot_longer(aic_bic_df, cols = c(AIC, BIC), names_to = "Criterion", values_to = "Value")

# Reorder Model to show DEGPD and ZIDEGPD on the left
aic_bic_long$Model <- factor(aic_bic_long$Model,
                             levels = c("DEGPD", "ZIDEGPD", "NB", "ZINB")
)

# Compute mean of the best model (lowest minimum) per criterion
min_per_model <- aic_bic_long %>%
  group_by(Model, Criterion) %>%
  summarise(Min_Value = min(Value, na.rm = TRUE), .groups = "drop")

lowest_model_per_criterion <- min_per_model %>%
  group_by(Criterion) %>%
  filter(Min_Value == min(Min_Value)) %>%
  select(Criterion, Model)

mean_value_of_best_model <- aic_bic_long %>%
  inner_join(lowest_model_per_criterion, by = c("Criterion", "Model")) %>%
  group_by(Criterion) %>%
  summarise(Line_Value = mean(Value, na.rm = TRUE), .groups = "drop")

# ======= 1. Flask-Shaped Boxplot (varwidth) =======
ggplot(aic_bic_long, aes(x = Model, y = Value, fill = Model)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, varwidth = TRUE) +
  #geom_hline(data = mean_value_of_best_model,
  #aes(yintercept = Line_Value),
  #color = "red", linetype = "dashed", linewidth = 0.5,
  # inherit.aes = FALSE) +
  facet_wrap(~ Criterion, scales = "free_y") +
  labs(title = "Flask-Shaped Boxplots (Varwidth)\nDashed Red Line = Mean of Best Model",
       x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# ======= 2. Violin + Boxplot =======
p<-ggplot(aic_bic_long, aes(x = Model, y = Value, fill = Model)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.4, outlier.shape = TRUE, fill = "white", color = "black", linewidth = 0.2, alpha=0.6) +
  #geom_hline(data = mean_value_of_best_model,
  #           aes(yintercept = Line_Value),
  #          color = "red", linetype = "dashed", linewidth = 0.3,inherit.aes = FALSE) +
  facet_wrap(~ Criterion, scales = "free_y") +
  labs(title = "AIC and BIC plots of the fitted models at 1% medium outliers",
       x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

p




# Save to PDF
ggsave("AIC_BIC_0.01_medium_0.50_zero_degpd.pdf", plot = p, width = 6, height = 5)




###. == high outliers with 0.01 ratio and 25% zeros. ===

# Initialize storage matrices
results_matrix <- matrix(NA, nrow = 1000, ncol = 10)
colnames(results_matrix) <- c("Poisson_AIC", "Poisson_BIC",
                              "NB_AIC", "NB_BIC",
                              "ZINB_AIC", "ZINB_BIC",
                              "DEGPD_AIC", "DEGPD_BIC",
                              "ZIDEGPD_AIC", "ZIDEGPD_BIC")

set.seed(123)

for (i in 1:1000) {
  # Simulate data once
  sim_data <- simulate_one_rep(outlier_level = "high", zero_inflation = 0.50, outlier_ratio = 0.01)
  
  # 1. Poisson Model
  poisson_model <- try(glm(Y ~ X1 + X2, family = poisson(link = "log"), data = sim_data), silent = TRUE)
  if (!inherits(poisson_model, "try-error")) {
    results_matrix[i, 1] <- AIC(poisson_model)
    results_matrix[i, 2] <- BIC(poisson_model)
  }
  
  # 2. Negative Binomial Model
  nb_model <- try(glm.nb(Y ~ X1 + X2, data = sim_data), silent = TRUE)
  if (!inherits(nb_model, "try-error")) {
    results_matrix[i, 3] <- AIC(nb_model)
    results_matrix[i, 4] <- BIC(nb_model)
  }
  
  # 3. ZINB Model
  zinb_model <- try(zeroinfl(Y ~ X1 + X2 | X1 + X2, data = sim_data, dist = "negbin"), silent = TRUE)
  if (!inherits(zinb_model, "try-error")) {
    results_matrix[i, 5] <- AIC(zinb_model)
    results_matrix[i, 6] <- BIC(zinb_model)
  }
  
  # 4. DEGPD Model
  fmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1)
  inits_degpd <- c(sqrt(sd(sim_data$Y)/2), sqrt(sd(sim_data$Y)/20), 0.3)
  degpd1 <- try(devgam(fmla_degpd1, data = sim_data, family = "degpd",
                       degpd.args = list(m = 1), trace = 0), silent = TRUE)
  if (!inherits(degpd1, "try-error")) {
    results_matrix[i, 7] <- AIC(degpd1)
    results_matrix[i, 8] <- BIC(degpd1)
  }
  
  # 5. ZIDEGPD Model
  zifmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1, logitpi = ~1)
  inits_zidegpd <- c(sd(sim_data$Y)/2, sd(sim_data$Y)/2, 0.2, mean(sim_data$Y == 0))
  zidegpd1 <- try(devgam(zifmla_degpd1, data = sim_data, family = "zidegpd",
                         degpd.args = list(m = 1), trace = 0, inits = inits_zidegpd), silent = TRUE)
  if (!inherits(zidegpd1, "try-error")) {
    results_matrix[i, 9] <- AIC(zidegpd1)
    results_matrix[i, 10] <- BIC(zidegpd1)
  }
  
  if (i %% 100 == 0) cat("Completed:", i, "\n")
}



# Compute means for each model
mean_aic_bic <- data.frame(
  Model = c("Poisson", "Negative Binomial", "ZINB", "DEGPD", "ZIDEGPD"),
  Mean_AIC = c(
    mean(results_matrix[, 1], na.rm = TRUE),
    mean(results_matrix[, 3], na.rm = TRUE),
    mean(results_matrix[, 5], na.rm = TRUE),
    mean(results_matrix[, 7], na.rm = TRUE),
    mean(results_matrix[, 9], na.rm = TRUE)
  ),
  Mean_BIC = c(
    mean(results_matrix[, 2], na.rm = TRUE),
    mean(results_matrix[, 4], na.rm = TRUE),
    mean(results_matrix[, 6], na.rm = TRUE),
    mean(results_matrix[, 8], na.rm = TRUE),
    mean(results_matrix[, 10], na.rm = TRUE)
  )
)

# Print summary table
print(mean_aic_bic, row.names = FALSE)

# Combine AIC/BIC matrices into a long-format dataframe
aic_bic_df <- data.frame(
  Model = rep(c( "NB", "ZINB", "DEGPD", "ZIDEGPD"), each = nrow(results_matrix)),
  AIC = c(#results_matrix[, 1],
    results_matrix[, 3],
    results_matrix[, 5],
    results_matrix[, 7],
    results_matrix[, 9]),
  BIC = c(#results_matrix[, 2],
    results_matrix[, 4],
    results_matrix[, 6],
    results_matrix[, 8],
    results_matrix[, 10])
)

# Convert to long format
aic_bic_long <- pivot_longer(aic_bic_df, cols = c(AIC, BIC), names_to = "Criterion", values_to = "Value")

# Reorder Model to show DEGPD and ZIDEGPD on the left
aic_bic_long$Model <- factor(aic_bic_long$Model,
                             levels = c("DEGPD", "ZIDEGPD", "NB", "ZINB")
)

# Compute mean of the best model (lowest minimum) per criterion
min_per_model <- aic_bic_long %>%
  group_by(Model, Criterion) %>%
  summarise(Min_Value = min(Value, na.rm = TRUE), .groups = "drop")

lowest_model_per_criterion <- min_per_model %>%
  group_by(Criterion) %>%
  filter(Min_Value == min(Min_Value)) %>%
  select(Criterion, Model)

mean_value_of_best_model <- aic_bic_long %>%
  inner_join(lowest_model_per_criterion, by = c("Criterion", "Model")) %>%
  group_by(Criterion) %>%
  summarise(Line_Value = mean(Value, na.rm = TRUE), .groups = "drop")

# ======= 1. Flask-Shaped Boxplot (varwidth) =======
ggplot(aic_bic_long, aes(x = Model, y = Value, fill = Model)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, varwidth = TRUE) +
  #geom_hline(data = mean_value_of_best_model,
  #aes(yintercept = Line_Value),
  #color = "red", linetype = "dashed", linewidth = 0.5,
  # inherit.aes = FALSE) +
  facet_wrap(~ Criterion, scales = "free_y") +
  labs(title = "Flask-Shaped Boxplots (Varwidth)\nDashed Red Line = Mean of Best Model",
       x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# ======= 2. Violin + Boxplot =======
p<-ggplot(aic_bic_long, aes(x = Model, y = Value, fill = Model)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.4, outlier.shape = TRUE, fill = "white", color = "black", linewidth = 0.2, alpha=0.6) +
  #geom_hline(data = mean_value_of_best_model,
  #           aes(yintercept = Line_Value),
  #          color = "red", linetype = "dashed", linewidth = 0.3,inherit.aes = FALSE) +
  facet_wrap(~ Criterion, scales = "free_y") +
  labs(title = "AIC and BIC plots of the fitted models at 1% high outliers",
       x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

p




# Save to PDF
ggsave("AIC_BIC_0.01_high_0.50_zero_degpd.pdf", plot = p, width = 6, height = 5)






##==============================================================================================================================================






###. == Low outliers with 5% ratio and 25% zeros

# Initialize storage matrices
results_matrix <- matrix(NA, nrow = 1000, ncol = 10)
colnames(results_matrix) <- c("Poisson_AIC", "Poisson_BIC",
                              "NB_AIC", "NB_BIC",
                              "ZINB_AIC", "ZINB_BIC",
                              "DEGPD_AIC", "DEGPD_BIC",
                              "ZIDEGPD_AIC", "ZIDEGPD_BIC")

set.seed(123)

for (i in 1:1000) {
  # Simulate data once
  sim_data <- simulate_one_rep(outlier_level = "low", zero_inflation = 0.50, outlier_ratio = 0.05)
  
  # 1. Poisson Model
  poisson_model <- try(glm(Y ~ X1 + X2, family = poisson(link = "log"), data = sim_data), silent = TRUE)
  if (!inherits(poisson_model, "try-error")) {
    results_matrix[i, 1] <- AIC(poisson_model)
    results_matrix[i, 2] <- BIC(poisson_model)
  }
  
  # 2. Negative Binomial Model
  nb_model <- try(glm.nb(Y ~ X1 + X2, data = sim_data), silent = TRUE)
  if (!inherits(nb_model, "try-error")) {
    results_matrix[i, 3] <- AIC(nb_model)
    results_matrix[i, 4] <- BIC(nb_model)
  }
  
  # 3. ZINB Model
  zinb_model <- try(zeroinfl(Y ~ X1 + X2 | X1 + X2, data = sim_data, dist = "negbin"), silent = TRUE)
  if (!inherits(zinb_model, "try-error")) {
    results_matrix[i, 5] <- AIC(zinb_model)
    results_matrix[i, 6] <- BIC(zinb_model)
  }
  
  # 4. DEGPD Model
  fmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1)
  inits_degpd <- c(sqrt(sd(sim_data$Y)/2), sqrt(sd(sim_data$Y)/20), 0.3)
  degpd1 <- try(devgam(fmla_degpd1, data = sim_data, family = "degpd",
                       degpd.args = list(m = 1), trace = 0), silent = TRUE)
  if (!inherits(degpd1, "try-error")) {
    results_matrix[i, 7] <- AIC(degpd1)
    results_matrix[i, 8] <- BIC(degpd1)
  }
  
  # 5. ZIDEGPD Model
  zifmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1, logitpi = ~1)
  inits_zidegpd <- c(sd(sim_data$Y)/2, sd(sim_data$Y)/2, 0.2, mean(sim_data$Y == 0))
  zidegpd1 <- try(devgam(zifmla_degpd1, data = sim_data, family = "zidegpd",
                         degpd.args = list(m = 1), trace = 0, inits = inits_zidegpd), silent = TRUE)
  if (!inherits(zidegpd1, "try-error")) {
    results_matrix[i, 9] <- AIC(zidegpd1)
    results_matrix[i, 10] <- BIC(zidegpd1)
  }
  
  if (i %% 100 == 0) cat("Completed:", i, "\n")
}







# Load required packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggdist)
library(ggbeeswarm)

# Compute means for each model
mean_aic_bic <- data.frame(
  Model = c("Poisson", "Negative Binomial", "ZINB", "DEGPD", "ZIDEGPD"),
  Mean_AIC = c(
    mean(results_matrix[, 1], na.rm = TRUE),
    mean(results_matrix[, 3], na.rm = TRUE),
    mean(results_matrix[, 5], na.rm = TRUE),
    mean(results_matrix[, 7], na.rm = TRUE),
    mean(results_matrix[, 9], na.rm = TRUE)
  ),
  Mean_BIC = c(
    mean(results_matrix[, 2], na.rm = TRUE),
    mean(results_matrix[, 4], na.rm = TRUE),
    mean(results_matrix[, 6], na.rm = TRUE),
    mean(results_matrix[, 8], na.rm = TRUE),
    mean(results_matrix[, 10], na.rm = TRUE)
  )
)

# Print summary table
print(mean_aic_bic, row.names = FALSE)

# Combine AIC/BIC matrices into a long-format dataframe
aic_bic_df <- data.frame(
  Model = rep(c( "NB", "ZINB", "DEGPD", "ZIDEGPD"), each = nrow(results_matrix)),
  AIC = c(#results_matrix[, 1],
    results_matrix[, 3],
    results_matrix[, 5],
    results_matrix[, 7],
    results_matrix[, 9]),
  BIC = c(#results_matrix[, 2],
    results_matrix[, 4],
    results_matrix[, 6],
    results_matrix[, 8],
    results_matrix[, 10])
)

# Convert to long format
aic_bic_long <- pivot_longer(aic_bic_df, cols = c(AIC, BIC), names_to = "Criterion", values_to = "Value")

# Reorder Model to show DEGPD and ZIDEGPD on the left
aic_bic_long$Model <- factor(aic_bic_long$Model,
                             levels = c("DEGPD", "ZIDEGPD", "NB", "ZINB")
)

# Compute mean of the best model (lowest minimum) per criterion
min_per_model <- aic_bic_long %>%
  group_by(Model, Criterion) %>%
  summarise(Min_Value = min(Value, na.rm = TRUE), .groups = "drop")

lowest_model_per_criterion <- min_per_model %>%
  group_by(Criterion) %>%
  filter(Min_Value == min(Min_Value)) %>%
  select(Criterion, Model)

mean_value_of_best_model <- aic_bic_long %>%
  inner_join(lowest_model_per_criterion, by = c("Criterion", "Model")) %>%
  group_by(Criterion) %>%
  summarise(Line_Value = mean(Value, na.rm = TRUE), .groups = "drop")

# ======= 1. Flask-Shaped Boxplot (varwidth) =======
ggplot(aic_bic_long, aes(x = Model, y = Value, fill = Model)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, varwidth = TRUE) +
  #geom_hline(data = mean_value_of_best_model,
  #aes(yintercept = Line_Value),
  #color = "red", linetype = "dashed", linewidth = 0.5,
  # inherit.aes = FALSE) +
  facet_wrap(~ Criterion, scales = "free_y") +
  labs(title = "Flask-Shaped Boxplots (Varwidth)\nDashed Red Line = Mean of Best Model",
       x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# ======= 2. Violin + Boxplot =======
p<-ggplot(aic_bic_long, aes(x = Model, y = Value, fill = Model)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.4, outlier.shape = TRUE, fill = "white", color = "black", linewidth = 0.2, alpha=0.6) +
  #geom_hline(data = mean_value_of_best_model,
  #           aes(yintercept = Line_Value),
  #          color = "red", linetype = "dashed", linewidth = 0.3,inherit.aes = FALSE) +
  facet_wrap(~ Criterion, scales = "free_y") +
  labs(title = "AIC and BIC plots of the fitted models at 5% low outliers",
       x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

p




# Save to PDF
ggsave("AIC_BIC_0.05_low_0.50_zero_degpd.pdf", plot = p, width = 6, height = 5)








###. == medium outliers with 0.05 ratio and 25% zeros

# Initialize storage matrices
results_matrix <- matrix(NA, nrow = 1000, ncol = 10)
colnames(results_matrix) <- c("Poisson_AIC", "Poisson_BIC",
                              "NB_AIC", "NB_BIC",
                              "ZINB_AIC", "ZINB_BIC",
                              "DEGPD_AIC", "DEGPD_BIC",
                              "ZIDEGPD_AIC", "ZIDEGPD_BIC")

set.seed(123)

for (i in 1:1000) {
  # Simulate data once
  sim_data <- simulate_one_rep(outlier_level = "medium", zero_inflation = 0.50, outlier_ratio = 0.05)
  
  # 1. Poisson Model
  poisson_model <- try(glm(Y ~ X1 + X2, family = poisson(link = "log"), data = sim_data), silent = TRUE)
  if (!inherits(poisson_model, "try-error")) {
    results_matrix[i, 1] <- AIC(poisson_model)
    results_matrix[i, 2] <- BIC(poisson_model)
  }
  
  # 2. Negative Binomial Model
  nb_model <- try(glm.nb(Y ~ X1 + X2, data = sim_data), silent = TRUE)
  if (!inherits(nb_model, "try-error")) {
    results_matrix[i, 3] <- AIC(nb_model)
    results_matrix[i, 4] <- BIC(nb_model)
  }
  
  # 3. ZINB Model
  zinb_model <- try(zeroinfl(Y ~ X1 + X2 | X1 + X2, data = sim_data, dist = "negbin"), silent = TRUE)
  if (!inherits(zinb_model, "try-error")) {
    results_matrix[i, 5] <- AIC(zinb_model)
    results_matrix[i, 6] <- BIC(zinb_model)
  }
  
  # 4. DEGPD Model
  fmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1)
  inits_degpd <- c(sqrt(sd(sim_data$Y)/2), sqrt(sd(sim_data$Y)/20), 0.3)
  degpd1 <- try(devgam(fmla_degpd1, data = sim_data, family = "degpd",
                       degpd.args = list(m = 1), trace = 0), silent = TRUE)
  if (!inherits(degpd1, "try-error")) {
    results_matrix[i, 7] <- AIC(degpd1)
    results_matrix[i, 8] <- BIC(degpd1)
  }
  
  # 5. ZIDEGPD Model
  zifmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1, logitpi = ~1)
  inits_zidegpd <- c(sd(sim_data$Y)/2, sd(sim_data$Y)/2, 0.2, mean(sim_data$Y == 0))
  zidegpd1 <- try(devgam(zifmla_degpd1, data = sim_data, family = "zidegpd",
                         degpd.args = list(m = 1), trace = 0, inits = inits_zidegpd), silent = TRUE)
  if (!inherits(zidegpd1, "try-error")) {
    results_matrix[i, 9] <- AIC(zidegpd1)
    results_matrix[i, 10] <- BIC(zidegpd1)
  }
  
  if (i %% 100 == 0) cat("Completed:", i, "\n")
}






# Compute means for each model
mean_aic_bic <- data.frame(
  Model = c("Poisson", "Negative Binomial", "ZINB", "DEGPD", "ZIDEGPD"),
  Mean_AIC = c(
    mean(results_matrix[, 1], na.rm = TRUE),
    mean(results_matrix[, 3], na.rm = TRUE),
    mean(results_matrix[, 5], na.rm = TRUE),
    mean(results_matrix[, 7], na.rm = TRUE),
    mean(results_matrix[, 9], na.rm = TRUE)
  ),
  Mean_BIC = c(
    mean(results_matrix[, 2], na.rm = TRUE),
    mean(results_matrix[, 4], na.rm = TRUE),
    mean(results_matrix[, 6], na.rm = TRUE),
    mean(results_matrix[, 8], na.rm = TRUE),
    mean(results_matrix[, 10], na.rm = TRUE)
  )
)

# Print summary table
print(mean_aic_bic, row.names = FALSE)

# Combine AIC/BIC matrices into a long-format dataframe
aic_bic_df <- data.frame(
  Model = rep(c( "NB", "ZINB", "DEGPD", "ZIDEGPD"), each = nrow(results_matrix)),
  AIC = c(#results_matrix[, 1],
    results_matrix[, 3],
    results_matrix[, 5],
    results_matrix[, 7],
    results_matrix[, 9]),
  BIC = c(#results_matrix[, 2],
    results_matrix[, 4],
    results_matrix[, 6],
    results_matrix[, 8],
    results_matrix[, 10])
)

# Convert to long format
aic_bic_long <- pivot_longer(aic_bic_df, cols = c(AIC, BIC), names_to = "Criterion", values_to = "Value")

# Reorder Model to show DEGPD and ZIDEGPD on the left
aic_bic_long$Model <- factor(aic_bic_long$Model,
                             levels = c("DEGPD", "ZIDEGPD", "NB", "ZINB")
)

# Compute mean of the best model (lowest minimum) per criterion
min_per_model <- aic_bic_long %>%
  group_by(Model, Criterion) %>%
  summarise(Min_Value = min(Value, na.rm = TRUE), .groups = "drop")

lowest_model_per_criterion <- min_per_model %>%
  group_by(Criterion) %>%
  filter(Min_Value == min(Min_Value)) %>%
  select(Criterion, Model)

mean_value_of_best_model <- aic_bic_long %>%
  inner_join(lowest_model_per_criterion, by = c("Criterion", "Model")) %>%
  group_by(Criterion) %>%
  summarise(Line_Value = mean(Value, na.rm = TRUE), .groups = "drop")

# ======= 1. Flask-Shaped Boxplot (varwidth) =======
ggplot(aic_bic_long, aes(x = Model, y = Value, fill = Model)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, varwidth = TRUE) +
  #geom_hline(data = mean_value_of_best_model,
  #aes(yintercept = Line_Value),
  #color = "red", linetype = "dashed", linewidth = 0.5,
  # inherit.aes = FALSE) +
  facet_wrap(~ Criterion, scales = "free_y") +
  labs(title = "Flask-Shaped Boxplots (Varwidth)\nDashed Red Line = Mean of Best Model",
       x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# ======= 2. Violin + Boxplot =======
p<-ggplot(aic_bic_long, aes(x = Model, y = Value, fill = Model)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.4, outlier.shape = TRUE, fill = "white", color = "black", linewidth = 0.2, alpha=0.6) +
  #geom_hline(data = mean_value_of_best_model,
  #           aes(yintercept = Line_Value),
  #          color = "red", linetype = "dashed", linewidth = 0.3,inherit.aes = FALSE) +
  facet_wrap(~ Criterion, scales = "free_y") +
  labs(title = "AIC and BIC plots of the fitted models at 5% medium outliers",
       x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

p




# Save to PDF
ggsave("AIC_BIC_0.05_medium_0.50_zero_degpd.pdf", plot = p, width = 6, height = 5)




###. == high outliers with 0.01 ratio and 25% zeros. ===

# Initialize storage matrices
results_matrix <- matrix(NA, nrow = 1000, ncol = 10)
colnames(results_matrix) <- c("Poisson_AIC", "Poisson_BIC",
                              "NB_AIC", "NB_BIC",
                              "ZINB_AIC", "ZINB_BIC",
                              "DEGPD_AIC", "DEGPD_BIC",
                              "ZIDEGPD_AIC", "ZIDEGPD_BIC")

set.seed(123)

for (i in 1:1000) {
  # Simulate data once
  sim_data <- simulate_one_rep(outlier_level = "high", zero_inflation = 0.50, outlier_ratio = 0.05)
  
  # 1. Poisson Model
  poisson_model <- try(glm(Y ~ X1 + X2, family = poisson(link = "log"), data = sim_data), silent = TRUE)
  if (!inherits(poisson_model, "try-error")) {
    results_matrix[i, 1] <- AIC(poisson_model)
    results_matrix[i, 2] <- BIC(poisson_model)
  }
  
  # 2. Negative Binomial Model
  nb_model <- try(glm.nb(Y ~ X1 + X2, data = sim_data), silent = TRUE)
  if (!inherits(nb_model, "try-error")) {
    results_matrix[i, 3] <- AIC(nb_model)
    results_matrix[i, 4] <- BIC(nb_model)
  }
  
  # 3. ZINB Model
  zinb_model <- try(zeroinfl(Y ~ X1 + X2 | X1 + X2, data = sim_data, dist = "negbin"), silent = TRUE)
  if (!inherits(zinb_model, "try-error")) {
    results_matrix[i, 5] <- AIC(zinb_model)
    results_matrix[i, 6] <- BIC(zinb_model)
  }
  
  # 4. DEGPD Model
  fmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1)
  inits_degpd <- c(sqrt(sd(sim_data$Y)/2), sqrt(sd(sim_data$Y)/20), 0.3)
  degpd1 <- try(devgam(fmla_degpd1, data = sim_data, family = "degpd",
                       degpd.args = list(m = 1), trace = 0), silent = TRUE)
  if (!inherits(degpd1, "try-error")) {
    results_matrix[i, 7] <- AIC(degpd1)
    results_matrix[i, 8] <- BIC(degpd1)
  }
  
  # 5. ZIDEGPD Model
  zifmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1, logitpi = ~1)
  inits_zidegpd <- c(sd(sim_data$Y)/2, sd(sim_data$Y)/2, 0.2, mean(sim_data$Y == 0))
  zidegpd1 <- try(devgam(zifmla_degpd1, data = sim_data, family = "zidegpd",
                         degpd.args = list(m = 1), trace = 0, inits = inits_zidegpd), silent = TRUE)
  if (!inherits(zidegpd1, "try-error")) {
    results_matrix[i, 9] <- AIC(zidegpd1)
    results_matrix[i, 10] <- BIC(zidegpd1)
  }
  
  if (i %% 100 == 0) cat("Completed:", i, "\n")
}



# Compute means for each model
mean_aic_bic <- data.frame(
  Model = c("Poisson", "Negative Binomial", "ZINB", "DEGPD", "ZIDEGPD"),
  Mean_AIC = c(
    mean(results_matrix[, 1], na.rm = TRUE),
    mean(results_matrix[, 3], na.rm = TRUE),
    mean(results_matrix[, 5], na.rm = TRUE),
    mean(results_matrix[, 7], na.rm = TRUE),
    mean(results_matrix[, 9], na.rm = TRUE)
  ),
  Mean_BIC = c(
    mean(results_matrix[, 2], na.rm = TRUE),
    mean(results_matrix[, 4], na.rm = TRUE),
    mean(results_matrix[, 6], na.rm = TRUE),
    mean(results_matrix[, 8], na.rm = TRUE),
    mean(results_matrix[, 10], na.rm = TRUE)
  )
)

# Print summary table
print(mean_aic_bic, row.names = FALSE)

# Combine AIC/BIC matrices into a long-format dataframe
aic_bic_df <- data.frame(
  Model = rep(c( "NB", "ZINB", "DEGPD", "ZIDEGPD"), each = nrow(results_matrix)),
  AIC = c(#results_matrix[, 1],
    results_matrix[, 3],
    results_matrix[, 5],
    results_matrix[, 7],
    results_matrix[, 9]),
  BIC = c(#results_matrix[, 2],
    results_matrix[, 4],
    results_matrix[, 6],
    results_matrix[, 8],
    results_matrix[, 10])
)

# Convert to long format
aic_bic_long <- pivot_longer(aic_bic_df, cols = c(AIC, BIC), names_to = "Criterion", values_to = "Value")

# Reorder Model to show DEGPD and ZIDEGPD on the left
aic_bic_long$Model <- factor(aic_bic_long$Model,
                             levels = c("DEGPD", "ZIDEGPD", "NB", "ZINB")
)

# Compute mean of the best model (lowest minimum) per criterion
min_per_model <- aic_bic_long %>%
  group_by(Model, Criterion) %>%
  summarise(Min_Value = min(Value, na.rm = TRUE), .groups = "drop")

lowest_model_per_criterion <- min_per_model %>%
  group_by(Criterion) %>%
  filter(Min_Value == min(Min_Value)) %>%
  select(Criterion, Model)

mean_value_of_best_model <- aic_bic_long %>%
  inner_join(lowest_model_per_criterion, by = c("Criterion", "Model")) %>%
  group_by(Criterion) %>%
  summarise(Line_Value = mean(Value, na.rm = TRUE), .groups = "drop")

# ======= 1. Flask-Shaped Boxplot (varwidth) =======
ggplot(aic_bic_long, aes(x = Model, y = Value, fill = Model)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, varwidth = TRUE) +
  #geom_hline(data = mean_value_of_best_model,
  #aes(yintercept = Line_Value),
  #color = "red", linetype = "dashed", linewidth = 0.5,
  # inherit.aes = FALSE) +
  facet_wrap(~ Criterion, scales = "free_y") +
  labs(title = "Flask-Shaped Boxplots (Varwidth)\nDashed Red Line = Mean of Best Model",
       x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# ======= 2. Violin + Boxplot =======
p<-ggplot(aic_bic_long, aes(x = Model, y = Value, fill = Model)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.4, outlier.shape = TRUE, fill = "white", color = "black", linewidth = 0.2, alpha=0.6) +
  #geom_hline(data = mean_value_of_best_model,
  #           aes(yintercept = Line_Value),
  #          color = "red", linetype = "dashed", linewidth = 0.3,inherit.aes = FALSE) +
  facet_wrap(~ Criterion, scales = "free_y") +
  labs(title = "AIC and BIC plots of the fitted models at 5% high outliers",
       x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

p




# Save to PDF
ggsave("AIC_BIC_0.05_high_0.50_zero_degpd.pdf", plot = p, width = 6, height = 5)






#done



##===================================================================================================================================================================



###. == Low outliers with 10% ratio and 25% zeros

# Initialize storage matrices
results_matrix <- matrix(NA, nrow = 1000, ncol = 10)
colnames(results_matrix) <- c("Poisson_AIC", "Poisson_BIC",
                              "NB_AIC", "NB_BIC",
                              "ZINB_AIC", "ZINB_BIC",
                              "DEGPD_AIC", "DEGPD_BIC",
                              "ZIDEGPD_AIC", "ZIDEGPD_BIC")

set.seed(123)

for (i in 1:1000) {
  # Simulate data once
  sim_data <- simulate_one_rep(outlier_level = "low", zero_inflation = 0.50, outlier_ratio = 0.10)
  
  # 1. Poisson Model
  poisson_model <- try(glm(Y ~ X1 + X2, family = poisson(link = "log"), data = sim_data), silent = TRUE)
  if (!inherits(poisson_model, "try-error")) {
    results_matrix[i, 1] <- AIC(poisson_model)
    results_matrix[i, 2] <- BIC(poisson_model)
  }
  
  # 2. Negative Binomial Model
  nb_model <- try(glm.nb(Y ~ X1 + X2, data = sim_data), silent = TRUE)
  if (!inherits(nb_model, "try-error")) {
    results_matrix[i, 3] <- AIC(nb_model)
    results_matrix[i, 4] <- BIC(nb_model)
  }
  
  # 3. ZINB Model
  zinb_model <- try(zeroinfl(Y ~ X1 + X2 | X1 + X2, data = sim_data, dist = "negbin"), silent = TRUE)
  if (!inherits(zinb_model, "try-error")) {
    results_matrix[i, 5] <- AIC(zinb_model)
    results_matrix[i, 6] <- BIC(zinb_model)
  }
  
  # 4. DEGPD Model
  fmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1)
  inits_degpd <- c(sqrt(sd(sim_data$Y)/2), sqrt(sd(sim_data$Y)/20), 0.3)
  degpd1 <- try(devgam(fmla_degpd1, data = sim_data, family = "degpd",
                       degpd.args = list(m = 1), trace = 0), silent = TRUE)
  if (!inherits(degpd1, "try-error")) {
    results_matrix[i, 7] <- AIC(degpd1)
    results_matrix[i, 8] <- BIC(degpd1)
  }
  
  # 5. ZIDEGPD Model
  zifmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1, logitpi = ~1)
  inits_zidegpd <- c(sd(sim_data$Y)/2, sd(sim_data$Y)/2, 0.2, mean(sim_data$Y == 0))
  zidegpd1 <- try(devgam(zifmla_degpd1, data = sim_data, family = "zidegpd",
                         degpd.args = list(m = 1), trace = 0, inits = inits_zidegpd), silent = TRUE)
  if (!inherits(zidegpd1, "try-error")) {
    results_matrix[i, 9] <- AIC(zidegpd1)
    results_matrix[i, 10] <- BIC(zidegpd1)
  }
  
  if (i %% 100 == 0) cat("Completed:", i, "\n")
}




# Compute means for each model
mean_aic_bic <- data.frame(
  Model = c("Poisson", "Negative Binomial", "ZINB", "DEGPD", "ZIDEGPD"),
  Mean_AIC = c(
    mean(results_matrix[, 1], na.rm = TRUE),
    mean(results_matrix[, 3], na.rm = TRUE),
    mean(results_matrix[, 5], na.rm = TRUE),
    mean(results_matrix[, 7], na.rm = TRUE),
    mean(results_matrix[, 9], na.rm = TRUE)
  ),
  Mean_BIC = c(
    mean(results_matrix[, 2], na.rm = TRUE),
    mean(results_matrix[, 4], na.rm = TRUE),
    mean(results_matrix[, 6], na.rm = TRUE),
    mean(results_matrix[, 8], na.rm = TRUE),
    mean(results_matrix[, 10], na.rm = TRUE)
  )
)

# Print summary table
print(mean_aic_bic, row.names = FALSE)

# Combine AIC/BIC matrices into a long-format dataframe
aic_bic_df <- data.frame(
  Model = rep(c( "NB", "ZINB", "DEGPD", "ZIDEGPD"), each = nrow(results_matrix)),
  AIC = c(#results_matrix[, 1],
    results_matrix[, 3],
    results_matrix[, 5],
    results_matrix[, 7],
    results_matrix[, 9]),
  BIC = c(#results_matrix[, 2],
    results_matrix[, 4],
    results_matrix[, 6],
    results_matrix[, 8],
    results_matrix[, 10])
)

# Convert to long format
aic_bic_long <- pivot_longer(aic_bic_df, cols = c(AIC, BIC), names_to = "Criterion", values_to = "Value")

# Reorder Model to show DEGPD and ZIDEGPD on the left
aic_bic_long$Model <- factor(aic_bic_long$Model,
                             levels = c("DEGPD", "ZIDEGPD", "NB", "ZINB")
)

# Compute mean of the best model (lowest minimum) per criterion
min_per_model <- aic_bic_long %>%
  group_by(Model, Criterion) %>%
  summarise(Min_Value = min(Value, na.rm = TRUE), .groups = "drop")

lowest_model_per_criterion <- min_per_model %>%
  group_by(Criterion) %>%
  filter(Min_Value == min(Min_Value)) %>%
  select(Criterion, Model)

mean_value_of_best_model <- aic_bic_long %>%
  inner_join(lowest_model_per_criterion, by = c("Criterion", "Model")) %>%
  group_by(Criterion) %>%
  summarise(Line_Value = mean(Value, na.rm = TRUE), .groups = "drop")

# ======= 1. Flask-Shaped Boxplot (varwidth) =======
ggplot(aic_bic_long, aes(x = Model, y = Value, fill = Model)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, varwidth = TRUE) +
  #geom_hline(data = mean_value_of_best_model,
  #aes(yintercept = Line_Value),
  #color = "red", linetype = "dashed", linewidth = 0.5,
  # inherit.aes = FALSE) +
  facet_wrap(~ Criterion, scales = "free_y") +
  labs(title = "Flask-Shaped Boxplots (Varwidth)\nDashed Red Line = Mean of Best Model",
       x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# ======= 2. Violin + Boxplot =======
p<-ggplot(aic_bic_long, aes(x = Model, y = Value, fill = Model)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.4, outlier.shape = TRUE, fill = "white", color = "black", linewidth = 0.2, alpha=0.6) +
  #geom_hline(data = mean_value_of_best_model,
  #           aes(yintercept = Line_Value),
  #          color = "red", linetype = "dashed", linewidth = 0.3,inherit.aes = FALSE) +
  facet_wrap(~ Criterion, scales = "free_y") +
  labs(title = "AIC and BIC plots of the fitted models at 10% low outliers",
       x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

p




# Save to PDF
ggsave("AIC_BIC_0.10_low_0.50_zero_degpd.pdf", plot = p, width = 6, height = 5)








###. == medium outliers with 0.10 ratio and 25% zeros

# Initialize storage matrices
results_matrix <- matrix(NA, nrow = 1000, ncol = 10)
colnames(results_matrix) <- c("Poisson_AIC", "Poisson_BIC",
                              "NB_AIC", "NB_BIC",
                              "ZINB_AIC", "ZINB_BIC",
                              "DEGPD_AIC", "DEGPD_BIC",
                              "ZIDEGPD_AIC", "ZIDEGPD_BIC")

set.seed(123)

for (i in 1:1000) {
  # Simulate data once
  sim_data <- simulate_one_rep(outlier_level = "medium", zero_inflation = 0.50, outlier_ratio = 0.10)
  
  # 1. Poisson Model
  poisson_model <- try(glm(Y ~ X1 + X2, family = poisson(link = "log"), data = sim_data), silent = TRUE)
  if (!inherits(poisson_model, "try-error")) {
    results_matrix[i, 1] <- AIC(poisson_model)
    results_matrix[i, 2] <- BIC(poisson_model)
  }
  
  # 2. Negative Binomial Model
  nb_model <- try(glm.nb(Y ~ X1 + X2, data = sim_data), silent = TRUE)
  if (!inherits(nb_model, "try-error")) {
    results_matrix[i, 3] <- AIC(nb_model)
    results_matrix[i, 4] <- BIC(nb_model)
  }
  
  # 3. ZINB Model
  zinb_model <- try(zeroinfl(Y ~ X1 + X2 | X1 + X2, data = sim_data, dist = "negbin"), silent = TRUE)
  if (!inherits(zinb_model, "try-error")) {
    results_matrix[i, 5] <- AIC(zinb_model)
    results_matrix[i, 6] <- BIC(zinb_model)
  }
  
  # 4. DEGPD Model
  fmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1)
  inits_degpd <- c(sqrt(sd(sim_data$Y)/2), sqrt(sd(sim_data$Y)/20), 0.3)
  degpd1 <- try(devgam(fmla_degpd1, data = sim_data, family = "degpd",
                       degpd.args = list(m = 1), trace = 0), silent = TRUE)
  if (!inherits(degpd1, "try-error")) {
    results_matrix[i, 7] <- AIC(degpd1)
    results_matrix[i, 8] <- BIC(degpd1)
  }
  
  # 5. ZIDEGPD Model
  zifmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1, logitpi = ~1)
  inits_zidegpd <- c(sd(sim_data$Y)/2, sd(sim_data$Y)/2, 0.2, mean(sim_data$Y == 0))
  zidegpd1 <- try(devgam(zifmla_degpd1, data = sim_data, family = "zidegpd",
                         degpd.args = list(m = 1), trace = 0, inits = inits_zidegpd), silent = TRUE)
  if (!inherits(zidegpd1, "try-error")) {
    results_matrix[i, 9] <- AIC(zidegpd1)
    results_matrix[i, 10] <- BIC(zidegpd1)
  }
  
  if (i %% 100 == 0) cat("Completed:", i, "\n")
}






# Compute means for each model
mean_aic_bic <- data.frame(
  Model = c("Poisson", "Negative Binomial", "ZINB", "DEGPD", "ZIDEGPD"),
  Mean_AIC = c(
    mean(results_matrix[, 1], na.rm = TRUE),
    mean(results_matrix[, 3], na.rm = TRUE),
    mean(results_matrix[, 5], na.rm = TRUE),
    mean(results_matrix[, 7], na.rm = TRUE),
    mean(results_matrix[, 9], na.rm = TRUE)
  ),
  Mean_BIC = c(
    mean(results_matrix[, 2], na.rm = TRUE),
    mean(results_matrix[, 4], na.rm = TRUE),
    mean(results_matrix[, 6], na.rm = TRUE),
    mean(results_matrix[, 8], na.rm = TRUE),
    mean(results_matrix[, 10], na.rm = TRUE)
  )
)

# Print summary table
print(mean_aic_bic, row.names = FALSE)

# Combine AIC/BIC matrices into a long-format dataframe
aic_bic_df <- data.frame(
  Model = rep(c( "NB", "ZINB", "DEGPD", "ZIDEGPD"), each = nrow(results_matrix)),
  AIC = c(#results_matrix[, 1],
    results_matrix[, 3],
    results_matrix[, 5],
    results_matrix[, 7],
    results_matrix[, 9]),
  BIC = c(#results_matrix[, 2],
    results_matrix[, 4],
    results_matrix[, 6],
    results_matrix[, 8],
    results_matrix[, 10])
)

# Convert to long format
aic_bic_long <- pivot_longer(aic_bic_df, cols = c(AIC, BIC), names_to = "Criterion", values_to = "Value")

# Reorder Model to show DEGPD and ZIDEGPD on the left
aic_bic_long$Model <- factor(aic_bic_long$Model,
                             levels = c("DEGPD", "ZIDEGPD", "NB", "ZINB")
)

# Compute mean of the best model (lowest minimum) per criterion
min_per_model <- aic_bic_long %>%
  group_by(Model, Criterion) %>%
  summarise(Min_Value = min(Value, na.rm = TRUE), .groups = "drop")

lowest_model_per_criterion <- min_per_model %>%
  group_by(Criterion) %>%
  filter(Min_Value == min(Min_Value)) %>%
  select(Criterion, Model)

mean_value_of_best_model <- aic_bic_long %>%
  inner_join(lowest_model_per_criterion, by = c("Criterion", "Model")) %>%
  group_by(Criterion) %>%
  summarise(Line_Value = mean(Value, na.rm = TRUE), .groups = "drop")

# ======= 1. Flask-Shaped Boxplot (varwidth) =======
ggplot(aic_bic_long, aes(x = Model, y = Value, fill = Model)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, varwidth = TRUE) +
  #geom_hline(data = mean_value_of_best_model,
  #aes(yintercept = Line_Value),
  #color = "red", linetype = "dashed", linewidth = 0.5,
  # inherit.aes = FALSE) +
  facet_wrap(~ Criterion, scales = "free_y") +
  labs(title = "Flask-Shaped Boxplots (Varwidth)\nDashed Red Line = Mean of Best Model",
       x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# ======= 2. Violin + Boxplot =======
p<-ggplot(aic_bic_long, aes(x = Model, y = Value, fill = Model)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.4, outlier.shape = TRUE, fill = "white", color = "black", linewidth = 0.2, alpha=0.6) +
  #geom_hline(data = mean_value_of_best_model,
  #           aes(yintercept = Line_Value),
  #          color = "red", linetype = "dashed", linewidth = 0.3,inherit.aes = FALSE) +
  facet_wrap(~ Criterion, scales = "free_y") +
  labs(title = "AIC and BIC plots of the fitted models at 10% medium outliers",
       x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

p




# Save to PDF
ggsave("AIC_BIC_0.10_medium_0.50_zero_degpd.pdf", plot = p, width = 6, height = 5)




###. == high outliers with 0.10 ratio and 25% zeros. ===

# Initialize storage matrices
results_matrix <- matrix(NA, nrow = 1000, ncol = 10)
colnames(results_matrix) <- c("Poisson_AIC", "Poisson_BIC",
                              "NB_AIC", "NB_BIC",
                              "ZINB_AIC", "ZINB_BIC",
                              "DEGPD_AIC", "DEGPD_BIC",
                              "ZIDEGPD_AIC", "ZIDEGPD_BIC")

set.seed(123)

for (i in 1:1000) {
  # Simulate data once
  sim_data <- simulate_one_rep(outlier_level = "high", zero_inflation = 0.50, outlier_ratio = 0.10)
  
  # 1. Poisson Model
  poisson_model <- try(glm(Y ~ X1 + X2, family = poisson(link = "log"), data = sim_data), silent = TRUE)
  if (!inherits(poisson_model, "try-error")) {
    results_matrix[i, 1] <- AIC(poisson_model)
    results_matrix[i, 2] <- BIC(poisson_model)
  }
  
  # 2. Negative Binomial Model
  nb_model <- try(glm.nb(Y ~ X1 + X2, data = sim_data), silent = TRUE)
  if (!inherits(nb_model, "try-error")) {
    results_matrix[i, 3] <- AIC(nb_model)
    results_matrix[i, 4] <- BIC(nb_model)
  }
  
  # 3. ZINB Model
  zinb_model <- try(zeroinfl(Y ~ X1 + X2 | X1 + X2, data = sim_data, dist = "negbin"), silent = TRUE)
  if (!inherits(zinb_model, "try-error")) {
    results_matrix[i, 5] <- AIC(zinb_model)
    results_matrix[i, 6] <- BIC(zinb_model)
  }
  
  # 4. DEGPD Model
  fmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1)
  inits_degpd <- c(sqrt(sd(sim_data$Y)/2), sqrt(sd(sim_data$Y)/20), 0.3)
  degpd1 <- try(devgam(fmla_degpd1, data = sim_data, family = "degpd",
                       degpd.args = list(m = 1), trace = 0), silent = TRUE)
  if (!inherits(degpd1, "try-error")) {
    results_matrix[i, 7] <- AIC(degpd1)
    results_matrix[i, 8] <- BIC(degpd1)
  }
  
  # 5. ZIDEGPD Model
  zifmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1, logitpi = ~1)
  inits_zidegpd <- c(sd(sim_data$Y)/2, sd(sim_data$Y)/2, 0.2, mean(sim_data$Y == 0))
  zidegpd1 <- try(devgam(zifmla_degpd1, data = sim_data, family = "zidegpd",
                         degpd.args = list(m = 1), trace = 0, inits = inits_zidegpd), silent = TRUE)
  if (!inherits(zidegpd1, "try-error")) {
    results_matrix[i, 9] <- AIC(zidegpd1)
    results_matrix[i, 10] <- BIC(zidegpd1)
  }
  
  if (i %% 100 == 0) cat("Completed:", i, "\n")
}



# Compute means for each model
mean_aic_bic <- data.frame(
  Model = c("Poisson", "Negative Binomial", "ZINB", "DEGPD", "ZIDEGPD"),
  Mean_AIC = c(
    mean(results_matrix[, 1], na.rm = TRUE),
    mean(results_matrix[, 3], na.rm = TRUE),
    mean(results_matrix[, 5], na.rm = TRUE),
    mean(results_matrix[, 7], na.rm = TRUE),
    mean(results_matrix[, 9], na.rm = TRUE)
  ),
  Mean_BIC = c(
    mean(results_matrix[, 2], na.rm = TRUE),
    mean(results_matrix[, 4], na.rm = TRUE),
    mean(results_matrix[, 6], na.rm = TRUE),
    mean(results_matrix[, 8], na.rm = TRUE),
    mean(results_matrix[, 10], na.rm = TRUE)
  )
)

# Print summary table
print(mean_aic_bic, row.names = FALSE)

# Combine AIC/BIC matrices into a long-format dataframe
aic_bic_df <- data.frame(
  Model = rep(c( "NB", "ZINB", "DEGPD", "ZIDEGPD"), each = nrow(results_matrix)),
  AIC = c(#results_matrix[, 1],
    results_matrix[, 3],
    results_matrix[, 5],
    results_matrix[, 7],
    results_matrix[, 9]),
  BIC = c(#results_matrix[, 2],
    results_matrix[, 4],
    results_matrix[, 6],
    results_matrix[, 8],
    results_matrix[, 10])
)

# Convert to long format
aic_bic_long <- pivot_longer(aic_bic_df, cols = c(AIC, BIC), names_to = "Criterion", values_to = "Value")

# Reorder Model to show DEGPD and ZIDEGPD on the left
aic_bic_long$Model <- factor(aic_bic_long$Model,
                             levels = c("DEGPD", "ZIDEGPD", "NB", "ZINB")
)

# Compute mean of the best model (lowest minimum) per criterion
min_per_model <- aic_bic_long %>%
  group_by(Model, Criterion) %>%
  summarise(Min_Value = min(Value, na.rm = TRUE), .groups = "drop")

lowest_model_per_criterion <- min_per_model %>%
  group_by(Criterion) %>%
  filter(Min_Value == min(Min_Value)) %>%
  select(Criterion, Model)

mean_value_of_best_model <- aic_bic_long %>%
  inner_join(lowest_model_per_criterion, by = c("Criterion", "Model")) %>%
  group_by(Criterion) %>%
  summarise(Line_Value = mean(Value, na.rm = TRUE), .groups = "drop")

# ======= 1. Flask-Shaped Boxplot (varwidth) =======
ggplot(aic_bic_long, aes(x = Model, y = Value, fill = Model)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, varwidth = TRUE) +
  #geom_hline(data = mean_value_of_best_model,
  #aes(yintercept = Line_Value),
  #color = "red", linetype = "dashed", linewidth = 0.5,
  # inherit.aes = FALSE) +
  facet_wrap(~ Criterion, scales = "free_y") +
  labs(title = "Flask-Shaped Boxplots (Varwidth)\nDashed Red Line = Mean of Best Model",
       x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# ======= 2. Violin + Boxplot =======
p<-ggplot(aic_bic_long, aes(x = Model, y = Value, fill = Model)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.4, outlier.shape = TRUE, fill = "white", color = "black", linewidth = 0.2, alpha=0.6) +
  #geom_hline(data = mean_value_of_best_model,
  #           aes(yintercept = Line_Value),
  #          color = "red", linetype = "dashed", linewidth = 0.3,inherit.aes = FALSE) +
  facet_wrap(~ Criterion, scales = "free_y") +
  labs(title = "AIC and BIC plots of the fitted models at 10% high outliers",
       x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

p




# Save to PDF
ggsave("AIC_BIC_0.10_high_0.50_zero_degpd.pdf", plot = p, width = 6, height = 5)























#done

##. 50 % zeros========================================================================================================================





###. == Low outliers with 0.01 ratio and 25% zeros

# Initialize storage matrices
results_matrix <- matrix(NA, nrow = 1000, ncol = 10)
colnames(results_matrix) <- c("Poisson_AIC", "Poisson_BIC",
                              "NB_AIC", "NB_BIC",
                              "ZINB_AIC", "ZINB_BIC",
                              "DEGPD_AIC", "DEGPD_BIC",
                              "ZIDEGPD_AIC", "ZIDEGPD_BIC")

set.seed(123)

for (i in 1:1000) {
  # Simulate data once
  sim_data <- simulate_one_rep(outlier_level = "low", zero_inflation = 0.70, outlier_ratio = 0.01)
  
  # 1. Poisson Model
  poisson_model <- try(glm(Y ~ X1 + X2, family = poisson(link = "log"), data = sim_data), silent = TRUE)
  if (!inherits(poisson_model, "try-error")) {
    results_matrix[i, 1] <- AIC(poisson_model)
    results_matrix[i, 2] <- BIC(poisson_model)
  }
  
  # 2. Negative Binomial Model
  nb_model <- try(glm.nb(Y ~ X1 + X2, data = sim_data), silent = TRUE)
  if (!inherits(nb_model, "try-error")) {
    results_matrix[i, 3] <- AIC(nb_model)
    results_matrix[i, 4] <- BIC(nb_model)
  }
  
  # 3. ZINB Model
  zinb_model <- try(zeroinfl(Y ~ X1 + X2 | X1 + X2, data = sim_data, dist = "negbin"), silent = TRUE)
  if (!inherits(zinb_model, "try-error")) {
    results_matrix[i, 5] <- AIC(zinb_model)
    results_matrix[i, 6] <- BIC(zinb_model)
  }
  
  # 4. DEGPD Model
  fmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1)
  inits_degpd <- c(sqrt(sd(sim_data$Y)/2), sqrt(sd(sim_data$Y)/20), 0.3)
  degpd1 <- try(devgam(fmla_degpd1, data = sim_data, family = "degpd",
                       degpd.args = list(m = 1), trace = 0), silent = TRUE)
  if (!inherits(degpd1, "try-error")) {
    results_matrix[i, 7] <- AIC(degpd1)
    results_matrix[i, 8] <- BIC(degpd1)
  }
  
  # 5. ZIDEGPD Model
  zifmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1, logitpi = ~1)
  inits_zidegpd <- c(sd(sim_data$Y)/2, sd(sim_data$Y)/2, 0.2, mean(sim_data$Y == 0))
  zidegpd1 <- try(devgam(zifmla_degpd1, data = sim_data, family = "zidegpd",
                         degpd.args = list(m = 1), trace = 0, inits = inits_zidegpd), silent = TRUE)
  if (!inherits(zidegpd1, "try-error")) {
    results_matrix[i, 9] <- AIC(zidegpd1)
    results_matrix[i, 10] <- BIC(zidegpd1)
  }
  
  if (i %% 100 == 0) cat("Completed:", i, "\n")
}






# Compute means for each model
mean_aic_bic <- data.frame(
  Model = c("Poisson", "Negative Binomial", "ZINB", "DEGPD", "ZIDEGPD"),
  Mean_AIC = c(
    mean(results_matrix[, 1], na.rm = TRUE),
    mean(results_matrix[, 3], na.rm = TRUE),
    mean(results_matrix[, 5], na.rm = TRUE),
    mean(results_matrix[, 7], na.rm = TRUE),
    mean(results_matrix[, 9], na.rm = TRUE)
  ),
  Mean_BIC = c(
    mean(results_matrix[, 2], na.rm = TRUE),
    mean(results_matrix[, 4], na.rm = TRUE),
    mean(results_matrix[, 6], na.rm = TRUE),
    mean(results_matrix[, 8], na.rm = TRUE),
    mean(results_matrix[, 10], na.rm = TRUE)
  )
)

# Print summary table
print(mean_aic_bic, row.names = FALSE)

# Combine AIC/BIC matrices into a long-format dataframe
aic_bic_df <- data.frame(
  Model = rep(c( "NB", "ZINB", "DEGPD", "ZIDEGPD"), each = nrow(results_matrix)),
  AIC = c(#results_matrix[, 1],
    results_matrix[, 3],
    results_matrix[, 5],
    results_matrix[, 7],
    results_matrix[, 9]),
  BIC = c(#results_matrix[, 2],
    results_matrix[, 4],
    results_matrix[, 6],
    results_matrix[, 8],
    results_matrix[, 10])
)

# Convert to long format
aic_bic_long <- pivot_longer(aic_bic_df, cols = c(AIC, BIC), names_to = "Criterion", values_to = "Value")

# Reorder Model to show DEGPD and ZIDEGPD on the left
aic_bic_long$Model <- factor(aic_bic_long$Model,
                             levels = c("DEGPD", "ZIDEGPD", "NB", "ZINB")
)

# Compute mean of the best model (lowest minimum) per criterion
min_per_model <- aic_bic_long %>%
  group_by(Model, Criterion) %>%
  summarise(Min_Value = min(Value, na.rm = TRUE), .groups = "drop")

lowest_model_per_criterion <- min_per_model %>%
  group_by(Criterion) %>%
  filter(Min_Value == min(Min_Value)) %>%
  select(Criterion, Model)

mean_value_of_best_model <- aic_bic_long %>%
  inner_join(lowest_model_per_criterion, by = c("Criterion", "Model")) %>%
  group_by(Criterion) %>%
  summarise(Line_Value = mean(Value, na.rm = TRUE), .groups = "drop")

# ======= 1. Flask-Shaped Boxplot (varwidth) =======
ggplot(aic_bic_long, aes(x = Model, y = Value, fill = Model)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, varwidth = TRUE) +
  #geom_hline(data = mean_value_of_best_model,
  #aes(yintercept = Line_Value),
  #color = "red", linetype = "dashed", linewidth = 0.5,
  # inherit.aes = FALSE) +
  facet_wrap(~ Criterion, scales = "free_y") +
  labs(title = "Flask-Shaped Boxplots (Varwidth)\nDashed Red Line = Mean of Best Model",
       x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# ======= 2. Violin + Boxplot =======
p<-ggplot(aic_bic_long, aes(x = Model, y = Value, fill = Model)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.4, outlier.shape = TRUE, fill = "white", color = "black", linewidth = 0.2, alpha=0.6) +
  #geom_hline(data = mean_value_of_best_model,
  #           aes(yintercept = Line_Value),
  #          color = "red", linetype = "dashed", linewidth = 0.3,inherit.aes = FALSE) +
  facet_wrap(~ Criterion, scales = "free_y") +
  labs(title = "AIC and BIC plots of the fitted models at 1% low outliers",
       x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

p




# Save to PDF
ggsave("AIC_BIC_0.01_low_0.70_zero_degpd.pdf", plot = p, width = 6, height = 5)








###. == medium outliers with 0.01 ratio and 25% zeros

# Initialize storage matrices
results_matrix <- matrix(NA, nrow = 1000, ncol = 10)
colnames(results_matrix) <- c("Poisson_AIC", "Poisson_BIC",
                              "NB_AIC", "NB_BIC",
                              "ZINB_AIC", "ZINB_BIC",
                              "DEGPD_AIC", "DEGPD_BIC",
                              "ZIDEGPD_AIC", "ZIDEGPD_BIC")

set.seed(123)

for (i in 1:1000) {
  # Simulate data once
  sim_data <- simulate_one_rep(outlier_level = "medium", zero_inflation = 0.70, outlier_ratio = 0.01)
  
  # 1. Poisson Model
  poisson_model <- try(glm(Y ~ X1 + X2, family = poisson(link = "log"), data = sim_data), silent = TRUE)
  if (!inherits(poisson_model, "try-error")) {
    results_matrix[i, 1] <- AIC(poisson_model)
    results_matrix[i, 2] <- BIC(poisson_model)
  }
  
  # 2. Negative Binomial Model
  nb_model <- try(glm.nb(Y ~ X1 + X2, data = sim_data), silent = TRUE)
  if (!inherits(nb_model, "try-error")) {
    results_matrix[i, 3] <- AIC(nb_model)
    results_matrix[i, 4] <- BIC(nb_model)
  }
  
  # 3. ZINB Model
  zinb_model <- try(zeroinfl(Y ~ X1 + X2 | X1 + X2, data = sim_data, dist = "negbin"), silent = TRUE)
  if (!inherits(zinb_model, "try-error")) {
    results_matrix[i, 5] <- AIC(zinb_model)
    results_matrix[i, 6] <- BIC(zinb_model)
  }
  
  # 4. DEGPD Model
  fmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1)
  inits_degpd <- c(sqrt(sd(sim_data$Y)/2), sqrt(sd(sim_data$Y)/20), 0.3)
  degpd1 <- try(devgam(fmla_degpd1, data = sim_data, family = "degpd",
                       degpd.args = list(m = 1), trace = 0), silent = TRUE)
  if (!inherits(degpd1, "try-error")) {
    results_matrix[i, 7] <- AIC(degpd1)
    results_matrix[i, 8] <- BIC(degpd1)
  }
  
  # 5. ZIDEGPD Model
  zifmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1, logitpi = ~1)
  inits_zidegpd <- c(sd(sim_data$Y)/2, sd(sim_data$Y)/2, 0.2, mean(sim_data$Y == 0))
  zidegpd1 <- try(devgam(zifmla_degpd1, data = sim_data, family = "zidegpd",
                         degpd.args = list(m = 1), trace = 0, inits = inits_zidegpd), silent = TRUE)
  if (!inherits(zidegpd1, "try-error")) {
    results_matrix[i, 9] <- AIC(zidegpd1)
    results_matrix[i, 10] <- BIC(zidegpd1)
  }
  
  if (i %% 100 == 0) cat("Completed:", i, "\n")
}






# Compute means for each model
mean_aic_bic <- data.frame(
  Model = c("Poisson", "Negative Binomial", "ZINB", "DEGPD", "ZIDEGPD"),
  Mean_AIC = c(
    mean(results_matrix[, 1], na.rm = TRUE),
    mean(results_matrix[, 3], na.rm = TRUE),
    mean(results_matrix[, 5], na.rm = TRUE),
    mean(results_matrix[, 7], na.rm = TRUE),
    mean(results_matrix[, 9], na.rm = TRUE)
  ),
  Mean_BIC = c(
    mean(results_matrix[, 2], na.rm = TRUE),
    mean(results_matrix[, 4], na.rm = TRUE),
    mean(results_matrix[, 6], na.rm = TRUE),
    mean(results_matrix[, 8], na.rm = TRUE),
    mean(results_matrix[, 10], na.rm = TRUE)
  )
)

# Print summary table
print(mean_aic_bic, row.names = FALSE)

# Combine AIC/BIC matrices into a long-format dataframe
aic_bic_df <- data.frame(
  Model = rep(c( "NB", "ZINB", "DEGPD", "ZIDEGPD"), each = nrow(results_matrix)),
  AIC = c(#results_matrix[, 1],
    results_matrix[, 3],
    results_matrix[, 5],
    results_matrix[, 7],
    results_matrix[, 9]),
  BIC = c(#results_matrix[, 2],
    results_matrix[, 4],
    results_matrix[, 6],
    results_matrix[, 8],
    results_matrix[, 10])
)

# Convert to long format
aic_bic_long <- pivot_longer(aic_bic_df, cols = c(AIC, BIC), names_to = "Criterion", values_to = "Value")

# Reorder Model to show DEGPD and ZIDEGPD on the left
aic_bic_long$Model <- factor(aic_bic_long$Model,
                             levels = c("DEGPD", "ZIDEGPD", "NB", "ZINB")
)

# Compute mean of the best model (lowest minimum) per criterion
min_per_model <- aic_bic_long %>%
  group_by(Model, Criterion) %>%
  summarise(Min_Value = min(Value, na.rm = TRUE), .groups = "drop")

lowest_model_per_criterion <- min_per_model %>%
  group_by(Criterion) %>%
  filter(Min_Value == min(Min_Value)) %>%
  select(Criterion, Model)

mean_value_of_best_model <- aic_bic_long %>%
  inner_join(lowest_model_per_criterion, by = c("Criterion", "Model")) %>%
  group_by(Criterion) %>%
  summarise(Line_Value = mean(Value, na.rm = TRUE), .groups = "drop")

# ======= 1. Flask-Shaped Boxplot (varwidth) =======
ggplot(aic_bic_long, aes(x = Model, y = Value, fill = Model)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, varwidth = TRUE) +
  #geom_hline(data = mean_value_of_best_model,
  #aes(yintercept = Line_Value),
  #color = "red", linetype = "dashed", linewidth = 0.5,
  # inherit.aes = FALSE) +
  facet_wrap(~ Criterion, scales = "free_y") +
  labs(title = "Flask-Shaped Boxplots (Varwidth)\nDashed Red Line = Mean of Best Model",
       x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# ======= 2. Violin + Boxplot =======
p<-ggplot(aic_bic_long, aes(x = Model, y = Value, fill = Model)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.4, outlier.shape = TRUE, fill = "white", color = "black", linewidth = 0.2, alpha=0.6) +
  #geom_hline(data = mean_value_of_best_model,
  #           aes(yintercept = Line_Value),
  #          color = "red", linetype = "dashed", linewidth = 0.3,inherit.aes = FALSE) +
  facet_wrap(~ Criterion, scales = "free_y") +
  labs(title = "AIC and BIC plots of the fitted models at 1% medium outliers",
       x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

p




# Save to PDF
ggsave("AIC_BIC_0.01_medium_0.70_zero_degpd.pdf", plot = p, width = 6, height = 5)




###. == high outliers with 0.01 ratio and 25% zeros. ===

# Initialize storage matrices
results_matrix <- matrix(NA, nrow = 1000, ncol = 10)
colnames(results_matrix) <- c("Poisson_AIC", "Poisson_BIC",
                              "NB_AIC", "NB_BIC",
                              "ZINB_AIC", "ZINB_BIC",
                              "DEGPD_AIC", "DEGPD_BIC",
                              "ZIDEGPD_AIC", "ZIDEGPD_BIC")

set.seed(123)

for (i in 1:1000) {
  # Simulate data once
  sim_data <- simulate_one_rep(outlier_level = "high", zero_inflation = 0.70, outlier_ratio = 0.01)
  
  # 1. Poisson Model
  poisson_model <- try(glm(Y ~ X1 + X2, family = poisson(link = "log"), data = sim_data), silent = TRUE)
  if (!inherits(poisson_model, "try-error")) {
    results_matrix[i, 1] <- AIC(poisson_model)
    results_matrix[i, 2] <- BIC(poisson_model)
  }
  
  # 2. Negative Binomial Model
  nb_model <- try(glm.nb(Y ~ X1 + X2, data = sim_data), silent = TRUE)
  if (!inherits(nb_model, "try-error")) {
    results_matrix[i, 3] <- AIC(nb_model)
    results_matrix[i, 4] <- BIC(nb_model)
  }
  
  # 3. ZINB Model
  zinb_model <- try(zeroinfl(Y ~ X1 + X2 | X1 + X2, data = sim_data, dist = "negbin"), silent = TRUE)
  if (!inherits(zinb_model, "try-error")) {
    results_matrix[i, 5] <- AIC(zinb_model)
    results_matrix[i, 6] <- BIC(zinb_model)
  }
  
  # 4. DEGPD Model
  fmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1)
  inits_degpd <- c(sqrt(sd(sim_data$Y)/2), sqrt(sd(sim_data$Y)/20), 0.3)
  degpd1 <- try(devgam(fmla_degpd1, data = sim_data, family = "degpd",
                       degpd.args = list(m = 1), trace = 0), silent = TRUE)
  if (!inherits(degpd1, "try-error")) {
    results_matrix[i, 7] <- AIC(degpd1)
    results_matrix[i, 8] <- BIC(degpd1)
  }
  
  # 5. ZIDEGPD Model
  zifmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1, logitpi = ~1)
  inits_zidegpd <- c(sd(sim_data$Y)/2, sd(sim_data$Y)/2, 0.2, mean(sim_data$Y == 0))
  zidegpd1 <- try(devgam(zifmla_degpd1, data = sim_data, family = "zidegpd",
                         degpd.args = list(m = 1), trace = 0, inits = inits_zidegpd), silent = TRUE)
  if (!inherits(zidegpd1, "try-error")) {
    results_matrix[i, 9] <- AIC(zidegpd1)
    results_matrix[i, 10] <- BIC(zidegpd1)
  }
  
  if (i %% 100 == 0) cat("Completed:", i, "\n")
}



# Compute means for each model
mean_aic_bic <- data.frame(
  Model = c("Poisson", "Negative Binomial", "ZINB", "DEGPD", "ZIDEGPD"),
  Mean_AIC = c(
    mean(results_matrix[, 1], na.rm = TRUE),
    mean(results_matrix[, 3], na.rm = TRUE),
    mean(results_matrix[, 5], na.rm = TRUE),
    mean(results_matrix[, 7], na.rm = TRUE),
    mean(results_matrix[, 9], na.rm = TRUE)
  ),
  Mean_BIC = c(
    mean(results_matrix[, 2], na.rm = TRUE),
    mean(results_matrix[, 4], na.rm = TRUE),
    mean(results_matrix[, 6], na.rm = TRUE),
    mean(results_matrix[, 8], na.rm = TRUE),
    mean(results_matrix[, 10], na.rm = TRUE)
  )
)

# Print summary table
print(mean_aic_bic, row.names = FALSE)

# Combine AIC/BIC matrices into a long-format dataframe
aic_bic_df <- data.frame(
  Model = rep(c( "NB", "ZINB", "DEGPD", "ZIDEGPD"), each = nrow(results_matrix)),
  AIC = c(#results_matrix[, 1],
    results_matrix[, 3],
    results_matrix[, 5],
    results_matrix[, 7],
    results_matrix[, 9]),
  BIC = c(#results_matrix[, 2],
    results_matrix[, 4],
    results_matrix[, 6],
    results_matrix[, 8],
    results_matrix[, 10])
)

# Convert to long format
aic_bic_long <- pivot_longer(aic_bic_df, cols = c(AIC, BIC), names_to = "Criterion", values_to = "Value")

# Reorder Model to show DEGPD and ZIDEGPD on the left
aic_bic_long$Model <- factor(aic_bic_long$Model,
                             levels = c("DEGPD", "ZIDEGPD", "NB", "ZINB")
)

# Compute mean of the best model (lowest minimum) per criterion
min_per_model <- aic_bic_long %>%
  group_by(Model, Criterion) %>%
  summarise(Min_Value = min(Value, na.rm = TRUE), .groups = "drop")

lowest_model_per_criterion <- min_per_model %>%
  group_by(Criterion) %>%
  filter(Min_Value == min(Min_Value)) %>%
  select(Criterion, Model)

mean_value_of_best_model <- aic_bic_long %>%
  inner_join(lowest_model_per_criterion, by = c("Criterion", "Model")) %>%
  group_by(Criterion) %>%
  summarise(Line_Value = mean(Value, na.rm = TRUE), .groups = "drop")

# ======= 1. Flask-Shaped Boxplot (varwidth) =======
ggplot(aic_bic_long, aes(x = Model, y = Value, fill = Model)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, varwidth = TRUE) +
  #geom_hline(data = mean_value_of_best_model,
  #aes(yintercept = Line_Value),
  #color = "red", linetype = "dashed", linewidth = 0.5,
  # inherit.aes = FALSE) +
  facet_wrap(~ Criterion, scales = "free_y") +
  labs(title = "Flask-Shaped Boxplots (Varwidth)\nDashed Red Line = Mean of Best Model",
       x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# ======= 2. Violin + Boxplot =======
p<-ggplot(aic_bic_long, aes(x = Model, y = Value, fill = Model)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.4, outlier.shape = TRUE, fill = "white", color = "black", linewidth = 0.2, alpha=0.6) +
  #geom_hline(data = mean_value_of_best_model,
  #           aes(yintercept = Line_Value),
  #          color = "red", linetype = "dashed", linewidth = 0.3,inherit.aes = FALSE) +
  facet_wrap(~ Criterion, scales = "free_y") +
  labs(title = "AIC and BIC plots of the fitted models at 1% high outliers",
       x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

p




# Save to PDF
ggsave("AIC_BIC_0.01_high_0.70_zero_degpd.pdf", plot = p, width = 6, height = 5)






##==============================================================================================================================================






###. == Low outliers with 5% ratio and 25% zeros

# Initialize storage matrices
results_matrix <- matrix(NA, nrow = 1000, ncol = 10)
colnames(results_matrix) <- c("Poisson_AIC", "Poisson_BIC",
                              "NB_AIC", "NB_BIC",
                              "ZINB_AIC", "ZINB_BIC",
                              "DEGPD_AIC", "DEGPD_BIC",
                              "ZIDEGPD_AIC", "ZIDEGPD_BIC")

set.seed(123)

for (i in 1:1000) {
  # Simulate data once
  sim_data <- simulate_one_rep(outlier_level = "low", zero_inflation = 0.70, outlier_ratio = 0.05)
  
  # 1. Poisson Model
  poisson_model <- try(glm(Y ~ X1 + X2, family = poisson(link = "log"), data = sim_data), silent = TRUE)
  if (!inherits(poisson_model, "try-error")) {
    results_matrix[i, 1] <- AIC(poisson_model)
    results_matrix[i, 2] <- BIC(poisson_model)
  }
  
  # 2. Negative Binomial Model
  nb_model <- try(glm.nb(Y ~ X1 + X2, data = sim_data), silent = TRUE)
  if (!inherits(nb_model, "try-error")) {
    results_matrix[i, 3] <- AIC(nb_model)
    results_matrix[i, 4] <- BIC(nb_model)
  }
  
  # 3. ZINB Model
  zinb_model <- try(zeroinfl(Y ~ X1 + X2 | X1 + X2, data = sim_data, dist = "negbin"), silent = TRUE)
  if (!inherits(zinb_model, "try-error")) {
    results_matrix[i, 5] <- AIC(zinb_model)
    results_matrix[i, 6] <- BIC(zinb_model)
  }
  
  # 4. DEGPD Model
  fmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1)
  inits_degpd <- c(sqrt(sd(sim_data$Y)/2), sqrt(sd(sim_data$Y)/20), 0.3)
  degpd1 <- try(devgam(fmla_degpd1, data = sim_data, family = "degpd",
                       degpd.args = list(m = 1), trace = 0), silent = TRUE)
  if (!inherits(degpd1, "try-error")) {
    results_matrix[i, 7] <- AIC(degpd1)
    results_matrix[i, 8] <- BIC(degpd1)
  }
  
  # 5. ZIDEGPD Model
  zifmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1, logitpi = ~1)
  inits_zidegpd <- c(sd(sim_data$Y)/2, sd(sim_data$Y)/2, 0.2, mean(sim_data$Y == 0))
  zidegpd1 <- try(devgam(zifmla_degpd1, data = sim_data, family = "zidegpd",
                         degpd.args = list(m = 1), trace = 0, inits = inits_zidegpd), silent = TRUE)
  if (!inherits(zidegpd1, "try-error")) {
    results_matrix[i, 9] <- AIC(zidegpd1)
    results_matrix[i, 10] <- BIC(zidegpd1)
  }
  
  if (i %% 100 == 0) cat("Completed:", i, "\n")
}







# Load required packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggdist)
library(ggbeeswarm)

# Compute means for each model
mean_aic_bic <- data.frame(
  Model = c("Poisson", "Negative Binomial", "ZINB", "DEGPD", "ZIDEGPD"),
  Mean_AIC = c(
    mean(results_matrix[, 1], na.rm = TRUE),
    mean(results_matrix[, 3], na.rm = TRUE),
    mean(results_matrix[, 5], na.rm = TRUE),
    mean(results_matrix[, 7], na.rm = TRUE),
    mean(results_matrix[, 9], na.rm = TRUE)
  ),
  Mean_BIC = c(
    mean(results_matrix[, 2], na.rm = TRUE),
    mean(results_matrix[, 4], na.rm = TRUE),
    mean(results_matrix[, 6], na.rm = TRUE),
    mean(results_matrix[, 8], na.rm = TRUE),
    mean(results_matrix[, 10], na.rm = TRUE)
  )
)

# Print summary table
print(mean_aic_bic, row.names = FALSE)

# Combine AIC/BIC matrices into a long-format dataframe
aic_bic_df <- data.frame(
  Model = rep(c( "NB", "ZINB", "DEGPD", "ZIDEGPD"), each = nrow(results_matrix)),
  AIC = c(#results_matrix[, 1],
    results_matrix[, 3],
    results_matrix[, 5],
    results_matrix[, 7],
    results_matrix[, 9]),
  BIC = c(#results_matrix[, 2],
    results_matrix[, 4],
    results_matrix[, 6],
    results_matrix[, 8],
    results_matrix[, 10])
)

# Convert to long format
aic_bic_long <- pivot_longer(aic_bic_df, cols = c(AIC, BIC), names_to = "Criterion", values_to = "Value")

# Reorder Model to show DEGPD and ZIDEGPD on the left
aic_bic_long$Model <- factor(aic_bic_long$Model,
                             levels = c("DEGPD", "ZIDEGPD", "NB", "ZINB")
)

# Compute mean of the best model (lowest minimum) per criterion
min_per_model <- aic_bic_long %>%
  group_by(Model, Criterion) %>%
  summarise(Min_Value = min(Value, na.rm = TRUE), .groups = "drop")

lowest_model_per_criterion <- min_per_model %>%
  group_by(Criterion) %>%
  filter(Min_Value == min(Min_Value)) %>%
  select(Criterion, Model)

mean_value_of_best_model <- aic_bic_long %>%
  inner_join(lowest_model_per_criterion, by = c("Criterion", "Model")) %>%
  group_by(Criterion) %>%
  summarise(Line_Value = mean(Value, na.rm = TRUE), .groups = "drop")

# ======= 1. Flask-Shaped Boxplot (varwidth) =======
ggplot(aic_bic_long, aes(x = Model, y = Value, fill = Model)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, varwidth = TRUE) +
  #geom_hline(data = mean_value_of_best_model,
  #aes(yintercept = Line_Value),
  #color = "red", linetype = "dashed", linewidth = 0.5,
  # inherit.aes = FALSE) +
  facet_wrap(~ Criterion, scales = "free_y") +
  labs(title = "Flask-Shaped Boxplots (Varwidth)\nDashed Red Line = Mean of Best Model",
       x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# ======= 2. Violin + Boxplot =======
p<-ggplot(aic_bic_long, aes(x = Model, y = Value, fill = Model)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.4, outlier.shape = TRUE, fill = "white", color = "black", linewidth = 0.2, alpha=0.6) +
  #geom_hline(data = mean_value_of_best_model,
  #           aes(yintercept = Line_Value),
  #          color = "red", linetype = "dashed", linewidth = 0.3,inherit.aes = FALSE) +
  facet_wrap(~ Criterion, scales = "free_y") +
  labs(title = "AIC and BIC plots of the fitted models at 5% low outliers",
       x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

p




# Save to PDF
ggsave("AIC_BIC_0.05_low_0.70_zero_degpd.pdf", plot = p, width = 6, height = 5)








###. == medium outliers with 0.05 ratio and 25% zeros

# Initialize storage matrices
results_matrix <- matrix(NA, nrow = 1000, ncol = 10)
colnames(results_matrix) <- c("Poisson_AIC", "Poisson_BIC",
                              "NB_AIC", "NB_BIC",
                              "ZINB_AIC", "ZINB_BIC",
                              "DEGPD_AIC", "DEGPD_BIC",
                              "ZIDEGPD_AIC", "ZIDEGPD_BIC")

set.seed(123)

for (i in 1:1000) {
  # Simulate data once
  sim_data <- simulate_one_rep(outlier_level = "medium", zero_inflation = 0.70, outlier_ratio = 0.05)
  
  # 1. Poisson Model
  poisson_model <- try(glm(Y ~ X1 + X2, family = poisson(link = "log"), data = sim_data), silent = TRUE)
  if (!inherits(poisson_model, "try-error")) {
    results_matrix[i, 1] <- AIC(poisson_model)
    results_matrix[i, 2] <- BIC(poisson_model)
  }
  
  # 2. Negative Binomial Model
  nb_model <- try(glm.nb(Y ~ X1 + X2, data = sim_data), silent = TRUE)
  if (!inherits(nb_model, "try-error")) {
    results_matrix[i, 3] <- AIC(nb_model)
    results_matrix[i, 4] <- BIC(nb_model)
  }
  
  # 3. ZINB Model
  zinb_model <- try(zeroinfl(Y ~ X1 + X2 | X1 + X2, data = sim_data, dist = "negbin"), silent = TRUE)
  if (!inherits(zinb_model, "try-error")) {
    results_matrix[i, 5] <- AIC(zinb_model)
    results_matrix[i, 6] <- BIC(zinb_model)
  }
  
  # 4. DEGPD Model
  fmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1)
  inits_degpd <- c(sqrt(sd(sim_data$Y)/2), sqrt(sd(sim_data$Y)/20), 0.3)
  degpd1 <- try(devgam(fmla_degpd1, data = sim_data, family = "degpd",
                       degpd.args = list(m = 1), trace = 0), silent = TRUE)
  if (!inherits(degpd1, "try-error")) {
    results_matrix[i, 7] <- AIC(degpd1)
    results_matrix[i, 8] <- BIC(degpd1)
  }
  
  # 5. ZIDEGPD Model
  zifmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1, logitpi = ~1)
  inits_zidegpd <- c(sd(sim_data$Y)/2, sd(sim_data$Y)/2, 0.2, mean(sim_data$Y == 0))
  zidegpd1 <- try(devgam(zifmla_degpd1, data = sim_data, family = "zidegpd",
                         degpd.args = list(m = 1), trace = 0, inits = inits_zidegpd), silent = TRUE)
  if (!inherits(zidegpd1, "try-error")) {
    results_matrix[i, 9] <- AIC(zidegpd1)
    results_matrix[i, 10] <- BIC(zidegpd1)
  }
  
  if (i %% 100 == 0) cat("Completed:", i, "\n")
}






# Compute means for each model
mean_aic_bic <- data.frame(
  Model = c("Poisson", "Negative Binomial", "ZINB", "DEGPD", "ZIDEGPD"),
  Mean_AIC = c(
    mean(results_matrix[, 1], na.rm = TRUE),
    mean(results_matrix[, 3], na.rm = TRUE),
    mean(results_matrix[, 5], na.rm = TRUE),
    mean(results_matrix[, 7], na.rm = TRUE),
    mean(results_matrix[, 9], na.rm = TRUE)
  ),
  Mean_BIC = c(
    mean(results_matrix[, 2], na.rm = TRUE),
    mean(results_matrix[, 4], na.rm = TRUE),
    mean(results_matrix[, 6], na.rm = TRUE),
    mean(results_matrix[, 8], na.rm = TRUE),
    mean(results_matrix[, 10], na.rm = TRUE)
  )
)

# Print summary table
print(mean_aic_bic, row.names = FALSE)

# Combine AIC/BIC matrices into a long-format dataframe
aic_bic_df <- data.frame(
  Model = rep(c( "NB", "ZINB", "DEGPD", "ZIDEGPD"), each = nrow(results_matrix)),
  AIC = c(#results_matrix[, 1],
    results_matrix[, 3],
    results_matrix[, 5],
    results_matrix[, 7],
    results_matrix[, 9]),
  BIC = c(#results_matrix[, 2],
    results_matrix[, 4],
    results_matrix[, 6],
    results_matrix[, 8],
    results_matrix[, 10])
)

# Convert to long format
aic_bic_long <- pivot_longer(aic_bic_df, cols = c(AIC, BIC), names_to = "Criterion", values_to = "Value")

# Reorder Model to show DEGPD and ZIDEGPD on the left
aic_bic_long$Model <- factor(aic_bic_long$Model,
                             levels = c("DEGPD", "ZIDEGPD", "NB", "ZINB")
)

# Compute mean of the best model (lowest minimum) per criterion
min_per_model <- aic_bic_long %>%
  group_by(Model, Criterion) %>%
  summarise(Min_Value = min(Value, na.rm = TRUE), .groups = "drop")

lowest_model_per_criterion <- min_per_model %>%
  group_by(Criterion) %>%
  filter(Min_Value == min(Min_Value)) %>%
  select(Criterion, Model)

mean_value_of_best_model <- aic_bic_long %>%
  inner_join(lowest_model_per_criterion, by = c("Criterion", "Model")) %>%
  group_by(Criterion) %>%
  summarise(Line_Value = mean(Value, na.rm = TRUE), .groups = "drop")

# ======= 1. Flask-Shaped Boxplot (varwidth) =======
ggplot(aic_bic_long, aes(x = Model, y = Value, fill = Model)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, varwidth = TRUE) +
  #geom_hline(data = mean_value_of_best_model,
  #aes(yintercept = Line_Value),
  #color = "red", linetype = "dashed", linewidth = 0.5,
  # inherit.aes = FALSE) +
  facet_wrap(~ Criterion, scales = "free_y") +
  labs(title = "Flask-Shaped Boxplots (Varwidth)\nDashed Red Line = Mean of Best Model",
       x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# ======= 2. Violin + Boxplot =======
p<-ggplot(aic_bic_long, aes(x = Model, y = Value, fill = Model)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.4, outlier.shape = TRUE, fill = "white", color = "black", linewidth = 0.2, alpha=0.6) +
  #geom_hline(data = mean_value_of_best_model,
  #           aes(yintercept = Line_Value),
  #          color = "red", linetype = "dashed", linewidth = 0.3,inherit.aes = FALSE) +
  facet_wrap(~ Criterion, scales = "free_y") +
  labs(title = "AIC and BIC plots of the fitted models at 5% medium outliers",
       x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

p




# Save to PDF
ggsave("AIC_BIC_0.05_medium_0.70_zero_degpd.pdf", plot = p, width = 6, height = 5)




###. == high outliers with 0.01 ratio and 25% zeros. ===

# Initialize storage matrices
results_matrix <- matrix(NA, nrow = 1000, ncol = 10)
colnames(results_matrix) <- c("Poisson_AIC", "Poisson_BIC",
                              "NB_AIC", "NB_BIC",
                              "ZINB_AIC", "ZINB_BIC",
                              "DEGPD_AIC", "DEGPD_BIC",
                              "ZIDEGPD_AIC", "ZIDEGPD_BIC")

set.seed(123)

for (i in 1:1000) {
  # Simulate data once
  sim_data <- simulate_one_rep(outlier_level = "high", zero_inflation = 0.70, outlier_ratio = 0.05)
  
  # 1. Poisson Model
  poisson_model <- try(glm(Y ~ X1 + X2, family = poisson(link = "log"), data = sim_data), silent = TRUE)
  if (!inherits(poisson_model, "try-error")) {
    results_matrix[i, 1] <- AIC(poisson_model)
    results_matrix[i, 2] <- BIC(poisson_model)
  }
  
  # 2. Negative Binomial Model
  nb_model <- try(glm.nb(Y ~ X1 + X2, data = sim_data), silent = TRUE)
  if (!inherits(nb_model, "try-error")) {
    results_matrix[i, 3] <- AIC(nb_model)
    results_matrix[i, 4] <- BIC(nb_model)
  }
  
  # 3. ZINB Model
  zinb_model <- try(zeroinfl(Y ~ X1 + X2 | X1 + X2, data = sim_data, dist = "negbin"), silent = TRUE)
  if (!inherits(zinb_model, "try-error")) {
    results_matrix[i, 5] <- AIC(zinb_model)
    results_matrix[i, 6] <- BIC(zinb_model)
  }
  
  # 4. DEGPD Model
  fmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1)
  inits_degpd <- c(sqrt(sd(sim_data$Y)/2), sqrt(sd(sim_data$Y)/20), 0.3)
  degpd1 <- try(devgam(fmla_degpd1, data = sim_data, family = "degpd",
                       degpd.args = list(m = 1), trace = 0), silent = TRUE)
  if (!inherits(degpd1, "try-error")) {
    results_matrix[i, 7] <- AIC(degpd1)
    results_matrix[i, 8] <- BIC(degpd1)
  }
  
  # 5. ZIDEGPD Model
  zifmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1, logitpi = ~1)
  inits_zidegpd <- c(sd(sim_data$Y)/2, sd(sim_data$Y)/2, 0.2, mean(sim_data$Y == 0))
  zidegpd1 <- try(devgam(zifmla_degpd1, data = sim_data, family = "zidegpd",
                         degpd.args = list(m = 1), trace = 0, inits = inits_zidegpd), silent = TRUE)
  if (!inherits(zidegpd1, "try-error")) {
    results_matrix[i, 9] <- AIC(zidegpd1)
    results_matrix[i, 10] <- BIC(zidegpd1)
  }
  
  if (i %% 100 == 0) cat("Completed:", i, "\n")
}



# Compute means for each model
mean_aic_bic <- data.frame(
  Model = c("Poisson", "Negative Binomial", "ZINB", "DEGPD", "ZIDEGPD"),
  Mean_AIC = c(
    mean(results_matrix[, 1], na.rm = TRUE),
    mean(results_matrix[, 3], na.rm = TRUE),
    mean(results_matrix[, 5], na.rm = TRUE),
    mean(results_matrix[, 7], na.rm = TRUE),
    mean(results_matrix[, 9], na.rm = TRUE)
  ),
  Mean_BIC = c(
    mean(results_matrix[, 2], na.rm = TRUE),
    mean(results_matrix[, 4], na.rm = TRUE),
    mean(results_matrix[, 6], na.rm = TRUE),
    mean(results_matrix[, 8], na.rm = TRUE),
    mean(results_matrix[, 10], na.rm = TRUE)
  )
)

# Print summary table
print(mean_aic_bic, row.names = FALSE)

# Combine AIC/BIC matrices into a long-format dataframe
aic_bic_df <- data.frame(
  Model = rep(c( "NB", "ZINB", "DEGPD", "ZIDEGPD"), each = nrow(results_matrix)),
  AIC = c(#results_matrix[, 1],
    results_matrix[, 3],
    results_matrix[, 5],
    results_matrix[, 7],
    results_matrix[, 9]),
  BIC = c(#results_matrix[, 2],
    results_matrix[, 4],
    results_matrix[, 6],
    results_matrix[, 8],
    results_matrix[, 10])
)

# Convert to long format
aic_bic_long <- pivot_longer(aic_bic_df, cols = c(AIC, BIC), names_to = "Criterion", values_to = "Value")

# Reorder Model to show DEGPD and ZIDEGPD on the left
aic_bic_long$Model <- factor(aic_bic_long$Model,
                             levels = c("DEGPD", "ZIDEGPD", "NB", "ZINB")
)

# Compute mean of the best model (lowest minimum) per criterion
min_per_model <- aic_bic_long %>%
  group_by(Model, Criterion) %>%
  summarise(Min_Value = min(Value, na.rm = TRUE), .groups = "drop")

lowest_model_per_criterion <- min_per_model %>%
  group_by(Criterion) %>%
  filter(Min_Value == min(Min_Value)) %>%
  select(Criterion, Model)

mean_value_of_best_model <- aic_bic_long %>%
  inner_join(lowest_model_per_criterion, by = c("Criterion", "Model")) %>%
  group_by(Criterion) %>%
  summarise(Line_Value = mean(Value, na.rm = TRUE), .groups = "drop")

# ======= 1. Flask-Shaped Boxplot (varwidth) =======
ggplot(aic_bic_long, aes(x = Model, y = Value, fill = Model)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, varwidth = TRUE) +
  #geom_hline(data = mean_value_of_best_model,
  #aes(yintercept = Line_Value),
  #color = "red", linetype = "dashed", linewidth = 0.5,
  # inherit.aes = FALSE) +
  facet_wrap(~ Criterion, scales = "free_y") +
  labs(title = "Flask-Shaped Boxplots (Varwidth)\nDashed Red Line = Mean of Best Model",
       x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# ======= 2. Violin + Boxplot =======
p<-ggplot(aic_bic_long, aes(x = Model, y = Value, fill = Model)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.4, outlier.shape = TRUE, fill = "white", color = "black", linewidth = 0.2, alpha=0.6) +
  #geom_hline(data = mean_value_of_best_model,
  #           aes(yintercept = Line_Value),
  #          color = "red", linetype = "dashed", linewidth = 0.3,inherit.aes = FALSE) +
  facet_wrap(~ Criterion, scales = "free_y") +
  labs(title = "AIC and BIC plots of the fitted models at 5% high outliers",
       x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

p




# Save to PDF
ggsave("AIC_BIC_0.05_high_0.70_zero_degpd.pdf", plot = p, width = 6, height = 5)










##===================================================================================================================================================================



###. == Low outliers with 10% ratio and 25% zeros

# Initialize storage matrices
results_matrix <- matrix(NA, nrow = 1000, ncol = 10)
colnames(results_matrix) <- c("Poisson_AIC", "Poisson_BIC",
                              "NB_AIC", "NB_BIC",
                              "ZINB_AIC", "ZINB_BIC",
                              "DEGPD_AIC", "DEGPD_BIC",
                              "ZIDEGPD_AIC", "ZIDEGPD_BIC")

set.seed(123)

for (i in 1:1000) {
  # Simulate data once
  sim_data <- simulate_one_rep(outlier_level = "low", zero_inflation = 0.70, outlier_ratio = 0.10)
  
  # 1. Poisson Model
  poisson_model <- try(glm(Y ~ X1 + X2, family = poisson(link = "log"), data = sim_data), silent = TRUE)
  if (!inherits(poisson_model, "try-error")) {
    results_matrix[i, 1] <- AIC(poisson_model)
    results_matrix[i, 2] <- BIC(poisson_model)
  }
  
  # 2. Negative Binomial Model
  nb_model <- try(glm.nb(Y ~ X1 + X2, data = sim_data), silent = TRUE)
  if (!inherits(nb_model, "try-error")) {
    results_matrix[i, 3] <- AIC(nb_model)
    results_matrix[i, 4] <- BIC(nb_model)
  }
  
  # 3. ZINB Model
  zinb_model <- try(zeroinfl(Y ~ X1 + X2 | X1 + X2, data = sim_data, dist = "negbin"), silent = TRUE)
  if (!inherits(zinb_model, "try-error")) {
    results_matrix[i, 5] <- AIC(zinb_model)
    results_matrix[i, 6] <- BIC(zinb_model)
  }
  
  # 4. DEGPD Model
  fmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1)
  inits_degpd <- c(sqrt(sd(sim_data$Y)/2), sqrt(sd(sim_data$Y)/20), 0.3)
  degpd1 <- try(devgam(fmla_degpd1, data = sim_data, family = "degpd",
                       degpd.args = list(m = 1), trace = 0), silent = TRUE)
  if (!inherits(degpd1, "try-error")) {
    results_matrix[i, 7] <- AIC(degpd1)
    results_matrix[i, 8] <- BIC(degpd1)
  }
  
  # 5. ZIDEGPD Model
  zifmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1, logitpi = ~1)
  inits_zidegpd <- c(sd(sim_data$Y)/2, sd(sim_data$Y)/2, 0.2, mean(sim_data$Y == 0))
  zidegpd1 <- try(devgam(zifmla_degpd1, data = sim_data, family = "zidegpd",
                         degpd.args = list(m = 1), trace = 0, inits = inits_zidegpd), silent = TRUE)
  if (!inherits(zidegpd1, "try-error")) {
    results_matrix[i, 9] <- AIC(zidegpd1)
    results_matrix[i, 10] <- BIC(zidegpd1)
  }
  
  if (i %% 100 == 0) cat("Completed:", i, "\n")
}




# Compute means for each model
mean_aic_bic <- data.frame(
  Model = c("Poisson", "Negative Binomial", "ZINB", "DEGPD", "ZIDEGPD"),
  Mean_AIC = c(
    mean(results_matrix[, 1], na.rm = TRUE),
    mean(results_matrix[, 3], na.rm = TRUE),
    mean(results_matrix[, 5], na.rm = TRUE),
    mean(results_matrix[, 7], na.rm = TRUE),
    mean(results_matrix[, 9], na.rm = TRUE)
  ),
  Mean_BIC = c(
    mean(results_matrix[, 2], na.rm = TRUE),
    mean(results_matrix[, 4], na.rm = TRUE),
    mean(results_matrix[, 6], na.rm = TRUE),
    mean(results_matrix[, 8], na.rm = TRUE),
    mean(results_matrix[, 10], na.rm = TRUE)
  )
)

# Print summary table
print(mean_aic_bic, row.names = FALSE)

# Combine AIC/BIC matrices into a long-format dataframe
aic_bic_df <- data.frame(
  Model = rep(c( "NB", "ZINB", "DEGPD", "ZIDEGPD"), each = nrow(results_matrix)),
  AIC = c(#results_matrix[, 1],
    results_matrix[, 3],
    results_matrix[, 5],
    results_matrix[, 7],
    results_matrix[, 9]),
  BIC = c(#results_matrix[, 2],
    results_matrix[, 4],
    results_matrix[, 6],
    results_matrix[, 8],
    results_matrix[, 10])
)

# Convert to long format
aic_bic_long <- pivot_longer(aic_bic_df, cols = c(AIC, BIC), names_to = "Criterion", values_to = "Value")

# Reorder Model to show DEGPD and ZIDEGPD on the left
aic_bic_long$Model <- factor(aic_bic_long$Model,
                             levels = c("DEGPD", "ZIDEGPD", "NB", "ZINB")
)

# Compute mean of the best model (lowest minimum) per criterion
min_per_model <- aic_bic_long %>%
  group_by(Model, Criterion) %>%
  summarise(Min_Value = min(Value, na.rm = TRUE), .groups = "drop")

lowest_model_per_criterion <- min_per_model %>%
  group_by(Criterion) %>%
  filter(Min_Value == min(Min_Value)) %>%
  select(Criterion, Model)

mean_value_of_best_model <- aic_bic_long %>%
  inner_join(lowest_model_per_criterion, by = c("Criterion", "Model")) %>%
  group_by(Criterion) %>%
  summarise(Line_Value = mean(Value, na.rm = TRUE), .groups = "drop")

# ======= 1. Flask-Shaped Boxplot (varwidth) =======
ggplot(aic_bic_long, aes(x = Model, y = Value, fill = Model)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, varwidth = TRUE) +
  #geom_hline(data = mean_value_of_best_model,
  #aes(yintercept = Line_Value),
  #color = "red", linetype = "dashed", linewidth = 0.5,
  # inherit.aes = FALSE) +
  facet_wrap(~ Criterion, scales = "free_y") +
  labs(title = "Flask-Shaped Boxplots (Varwidth)\nDashed Red Line = Mean of Best Model",
       x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# ======= 2. Violin + Boxplot =======
p<-ggplot(aic_bic_long, aes(x = Model, y = Value, fill = Model)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.4, outlier.shape = TRUE, fill = "white", color = "black", linewidth = 0.2, alpha=0.6) +
  #geom_hline(data = mean_value_of_best_model,
  #           aes(yintercept = Line_Value),
  #          color = "red", linetype = "dashed", linewidth = 0.3,inherit.aes = FALSE) +
  facet_wrap(~ Criterion, scales = "free_y") +
  labs(title = "AIC and BIC plots of the fitted models at 10% low outliers",
       x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

p




# Save to PDF
ggsave("AIC_BIC_0.10_low_0.70_zero_degpd.pdf", plot = p, width = 6, height = 5)








###. == medium outliers with 0.10 ratio and 25% zeros

# Initialize storage matrices
results_matrix <- matrix(NA, nrow = 1000, ncol = 10)
colnames(results_matrix) <- c("Poisson_AIC", "Poisson_BIC",
                              "NB_AIC", "NB_BIC",
                              "ZINB_AIC", "ZINB_BIC",
                              "DEGPD_AIC", "DEGPD_BIC",
                              "ZIDEGPD_AIC", "ZIDEGPD_BIC")

set.seed(123)

for (i in 1:1000) {
  # Simulate data once
  sim_data <- simulate_one_rep(outlier_level = "medium", zero_inflation = 0.70, outlier_ratio = 0.10)
  
  # 1. Poisson Model
  poisson_model <- try(glm(Y ~ X1 + X2, family = poisson(link = "log"), data = sim_data), silent = TRUE)
  if (!inherits(poisson_model, "try-error")) {
    results_matrix[i, 1] <- AIC(poisson_model)
    results_matrix[i, 2] <- BIC(poisson_model)
  }
  
  # 2. Negative Binomial Model
  nb_model <- try(glm.nb(Y ~ X1 + X2, data = sim_data), silent = TRUE)
  if (!inherits(nb_model, "try-error")) {
    results_matrix[i, 3] <- AIC(nb_model)
    results_matrix[i, 4] <- BIC(nb_model)
  }
  
  # 3. ZINB Model
  zinb_model <- try(zeroinfl(Y ~ X1 + X2 | X1 + X2, data = sim_data, dist = "negbin"), silent = TRUE)
  if (!inherits(zinb_model, "try-error")) {
    results_matrix[i, 5] <- AIC(zinb_model)
    results_matrix[i, 6] <- BIC(zinb_model)
  }
  
  # 4. DEGPD Model
  fmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1)
  inits_degpd <- c(sqrt(sd(sim_data$Y)/2), sqrt(sd(sim_data$Y)/20), 0.3)
  degpd1 <- try(devgam(fmla_degpd1, data = sim_data, family = "degpd",
                       degpd.args = list(m = 1), trace = 0), silent = TRUE)
  if (!inherits(degpd1, "try-error")) {
    results_matrix[i, 7] <- AIC(degpd1)
    results_matrix[i, 8] <- BIC(degpd1)
  }
  
  # 5. ZIDEGPD Model
  zifmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1, logitpi = ~1)
  inits_zidegpd <- c(sd(sim_data$Y)/2, sd(sim_data$Y)/2, 0.2, mean(sim_data$Y == 0))
  zidegpd1 <- try(devgam(zifmla_degpd1, data = sim_data, family = "zidegpd",
                         degpd.args = list(m = 1), trace = 0, inits = inits_zidegpd), silent = TRUE)
  if (!inherits(zidegpd1, "try-error")) {
    results_matrix[i, 9] <- AIC(zidegpd1)
    results_matrix[i, 10] <- BIC(zidegpd1)
  }
  
  if (i %% 100 == 0) cat("Completed:", i, "\n")
}






# Compute means for each model
mean_aic_bic <- data.frame(
  Model = c("Poisson", "Negative Binomial", "ZINB", "DEGPD", "ZIDEGPD"),
  Mean_AIC = c(
    mean(results_matrix[, 1], na.rm = TRUE),
    mean(results_matrix[, 3], na.rm = TRUE),
    mean(results_matrix[, 5], na.rm = TRUE),
    mean(results_matrix[, 7], na.rm = TRUE),
    mean(results_matrix[, 9], na.rm = TRUE)
  ),
  Mean_BIC = c(
    mean(results_matrix[, 2], na.rm = TRUE),
    mean(results_matrix[, 4], na.rm = TRUE),
    mean(results_matrix[, 6], na.rm = TRUE),
    mean(results_matrix[, 8], na.rm = TRUE),
    mean(results_matrix[, 10], na.rm = TRUE)
  )
)

# Print summary table
print(mean_aic_bic, row.names = FALSE)

# Combine AIC/BIC matrices into a long-format dataframe
aic_bic_df <- data.frame(
  Model = rep(c( "NB", "ZINB", "DEGPD", "ZIDEGPD"), each = nrow(results_matrix)),
  AIC = c(#results_matrix[, 1],
    results_matrix[, 3],
    results_matrix[, 5],
    results_matrix[, 7],
    results_matrix[, 9]),
  BIC = c(#results_matrix[, 2],
    results_matrix[, 4],
    results_matrix[, 6],
    results_matrix[, 8],
    results_matrix[, 10])
)

# Convert to long format
aic_bic_long <- pivot_longer(aic_bic_df, cols = c(AIC, BIC), names_to = "Criterion", values_to = "Value")

# Reorder Model to show DEGPD and ZIDEGPD on the left
aic_bic_long$Model <- factor(aic_bic_long$Model,
                             levels = c("DEGPD", "ZIDEGPD", "NB", "ZINB")
)

# Compute mean of the best model (lowest minimum) per criterion
min_per_model <- aic_bic_long %>%
  group_by(Model, Criterion) %>%
  summarise(Min_Value = min(Value, na.rm = TRUE), .groups = "drop")

lowest_model_per_criterion <- min_per_model %>%
  group_by(Criterion) %>%
  filter(Min_Value == min(Min_Value)) %>%
  select(Criterion, Model)

mean_value_of_best_model <- aic_bic_long %>%
  inner_join(lowest_model_per_criterion, by = c("Criterion", "Model")) %>%
  group_by(Criterion) %>%
  summarise(Line_Value = mean(Value, na.rm = TRUE), .groups = "drop")

# ======= 1. Flask-Shaped Boxplot (varwidth) =======
ggplot(aic_bic_long, aes(x = Model, y = Value, fill = Model)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, varwidth = TRUE) +
  #geom_hline(data = mean_value_of_best_model,
  #aes(yintercept = Line_Value),
  #color = "red", linetype = "dashed", linewidth = 0.5,
  # inherit.aes = FALSE) +
  facet_wrap(~ Criterion, scales = "free_y") +
  labs(title = "Flask-Shaped Boxplots (Varwidth)\nDashed Red Line = Mean of Best Model",
       x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# ======= 2. Violin + Boxplot =======
p<-ggplot(aic_bic_long, aes(x = Model, y = Value, fill = Model)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.4, outlier.shape = TRUE, fill = "white", color = "black", linewidth = 0.2, alpha=0.6) +
  #geom_hline(data = mean_value_of_best_model,
  #           aes(yintercept = Line_Value),
  #          color = "red", linetype = "dashed", linewidth = 0.3,inherit.aes = FALSE) +
  facet_wrap(~ Criterion, scales = "free_y") +
  labs(title = "AIC and BIC plots of the fitted models at 10% medium outliers",
       x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

p




# Save to PDF
ggsave("AIC_BIC_0.10_medium_0.70_zero_degpd.pdf", plot = p, width = 6, height = 5)




###. == high outliers with 0.10 ratio and 25% zeros. ===

# Initialize storage matrices
results_matrix <- matrix(NA, nrow = 1000, ncol = 10)
colnames(results_matrix) <- c("Poisson_AIC", "Poisson_BIC",
                              "NB_AIC", "NB_BIC",
                              "ZINB_AIC", "ZINB_BIC",
                              "DEGPD_AIC", "DEGPD_BIC",
                              "ZIDEGPD_AIC", "ZIDEGPD_BIC")

set.seed(123)

for (i in 1:1000) {
  # Simulate data once
  sim_data <- simulate_one_rep(outlier_level = "high", zero_inflation = 0.70, outlier_ratio = 0.10)
  
  # 1. Poisson Model
  poisson_model <- try(glm(Y ~ X1 + X2, family = poisson(link = "log"), data = sim_data), silent = TRUE)
  if (!inherits(poisson_model, "try-error")) {
    results_matrix[i, 1] <- AIC(poisson_model)
    results_matrix[i, 2] <- BIC(poisson_model)
  }
  
  # 2. Negative Binomial Model
  nb_model <- try(glm.nb(Y ~ X1 + X2, data = sim_data), silent = TRUE)
  if (!inherits(nb_model, "try-error")) {
    results_matrix[i, 3] <- AIC(nb_model)
    results_matrix[i, 4] <- BIC(nb_model)
  }
  
  # 3. ZINB Model
  zinb_model <- try(zeroinfl(Y ~ X1 + X2 | X1 + X2, data = sim_data, dist = "negbin"), silent = TRUE)
  if (!inherits(zinb_model, "try-error")) {
    results_matrix[i, 5] <- AIC(zinb_model)
    results_matrix[i, 6] <- BIC(zinb_model)
  }
  
  # 4. DEGPD Model
  fmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1)
  # inits_degpd <- c(sqrt(sd(sim_data$Y)/2), sqrt(sd(sim_data$Y)/20), 0.3)
  degpd1 <- try(devgam(fmla_degpd1, data = sim_data, family = "degpd",
                       degpd.args = list(m = 1), trace = 0), silent = TRUE)
  if (!inherits(degpd1, "try-error")) {
    results_matrix[i, 7] <- AIC(degpd1)
    results_matrix[i, 8] <- BIC(degpd1)
  }
  
  # 5. ZIDEGPD Model
  zifmla_degpd1 <- list(lsigma = Y ~ X1 + X2, lxi = ~1, lkappa = ~1, logitpi = ~1)
  inits_zidegpd <- c(sd(sim_data$Y)/2, sd(sim_data$Y)/2, 0.2, mean(sim_data$Y == 0))
  zidegpd1 <- try(devgam(zifmla_degpd1, data = sim_data, family = "zidegpd",
                         degpd.args = list(m = 1), trace = 0, inits = inits_zidegpd), silent = TRUE)
  if (!inherits(zidegpd1, "try-error")) {
    results_matrix[i, 9] <- AIC(zidegpd1)
    results_matrix[i, 10] <- BIC(zidegpd1)
  }
  
  if (i %% 100 == 0) cat("Completed:", i, "\n")
}



# Compute means for each model
mean_aic_bic <- data.frame(
  Model = c("Poisson", "Negative Binomial", "ZINB", "DEGPD", "ZIDEGPD"),
  Mean_AIC = c(
    mean(results_matrix[, 1], na.rm = TRUE),
    mean(results_matrix[, 3], na.rm = TRUE),
    mean(results_matrix[, 5], na.rm = TRUE),
    mean(results_matrix[, 7], na.rm = TRUE),
    mean(results_matrix[, 9], na.rm = TRUE)
  ),
  Mean_BIC = c(
    mean(results_matrix[, 2], na.rm = TRUE),
    mean(results_matrix[, 4], na.rm = TRUE),
    mean(results_matrix[, 6], na.rm = TRUE),
    mean(results_matrix[, 8], na.rm = TRUE),
    mean(results_matrix[, 10], na.rm = TRUE)
  )
)

# Print summary table
print(mean_aic_bic, row.names = FALSE)

# Combine AIC/BIC matrices into a long-format dataframe
aic_bic_df <- data.frame(
  Model = rep(c( "NB", "ZINB", "DEGPD", "ZIDEGPD"), each = nrow(results_matrix)),
  AIC = c(#results_matrix[, 1],
    results_matrix[, 3],
    results_matrix[, 5],
    results_matrix[, 7],
    results_matrix[, 9]),
  BIC = c(#results_matrix[, 2],
    results_matrix[, 4],
    results_matrix[, 6],
    results_matrix[, 8],
    results_matrix[, 10])
)

# Convert to long format
aic_bic_long <- pivot_longer(aic_bic_df, cols = c(AIC, BIC), names_to = "Criterion", values_to = "Value")

# Reorder Model to show DEGPD and ZIDEGPD on the left
aic_bic_long$Model <- factor(aic_bic_long$Model,
                             levels = c("DEGPD", "ZIDEGPD", "NB", "ZINB")
)

# Compute mean of the best model (lowest minimum) per criterion
min_per_model <- aic_bic_long %>%
  group_by(Model, Criterion) %>%
  summarise(Min_Value = min(Value, na.rm = TRUE), .groups = "drop")

lowest_model_per_criterion <- min_per_model %>%
  group_by(Criterion) %>%
  filter(Min_Value == min(Min_Value)) %>%
  select(Criterion, Model)

mean_value_of_best_model <- aic_bic_long %>%
  inner_join(lowest_model_per_criterion, by = c("Criterion", "Model")) %>%
  group_by(Criterion) %>%
  summarise(Line_Value = mean(Value, na.rm = TRUE), .groups = "drop")

# ======= 1. Flask-Shaped Boxplot (varwidth) =======
ggplot(aic_bic_long, aes(x = Model, y = Value, fill = Model)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, varwidth = TRUE) +
  #geom_hline(data = mean_value_of_best_model,
  #aes(yintercept = Line_Value),
  #color = "red", linetype = "dashed", linewidth = 0.5,
  # inherit.aes = FALSE) +
  facet_wrap(~ Criterion, scales = "free_y") +
  labs(title = "Flask-Shaped Boxplots (Varwidth)\nDashed Red Line = Mean of Best Model",
       x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# ======= 2. Violin + Boxplot =======
p<-ggplot(aic_bic_long, aes(x = Model, y = Value, fill = Model)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.4, outlier.shape = TRUE, fill = "white", color = "black", linewidth = 0.2, alpha=0.6) +
  #geom_hline(data = mean_value_of_best_model,
  #           aes(yintercept = Line_Value),
  #          color = "red", linetype = "dashed", linewidth = 0.3,inherit.aes = FALSE) +
  facet_wrap(~ Criterion, scales = "free_y") +
  labs(title = "AIC and BIC plots of the fitted models at 10% high outliers",
       x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

p




# Save to PDF
ggsave("AIC_BIC_0.10_high_0.70_zero_degpd.pdf", plot = p, width = 6, height = 5)


