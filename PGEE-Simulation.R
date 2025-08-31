# ===================================================================
# 1. Load Required Libraries
# ===================================================================
# The PGEE package is used for the penalized GEE model.
# The geepack package is for the standard and oracle GEE models.
# The MASS package is used for generating multivariate normal data.
# The dplyr package is used for easy data manipulation and summarization.
# ===================================================================
library(PGEE) 
library(geepack) 
library(MASS) 
library(dplyr) 
library(ggplot2)
library(gridExtra)
library(reshape2)
library(tidyr)


# ===================================================================
# 2. Define Simulation Parameters
# ===================================================================
set.seed(2025) # for reproducibility

nsim <- 100 # Number of simulation repetitions
p <- 150 # Total number of covariates
n <- 200 # Number of subjects (sample size)
m <- 4 # Number of repeated observations per subject

# --- True regression coefficients ---
# Define the true non-zero coefficients. The rest will be zero.
beta_true_cont <- c(2.5, -2.0, 1.5, -1.0, 0.5)

# Create the full vector of true coefficients (p-dimensional)
beta_true_full <- c(beta_true_cont, rep(0, p - length(beta_true_cont)))

# Identify the positions of the true non-zero coefficients for the oracle model
signal_idx <- which(beta_true_full != 0) # Correctly identify non-zero signals

# --- Define true correlation values to test ---
rho_values <- c(0.4, 0.8)

# --- Candidate values for the PGEE tuning parameter ---
lambda.vec <- seq(0.05, 0.2, length.out = 20)



# ===================================================================
# 3. Define the Working Correlation Structures to Test
# ===================================================================
# This vector controls the new loop. We will fit models using each of these.
working_correlations <- c("independence", "exchangeable", "ar1")







# ===================================================================
# 4. Performance Evaluation Helper Function
# ===================================================================
evaluate_performance <- function(est_betas, true_betas, threshold = 1e-3) {
  # Identify the true non-zero coefficient positions
  signal_idx <- which(true_betas != 0)
  
  # Calculate Mean Squared Error (MSE)
  mse_val <- mean((est_betas - true_betas)^2)
  
  # Identify the positions of coefficients selected by the model
  selected <- which(abs(est_betas) > threshold)
  
  # Calculate True Positives (TP) and False Positives (FP)
  tp_val <- sum(selected %in% signal_idx)
  fp_val <- sum(!(selected %in% signal_idx))
  n_selected <- length(selected)
  n_true_signals <- length(signal_idx)
  
  # U (Underfitting): The model selected fewer variables than the true number of signals
  u_val <- as.numeric(n_selected < n_true_signals)
  
  # O (Overfitting): The model selected more variables than the true number of signals
  o_val <- as.numeric(n_selected > n_true_signals)
  
  # EXACT: The set of selected variables is identical to the true set
  exact_val <- as.numeric(setequal(selected, signal_idx))
  
  # Return a list containing all performance metrics
  return(list(MSE = mse_val, U = u_val, O = o_val, EXACT = exact_val, TP = tp_val, FP = fp_val))
}

# Example
#Beta_hat=c(2,1,1); Beta_real= c(2.5,0,1)
#evaluate_performance(Beta_hat,Beta_real)  






# ===================================================================
# 5. Main Simulation Loop
# ===================================================================

# --- Initialize lists to store results from all simulations ---
results_list_table1 <- list()
pgee_results_for_table2 <- list() # New list for Table 2 data
lambda_record <- c()



# --- Outermost loop for each true rho value ---
for (rho_val in rho_values) {
  
  cat(sprintf("\n\n#################### STARTING SIMULATIONS FOR RHO = %.1f ####################\n", rho_val))
  
  # --- Set the true correlation for the current set of simulations ---
  rho_eps_ar1 <- rho_val
  
  # --- Start the main loop for each simulation run ---
  for (sim_num in 1:nsim) {
    
    cat(sprintf("\n\n=============== Starting Simulation %d / %d (for rho = %.1f) ===============\n", sim_num, nsim, rho_eps_ar1))
    
    # --- Generate one simulated dataset ---
    
    # --- MODIFIED: Generate Covariates with Structure as specified ---
    # Correlation coefficient for the AR(1) structure of continuous covariates
    rho_x <- 0.5
    # Define the AR(1) covariance matrix for the p-1 continuous covariates
    cov_x <- outer(1:(p-1), 1:(p-1), function(i, j) rho_x^abs(i-j))
    
    # Generate the p-1 correlated continuous covariates for n subjects 
    X_continuous_wide <- mvrnorm(n = n, mu = rep(0, p - 1), Sigma = cov_x)
    
    # Generate the one binary covariate for n subjects 
    X_binary_wide <- rbinom(n = n, size = 1, prob = 0.5)
    
    # Combine into the final n x p time-invariant covariate matrix X
    # The first column is binary, the rest are correlated continuous
    X <- cbind(X_binary_wide, X_continuous_wide)
    # --- End of Covariate Generation Modification ---

    # Initialize the response matrix Y
    y_matrix <- matrix(0, nrow = n, ncol = m)
    
    # Define the true AR(1) correlation structure for the error term
    R_eps <- outer(1:m, 1:m, function(i, j) rho_eps_ar1^abs(i - j))
    
    # Generate data for each subject
    for (i in 1:n) {
      lin_pred <- X[i, ] %*% beta_true_full
      errors <- mvrnorm(n = 1, mu = rep(0, m), Sigma = R_eps)
      y_matrix[i, ] <- lin_pred + errors
    }
    
    # Convert data to the long format required by GEE functions
    colnames(X) <- paste0("X", 1:p)
    current_data <- data.frame(
      id = rep(1:n, each = m),
      y = as.vector(t(y_matrix)),
      X[rep(1:n, each = m), ]
    )
    
    # --- Define model formulas ---
    # Full model formula: y ~ X1 + X2 + ... + Xp
    formula_obj <- as.formula(paste("y ~", paste(paste0("X", 1:p), collapse = " + ")))
    # Oracle model formula: y ~ X1 + ... + X5 (only true signals)
    oracle_formula <- as.formula(paste("y ~", paste(paste0("X", signal_idx), collapse = " + ")))
    
    # --- Find the best lambda for PGEE using cross-validation (once per dataset) ---
    # As per the documentation, CVfit uses a working independence assumption.
    cv_result <- CVfit(
      formula = formula_obj,
      id = id,
      data = current_data,
      family = gaussian(),
      lambda.vec = lambda.vec,
      fold = 5
    )
    best_lambda <- cv_result$lam.opt
    lambda_record <- c(lambda_record, best_lambda)
    cat(sprintf("  - Cross-validation complete. Best lambda found: %.4f\n", best_lambda))
    
    
    # --- Loop over each specified working correlation structure ---
    for (w_cor in working_correlations) {
      
      cat(sprintf("  - Fitting models with '%s' working correlation...\n", w_cor))
      
      # --- Model 1: Standard GEE ---
      gee_fit <- geeglm(
        formula = formula_obj,
        id = id,
        data = current_data,
        family = gaussian(),
        corstr = w_cor # Use the current working correlation
      )
      # Extract coefficients (excluding the intercept)
      gee_betas <- coef(gee_fit)[-1]
      #print(gee_betas) 
      
      
      # --- Model 2: Oracle GEE ---
      oracle_fit <- geeglm(
        formula = oracle_formula,
        id = id,
        data = current_data,
        family = gaussian(),
        corstr = w_cor # Use the current working correlation
      )
      # Place the estimated coefficients back into a full p-dimensional vector
      oracle_betas_short <- coef(oracle_fit)[-1]
      oracle_betas <- numeric(p)
      oracle_betas[signal_idx] <- oracle_betas_short
      #print(oracle_betas)
      
      
      
      # --- Model 3: Penalized GEE ---
      # The PGEE package uses "AR-1" instead of "ar1".
      pgee_corstr <- if (w_cor == "ar1") "AR-1" else w_cor
      
      pgee_fit <- PGEE(
        formula = formula_obj,
        id = id,
        data = current_data,
        family = gaussian(),
        corstr = pgee_corstr, # Use the current working correlation
        lambda = best_lambda,  # Use the lambda from CV
        pindex = 1             # Do not penalize the intercept
      )
      pgee_betas <- coef(pgee_fit)[-1]
      #print(pgee_betas)
      
      
      
      # --- Evaluate and store results for the current working correlation ---
      perf_gee <- evaluate_performance(gee_betas, beta_true_full)
      perf_oracle <- evaluate_performance(oracle_betas, beta_true_full)
      perf_pgee <- evaluate_performance(pgee_betas, beta_true_full)
      
      # Add results to the list for Table 1
      results_list_table1[[length(results_list_table1) + 1]] <- data.frame(sim = sim_num, Method = "GEE", WorkingCorr = w_cor, TrueRho = rho_val, perf_gee)
      results_list_table1[[length(results_list_table1) + 1]] <- data.frame(sim = sim_num, Method = "Oracle", WorkingCorr = w_cor, TrueRho = rho_val, perf_oracle)
      results_list_table1[[length(results_list_table1) + 1]] <- data.frame(sim = sim_num, Method = "PGEE", WorkingCorr = w_cor, TrueRho = rho_val, perf_pgee)
      
      # --- NEW: Store PGEE results for Table 2 ---
      pgee_summary <- summary(pgee_fit)
      # summary's coefficient table includes the intercept, so signal indices are shifted by 1
      pgee_coeffs_summary <- pgee_summary$coefficients[(signal_idx + 1), ]
      
      pgee_results_for_table2[[length(pgee_results_for_table2) + 1]] <- data.frame(
        sim = sim_num, TrueRho = rho_val, WorkingCorr = w_cor, beta_idx = signal_idx,
        Estimate = pgee_coeffs_summary[, "Estimate"],
        RobustSE = pgee_coeffs_summary[, "Robust S.E."]
      )
      
    } # End of the working correlation loop
    
  } # End of the main simulation loop
  
} # End of the new outermost rho loop






# ===================================================================
# 6. Summarize and Print Final Results (Table 1)
# ===================================================================

cat("\n\n=============== All simulations complete. Summarizing results. ===============\n")

# Combine all data frames in the list into a single data frame
final_results_df_table1 <- do.call(rbind, results_list_table1)

# Group by Method, WorkingCorr, and TrueRho, then calculate the average of each metric
summary_table1 <- final_results_df_table1 %>%
  group_by(Method, WorkingCorr, TrueRho) %>%
  summarise(
    MSE = mean(MSE),
    U = mean(U),
    O = mean(O),
    EXACT = mean(EXACT),
    TP = mean(TP),
    FP = mean(FP),
    .groups = 'drop' # Recommended to avoid issues with further operations
  )

# --- Print the final summary table ---
cat("------------------------------------------------------------------\n")
cat(sprintf("Average Performance Over %d Simulations (Table 1):\n", nsim))
cat("------------------------------------------------------------------\n")
print(summary_table1)
cat("------------------------------------------------------------------\n")


saveRDS(summary_table1, file = "C:/Users/Lenovo/Desktop/summary_table1.rds")




# ===================================================================
# 7.Summarize and Print Final Results (Table 2 Format)
# ===================================================================
table2_df <- do.call(rbind, pgee_results_for_table2)
table2_summary <- table2_df %>%
  group_by(TrueRho, WorkingCorr, beta_idx) %>%
  summarise(
    Bias = abs(mean(Estimate) - beta_true_full[first(beta_idx)]),
    SD1 = mean(RobustSE),
    SD2 = sd(Estimate, na.rm = TRUE),
    CP = mean((Estimate - 1.96 * RobustSE <= beta_true_full[first(beta_idx)]) & 
                (Estimate + 1.96 * RobustSE >= beta_true_full[first(beta_idx)])),
    .groups = 'drop'
  )

cat("\n\n------------------------------------------------------------------\n")
cat("PGEE Detailed Performance Metrics (Table 2 Format):\n")
cat("------------------------------------------------------------------\n")


# --- ADD THIS LINE TO SAVE TABLE 2 ---
saveRDS(table2_summary, file = "C:/Users/Lenovo/Desktop/table2_summary.rds")




# Loop through each rho value to print its own table
for (rho in unique(table2_summary$TrueRho)) {
  cat(sprintf("\n----- True Rho = %.1f -----\n", rho))
  
  rho_subset <- table2_summary %>% filter(TrueRho == rho)
  
  # Create the header for the table
  header <- sprintf("%-15s %-8s", "Correlation", "Metric")
  for(b_idx in signal_idx) {
    header <- paste(header, sprintf("%-10s", paste0("\u03B2_", b_idx)))
  }
  cat(header, "\n")
  
  # Loop through each working correlation to print its rows
  for (corr in c("independence", "exchangeable", "ar1")) {
    # Check if data for this correlation exists
    if(!(corr %in% rho_subset$WorkingCorr)) next
    
    corr_subset <- rho_subset %>% filter(WorkingCorr == corr)
    
    # Prepare the rows for each metric
    bias_row <- sprintf("%-15s %-8s", corr, "Bias")
    sd1_row  <- sprintf("%-15s %-8s", "", "SD1")
    sd2_row  <- sprintf("%-15s %-8s", "", "SD2")
    cp_row   <- sprintf("%-15s %-8s", "", "CP")
    
    # Populate the rows with data for each beta
    for(b_idx in signal_idx) {
      vals <- corr_subset %>% filter(beta_idx == b_idx)
      if(nrow(vals) > 0){
        bias_row <- paste(bias_row, sprintf("%-10.4f", vals$Bias))
        sd1_row  <- paste(sd1_row,  sprintf("%-10.4f", vals$SD1))
        sd2_row  <- paste(sd2_row,  sprintf("%-10.4f", vals$SD2))
        cp_row   <- paste(cp_row,   sprintf("%-10.4f", vals$CP))
      } else {
        # In case a combination is missing, fill with NA
        bias_row <- paste(bias_row, sprintf("%-10s", "NA"))
        sd1_row  <- paste(sd1_row,  sprintf("%-10s", "NA"))
        sd2_row  <- paste(sd2_row,  sprintf("%-10s", "NA"))
        cp_row   <- paste(cp_row,   sprintf("%-10s", "NA"))
      }
    }
    
    # Print the formatted rows
    cat(bias_row, "\n")
    cat(sd1_row, "\n")
    cat(sd2_row, "\n")
    cat(cp_row, "\n")
  }
}
cat("------------------------------------------------------------------\n")



save(summary_table1, table2_summary, lambda_record, cv_result, signal_idx, beta_true_full,
     file = "C:/Users/Lenovo/Desktop/All_Simulation_Results.RData")




print(lambda_record)




my_color_palette <- c(
  "ar1" = "#fc8d62",          
  "exchangeable" = "#66c2a5",  
  "independence" = "#8da0cb"   
)

#MSE
p1 <- ggplot(summary_table1, aes(x = Method, y = MSE, fill = WorkingCorr)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~TrueRho, labeller = label_bquote(cols = "True " * rho * " = " * .(TrueRho))) +
  labs(
    title = "Mean Squared Error Comparison", 
    x = "Method", 
    y = "MSE", 
    fill = "Working Correlation"
  ) +
  theme_bw(base_size = 18) + 
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 18), 
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 16),
    
  ) +
  scale_fill_manual(values = my_color_palette)
ggsave(
  filename = "C:/Users/Lenovo/Desktop/plot_1_MSE.pdf", 
  plot = p1, 
  width = 10, 
  height = 8, 
  units = "in"
)







#TP VS FP
p3 <- ggplot(summary_table1, aes(x = FP, y = TP, color = Method, shape = WorkingCorr)) +
  geom_point(size = 8, alpha = 0.9) +   
  facet_wrap(~TrueRho, ncol = 1, labeller = label_bquote("True " * rho * " = " * .(TrueRho))) +
  geom_hline(yintercept = 5, linetype = "dotted", color = "grey60", size = 0.7) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey60", size = 0.7) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  scale_y_continuous(limits = c(4.97, 5.03)) +  
  labs(
    title = "True Positives vs False Positives",
    x = "False Positives", y = "True Positives",
    color = "", shape = ""
  ) +
  theme_bw(base_size = 18) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 22),
    strip.text = element_text(face = "bold", size = 20),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18),
    legend.position = "top",
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 10),
    panel.spacing = unit(0.3, "lines")
  ) +
  scale_color_brewer(type = "qual", palette = "Dark2")
# TP / FP
ggsave(
  filename = "C:/Users/Lenovo/Desktop/plot_3_TP_FP.pdf", 
  plot = p3, 
  width = 10, 
  height = 8, 
  units = "in"
)



# Bias
beta_labels <- c(expression(beta[1]), expression(beta[2]),
                 expression(beta[3]), expression(beta[4]),
                 expression(beta[5]))

p4 <- ggplot(table2_summary, aes(x = factor(beta_idx), y = Bias, fill = WorkingCorr)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~TrueRho, labeller = label_bquote("True " * rho * " = " * .(TrueRho)))+  
labs(
    title = "PGEE: Bias of Coefficient Estimates",
    x = "Coefficient Index",
    y = "Absolute Bias",
    fill = "Working Correlation"
  ) +
  theme_bw(base_size = 18) + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 22),  
    axis.title = element_text(size = 20),       
    axis.text = element_text(size = 18),       
    legend.title = element_text(size = 18),      
    legend.text = element_text(size = 16),       
    strip.text = element_text(face = "bold", size = 20),  
    legend.position = "right"
  ) +
  scale_x_discrete(labels = beta_labels) +  
  scale_fill_manual(values = my_color_palette)

# Bias
ggsave(
  filename = "C:/Users/Lenovo/Desktop/plot_4_Bias.pdf", 
  plot = p4, 
  width = 12, 
  height = 8, 
  units = "in"
)







library(cowplot)

row1 <- plot_grid(p1, p2, ncol = 2, align = 'hv')

row2 <- plot_grid(p3, p4, ncol = 2, align = 'hv')

row3 <- plot_grid(p5, p6, ncol = 2, align = 'hv')


final_aligned_plot <- plot_grid(row1, row2, row3, nrow = 3)



ggsave(
  filename = "C:/Users/Lenovo/Desktop/my_perfectly_aligned_plots.pdf.pdf",
  plot = final_aligned_plot,
  width = 12, 
  height = 15
)



summary_table1_copy <- summary_table1
summary_table1_copy$Condition <- paste(summary_table1_copy$Method, summary_table1_copy$WorkingCorr, sep = "_")


beta_labels <- c(expression(beta[1]), expression(beta[2]),
                 expression(beta[3]), expression(beta[4]),
                 expression(beta[5]))
p5 <- ggplot(table2_summary, aes(x = factor(beta_idx), y = CP, fill = WorkingCorr)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +  
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red", size = 1) +
  annotate("text", x = 1.2, y = 0.97, label = "95%", 
           color = "red", size = 4, hjust = 0, fontface = "bold") +  
  facet_wrap(~TrueRho, labeller = label_bquote(cols = "True " * rho * " = " * .(TrueRho)))+
  labs(
    title = "PGEE: 95% Confidence Interval Coverage Probability",
    x = "Coefficient", y = "Coverage Probability",
    fill = "Working Correlation"
  ) +
  scale_x_discrete(labels = parse(text = beta_labels)) +
  scale_fill_manual(values = my_color_palette)+
theme_bw(base_size = 18) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 22),
    strip.text = element_text(face = "bold", size = 20),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.position = "right"
  ) +
  ylim(0, 1)
ggsave(
  filename = "C:/Users/Lenovo/Desktop/plot_5_Coverage_Probability.pdf", 
  plot = p5, 
  width = 10, 
  height = 8, 
  units = "in"
)


                   
beta_labels <- c(expression(beta[1]), expression(beta[2]),
                 expression(beta[3]), expression(beta[4]),
                 expression(beta[5]))


beta_labels <- c(expression(beta[1]), expression(beta[2]),
                 expression(beta[3]), expression(beta[4]),
                 expression(beta[5]))


table2_sd_long <- table2_summary %>%
  select(TrueRho, WorkingCorr, beta_idx, SD1, SD2) %>%
  pivot_longer(cols = c("SD1", "SD2"),
               names_to = "SD_Type", values_to = "Value") %>%
  mutate(
    PanelLabel = case_when(
      TrueRho == 0.4 & WorkingCorr == "ar1" ~ "rho==0.4~','~AR(1)",
      TrueRho == 0.4 & WorkingCorr == "exchangeable" ~ "rho==0.4~','~Exchangeable",
      TrueRho == 0.4 & WorkingCorr == "independence" ~ "rho==0.4~','~Independence",
      TrueRho == 0.8 & WorkingCorr == "ar1" ~ "rho==0.8~','~AR(1)",
      TrueRho == 0.8 & WorkingCorr == "exchangeable" ~ "rho==0.8~','~Exchangeable",
      TrueRho == 0.8 & WorkingCorr == "independence" ~ "rho==0.8~','~Independence",
      TRUE ~ paste0("rho==", TrueRho, "~','~", WorkingCorr)
    )
  )


p_sd_comparison <- ggplot(table2_sd_long, aes(x = factor(beta_idx), y = Value, fill = SD_Type)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  facet_wrap(~ PanelLabel, labeller = label_parsed, nrow = 2) +
  scale_x_discrete(labels = beta_labels) +
  scale_fill_manual(
    values = c("SD1" = "#4E79A7", "SD2" = "#F28E2B"),
    labels = c("SD1 (Model SE)", "SD2 (Empirical SD)")
  ) +
  labs(
    title = "Comparison of Estimated SE vs. Empirical SD in PGEE",
    subtitle = "SD1 = Mean of Robust SE, SD2 = SD of Estimates. Closer is better.",
    x = "Coefficient Index (β)",
    y = "Standard Deviation",
    fill = "Type of SD"
  ) +
  theme_bw(base_size = 18) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 22),
    plot.subtitle = element_text(hjust = 0.5, size = 18),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    strip.text = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16)
  )


ggsave(
  filename = "C:/Users/Lenovo/Desktop/plot_6_SD_Comparison.pdf", 
  plot = p_sd_comparison, 
  width = 12, 
  height = 8, 
  units = "in"
)








                   
library(ggplot2)
df_lambda <- data.frame(lambda = lambda_record)

median_value <- median(df_lambda$lambda, na.rm = TRUE)


p_lambda <- ggplot(df_lambda, aes(x = "", y = lambda)) +
  geom_boxplot(fill = "lightblue", color = "darkblue", width = 0.3) +
  geom_hline(yintercept = median_value, color = "red", linetype = "dashed", size = 1.2) +
  annotate("text", x = 1.2, y = median_value+0.01,
           label = paste("Median =", round(median_value, 4)),
           color = "red", size = 5, hjust = 0) +
  labs(
    title = "Boxplot of Optimal Lambda",
    y = expression(paste("Optimal Lambda (", lambda, ")")),
    x = NULL
  ) +
  theme_bw(base_size = 18) +
  theme(
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

ggsave(
  filename = "C:/Users/Lenovo/Desktop/plot_7_boxplot_lambda.pdf",
  plot = p_lambda,
  width = 8,
  height = 6,
  units = "in"
)







tuning_data <- data.frame(
  lambda = cv_result$lam.vect,
  cv_error = cv_result$cv.vect
)

n_total <- n * m  
tuning_data$cv_error <- tuning_data$cv_error / n_total
min_error <- min(tuning_data$cv_error)

lam_opt <- cv_result$lam.opt


p_sensitivity <- ggplot(tuning_data, aes(x = lambda, y = cv_error)) +
  geom_line(color = "black", size = 1.5) +
  geom_point(color = "black", size = 3) +
  geom_vline(xintercept = lam_opt, linetype = "dashed", color = "red", size = 1.2) +
  geom_hline(yintercept = min_error, linetype = "dotted", color = "grey40", size = 1) +
  annotate(
    "text",
    x = lam_opt,
    y = min_error + 0.05 * (max(tuning_data$cv_error) - min_error),
    label = bquote("Optimal " * lambda * " = " * .(round(lam_opt, 4))),
    color = "red", size = 6, hjust = -0.1, fontface = "bold"
  ) +
  labs(
    title = bquote("PGEE Performance Sensitivity to " * lambda),
    subtitle = "Cross-Validation Error vs. Lambda Tuning Parameter",
    x = expression("Lambda (" * lambda * ")"),
    y = "Average Cross-Validation Error per Observation"  
  ) +
  theme_bw(base_size = 18) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    plot.subtitle = element_text(hjust = 0.5, size = 16),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16)
  )

ggsave(
  filename = "C:/Users/Lenovo/Desktop/plot_8_PGEE_Performance_Sensitivity_to_Lambda.pdf",
  plot = p_sensitivity,
  width = 8, height = 6, units = "in"
)


summary(cv_result$cv.vect)











