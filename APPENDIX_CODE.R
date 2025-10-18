### === LOAD NECESSARY LIBRARIES =============================================== ###
load_libraries <- function() {
  pkgs <- c("arothron", "dplyr", "tidyr", "ggplot2", "geomorph", "Morpho", "RRPP", "Rvcg",
            "rgl", "scatterplot3d", "car", "MASS", "mclust", "boot", "vegan",
            "progress", "ggpubr", "jsonlite", "abind")
  lapply(pkgs, require, character.only = TRUE)
}
load_libraries()
### === FUNCTION TO CONVERT AND READ LANDMARK DATA FROM 3DSLICER =============== ###
read.markups.json = function(file=NULL){
  dat = fromJSON(file, flatten=TRUE)
  n = length(dat$markups$controlPoints[[1]]$position)
  labels = dat$markups$controlPoints[[1]]$label
  temp = array(dim = c(n, 3), dimnames = list(labels, c("X", "Y", "Z")))
  for (i in 1:n) temp[i,] = dat$markups$controlPoints[[1]]$position[[i]]
  return(temp)
}
### === LOAD METADATA AND LANDMARKS ============================================ ###
metadata <- read.csv("YOUR FILE", stringsAsFactors = FALSE)
metadata$deidentified_record_number <- as.character(metadata$deidentified_record_number)
input_dir <- "YOUR .JSON LANDMARK DIRECTORY"
output_dir <- file.path(dirname(input_dir), "txt_converted")
if (!dir.exists(output_dir)) dir.create(output_dir)
files <- list.files(input_dir, pattern = ".mrk.json$", full.names = TRUE)
ids_raw <- tools::file_path_sans_ext(basename(files))
ids <- gsub("_.*", "", ids_raw) 
landmarks <- list()
group <- c()
failed_ids <- c()
n_landmarks_expected <- NULL
### === Process each landmark file ============================================= ###
for (i in seq_along(files)) {
  id <- ids[i]
  
  if (!(id %in% meta$deidentified_record_number)) {
    warning(paste("ID not found in metadata:", id))
    failed_ids <- c(failed_ids, id)
    next
  }
  
  sex <- meta$sex_code[meta$deidentified_record_number == id]
  lm <- tryCatch(read.markups.json(files[i]), error = function(e) {
    warning(paste("Failed to read file:", files[i]))
    return(NULL)
  })
  
  if (is.null(lm)) {
    failed_ids <- c(failed_ids, id)
    next
  }
  
  if (is.null(n_landmarks_expected)) {
    n_landmarks_expected <- nrow(lm)
  }
  if (nrow(lm) != n_landmarks_expected) {
    warning(paste("Landmark count mismatch for", id, "- Found:", nrow(lm),
                  "Expected:", n_landmarks_expected))
    failed_ids <- c(failed_ids, id)
    next
  }
  
  # Save .txt version
  write.table(lm,
              file = file.path(output_dir, paste0(id, ".txt")),
              sep = ",", col.names = FALSE, row.names = FALSE)
  
  # Store for analysis
  landmarks[[id]] <- lm
  group <- c(group, sex)
}
# Convert to array and group to factor
landmark_array <- abind(lapply(landmarks, as.matrix), along = 3)
group <- factor(group, levels = c("Male", "Female"))
### === GENERALISED PROCRUSTES ANALYSIS (GPA) ================================== ###
gpa_result <- gpagen(landmark_array, print.progress = TRUE)
### === PRINCIPLE COMPONENTS ANALYSIS (PCA) ==================================== ###
pca_result <- gm.prcomp(gpa_result$coords)
### === PC variance explained ================================================== ###
variance_explained <- (pca_result$sdev^2) / sum(pca_result$sdev^2) * 100
variance_explained #Lists variance explained 
cumulative_variance <- cumsum(variance_explained)
plot(variance_explained, type = "b", pch = 16,
     xlab = "Principal Component",
     ylab = "Variance Explained (%)",
     main = "Scree Plot")
plot(cumulative_variance, type = "b", pch = 16,
     xlab = "Principal Component",
     ylab = "Cumulative Variance Explained (%)",
     main = "Cumulative Variance Plot")
### === Define sex-group colours for plots === ###
female_color <- "#FF8A9A"  
male_color <- "#4FC3F7" 
### === Plot PCA with convex hulls === ###
plot(pca_result,
     pch = 21,
     bg = c(male_color, female_color)[group],
     col = "black",
     cex = 1.4,
     lwd = 1.2,
     asp = 1,
     main = "PCA of Basicranial Shape: All Landmarks")
pca_scores <- pca_result$x
male_scores <- pca_scores[group == "Male", 1:2]
female_scores <- pca_scores[group == "Female", 1:2]
polygon(male_scores[chull(male_scores), ],                                     # Repeat polygon for "female_scores"
        col = adjustcolor(male_color, alpha.f = 0.2),
        border = male_color, lwd = 1)
legend("topright",
       legend = levels(group),
       pch = 21,
       pt.bg = c(male_color, female_color), 
       pt.cex = 1.3)
### === calculate centroid size (CS) =========================================== ###
centroid_size <- gpa_result$Csize
male_coords <- gpa_result$coords[, , group == "Male"]                          # Repeat for "female_coords"
male_cs <- centroid_size[group == "Male"]                                      # Repeat for "female_cs"
### === calculate means by group === ###
mean_male_cs <- round(mean(centroid_size[group == "Male"], na.rm = TRUE), 2)   # Repeat for "mean_female_cs"
### === t-test for CS dimorphism === ###
t_test_centroid <- t.test(centroid_size ~ group)
p_val_cs <- ifelse(t_test_centroid$p.value < 0.001, "<0.001", round(t_test_centroid$p.value, 3))
### === CS violin plot with subtitle === ###  
subtitle_cs <- paste0("Mean Male: ", mean_male_cs, " | Mean Female: ", mean_female_cs, " | p = ", p_val_cs)
ggplot(data.frame(centroid_size, group), aes(x = group, y = centroid_size, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
  geom_jitter(width = 0.1, alpha = 0.4) +
  scale_fill_manual(values = c("Male" = "#4FC3F7", "Female" = "#FF8A9A")) +
  labs(
    title = "Centroid Size by Sex",
    subtitle = subtitle_cs,
    x = "Sex",
    y = "Centroid Size"
  ) +
  theme_minimal() +
  theme(legend.position = "none")
### === Split CS by sex === ###
cs_split <- split(centroid_size, group)
log_CS <- log(gpa_result$Csize)
cs_data <- data.frame(
  group = factor(group, levels = c("Male", "Female")),
  centroid_size = centroid_size
)
### === Simple allometry model === ###
shape_allo_sex <- procD.lm(gpa_result$coords ~ log(centroid_size) + group, iter = 999)
summary(shape_allo_sex)
### === calculate size-shape regression scores === ###
male_model <- procD.lm(male_coords ~ male_cs, iter = 999)                                       # Repeat for "female_model" using "female_coords"
log_cs_male <- log(centroid_size[group == "Male"])                                              # Repeat for "log_cs_female"
male_output <- plotAllometry(male_model, size = log_cs_male, method = "RegScore", plot = FALSE) # Repeat for "female_output"
male_scores <- male_output$RegScore
female_scores <- female_output$RegScore
df_male <- data.frame(
  Log_CS = log_cs_male,
  RegScore = male_scores,
  Sex = "Male"
)
                                                                               # Repeat for "df_female"
### === combine into one data frame === ###
regression_df <- rbind(df_male, df_female)
df_allo <- rbind(df_male, df_female)
### === Create scatterplot using ggplot2 (you will need to extract R-squared and p-values for the plot) === ###
### === RMSE FOR INTRAOBSERVER-ERROR =========================================== ###
mean_shape <- mshape(gpa_result$coords)
compute_rmse <- function(specimen_coords, mean_coords) {
  sqrt(mean((specimen_coords - mean_coords)^2))
}
rmse_values <- apply(gpa_result$coords, 3, compute_rmse, mean_coords = mean_shape)
rmse_df <- data.frame(ID = confirmed_ids, RMSE = rmse_values, Sex = confirmed_sex)
### === RMSE FOR INTER-OBSERVER ERROR ========================================== ###
### === Load and convert landmark data for Inter-observer RMSE using above code === ###
### === Perform GPA on new landmark files === ###
mean_shape <- mshape(gpa_result$coords)
compute_rmse <- function(specimen_coords, mean_coords) {
  sqrt(mean((specimen_coords - mean_coords)^2))
}
rmse_values <- apply(gpa_result$coords, 3, compute_rmse, mean_coords = mean_shape)
rmse_df <- data.frame(ID = confirmed_ids, RMSE = rmse_values, Sex = confirmed_sex)
### === Identify potential outliers (e.g., above 95th percentile) === ###
threshold <- quantile(rmse_values, 0.95)
outliers <- rmse_df[rmse_df$RMSE > threshold, ]
### === ALLOMETRY ANALYSIS ===================================================== ### 
### === Extract PCs === ###
pc1_scores <- pca_scores[,1]                                                   # Repeat to extract "pc2_scores" using "[, 2]" and "pc3_scores" using "[, 3]"
male_pc1 <- pc1_scores[group == "Male"]                                        # Repeat for "Female"
male_pc2 <- pc2_scores[group == "Male"]                                        # Repeat for "Female"
### === Sex regression for PC1 and PC2 === ###
male_pc1_model <- lm(male_pc1 ~ log_cs_male)                                   # Repeat for "Female"
male_pc2_model <- lm(male_pc2 ~ log_cs_male)                                   # Repeat for "Male"
### === R-squared and p-values === ###
r2_male_pc1 <- round(summary(male_pc1_model)$r.squared, 3)                     # Repeat for "Female"
p_male_pc1 <- round(summary(male_pc1_model)$coefficients[2, 4], 4)             # Repeat for "Female"
r2_male_pc2 <- round(summary(male_pc2_model)$r.squared, 3)                     # Repeat for "Female"
p_male_pc2 <- round(summary(male_pc2_model)$coefficients[2, 4], 4)             # Repeat for "Female"
### === Prepare data frames === ### 
df_male_pc1 <- data.frame(
  Log_CS = log_cs_male,
  PC_Score = male_pc1,
  Sex = "Male"
)                                                                              # Repeat for "df_female_pc1"
                                                                               # Repeat dataframes for "male_PC2" and "female_PC2"
### === Combine data for plotting === ###
df_pc1 <- rbind(df_male_pc1, df_female_pc1)                                    # Repeat for "df_pc2" using "df_male_pc2, df_female_pc2"
### === Plot PC1 vs log(CS) by sex === ###
ggplot(df_pc1, aes(x = Log_CS, y = PC_Score, color = Sex)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_color_manual(values = c("Male" = "#4FC3F7", "Female" = "#FF8A9A")) +
  labs(
    title = "Allometry: PC1 vs log(CS)",
    subtitle = paste0(
      "Male: R² = ", r2_male_pc1, ", p = ", p_male_pc1, " | ",
      "Female: R² = ", r2_female_pc1, ", p = ", p_female_pc1
    ),
    x = "log(CS)",
    y = "PC1"
  ) +
  theme_minimal()
                                                                               # Repeat allometry plot for PC2 vs log(CS) by sex using "df_PC2"
### === RESIDUALS ANALYSIS: PC1/PC2 vs log(CS) ================================= ###
### Assuming merged_df_clean contains PC1, PC2, log_CS, and sex_code (like merged_df for ECV)
### Fit combined interaction models (PC1 and PC2 ~ log_CS * sex_code)
lm_pc1_cs_interaction <- lm(PC1 ~ log_CS * sex_code, data = merged_df_clean)   # Repeat for "lm_pc2_cs_interaction"
### Extract residuals from interaction models and assign to dataframe
merged_df_clean$resid_pc1_cs <- resid(lm_pc1_cs_interaction)                   # Repeat for "resid_pc2_cs"
### === ALLOMETRY ANALYSIS: PCs vs ECV ========================================= ###
### === PC1 & PC2 regressions === ###
lm_male_pc1 <- lm(log_ECV ~ PC1, data = male_data)                             # Repeat for "female_data"
lm_male_pc2 <- lm(log_ECV ~ PC2, data = male_data)                             # Repeat for "male_data"
### === Extract data for plot === ###
male_pc1_r2 <- summary(lm_male_pc1)$r.squared                                  # Repeat for "female_data"
male_pc1_p  <- summary(lm_male_pc1)$coefficients[2, 4]                         # Repeat for "female_data"
male_pc2_r2 <- summary(lm_male_pc2)$r.squared                                  # Repeat for "female_data"
male_pc2_p  <- summary(lm_male_pc2)$coefficients[2, 4]                         # Repeat for "female_data"
### === Plotting, e.g. PC1 vs log(ECV) by sex === ###
ggplot(merged_df_clean, aes(x = log_ECV, y = PC1, color = sex_code)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", aes(group = sex_code), se = TRUE, linewidth = 1.1) +
  scale_color_manual(values = c("Male" = "#4FC3F7", "Female" = "#FF8A9A")) +
  labs(
    title = "Allometry: PC1 vs log(ECV)",
    subtitle = paste0(
      "Male: R² = ", round(male_pc1_r2, 3), 
      ", p = ", signif(male_pc1_p, 3), " | ",
      "Female: R² = ", round(female_pc1_r2, 3), 
      ", p = ", signif(female_pc1_p, 3)
    ),
    x = "log(ECV)",
    y = "PC1",
    color = "Sex"
  ) +
  theme_minimal()
### === Repeat plot for PC2 vs log(ECV) by sex === ###
### === Linear regression results (combined sex) === ###
lm_log_pc1 <- lm(PC1 ~ log_ECV, data = merged_df_clean) # Repeat for PC2
### === Combine with a "Sex" column for plotting === ###
resid_df <- bind_rows(
  male_data %>% mutate(Sex = "Male"),
  female_data %>% mutate(Sex = "Female")
)
### === AGE CORRELATION ======================================================== ###
lm_age <- lm(PC1 ~ age_years, data = merged_df)
### === If necessary, remove rows where 'sex_code' is NA === ###
merged_df_clean <- merged_df %>% filter(!is.na(sex_code))
### === RESIDUALS ANALYSIS: PC1/PC2 vs log(ECV) ================================ ###
lm_pc1_interaction <- lm(PC1 ~ log_ECV * sex_code, data = merged_df_clean)
merged_df_clean$resid_pc1 <- resid(lm_pc1_interaction)                         # Repeat for PC2 vs log(ECV)
### === Perform separate sex-specific regressions e.g. === ###
lm_male_pc2 <- lm(PC2 ~ log_ECV, data = male_data)                             # Repeat for "lm_female_pc2"
male_data$resid_pc2 <- resid(lm_male_pc2)                                      # Repeat for "female_data$resid_pc2 <- resid(lm_female_pc2)"
resid_df <- bind_rows(
  male_data %>% mutate(Sex = "Male"),
  female_data %>% mutate(Sex = "Female")
)
summary(lm(resid_pc2 ~ sex_code, data = merged_df_clean))
### === THIN PLATE SPLINE (TPS) MODELS ========================================= ###
par(mfrow = c(2, 2))
plotRefToTarget(gpa_result$coords[, , which.max(pca_scores[,1])],
                gpa_result$coords[, , which.min(pca_scores[,1])],
                method = "TPS", mag = 1.5, main = "PC1 Positive vs Negative")
plotRefToTarget(gpa_result$coords[, , which.max(pca_scores[,2])],
                gpa_result$coords[, , which.min(pca_scores[,2])],
                method = "TPS", mag = 1.5, main = "PC2 Positive vs Negative")
male_coords <- gpa_result$coords[, , group == "Male"]                          # Repeat for "female_coords"
mean_male <- apply(male_coords, c(1, 2), mean)                                 # Repeat for "mean_female"
plotRefToTarget(mean_male, mean_female, method = "TPS", mag = 1.5, main = "TPS Warp: Male → Female")
plotRefToTarget(mean_male, mean_female, method = "vector", mag = 1.5, main = "Deformation Vectors: Male → Female")
### === 3D PCA: PC1 vs PC2 vs PC3 ============================================== ###
colors <- ifelse(group == "Male", male_color,
                 ifelse(group == "Female", female_color, "gray"))
add_hull_3d <- function(scores, color, alpha = 0.2) {
  hull <- convhulln(scores, options = "FA")
  hull_faces <- hull$hull
  triangles3d(scores[hull_faces, ],
              col = adjustcolor(color, alpha.f = alpha),
              alpha = alpha,
              front = "fill", back = "fill")
}
# Prepare scores matrices for groups
male_scores_3d <- cbind(pc1[group == "Male"],
                        pc2[group == "Male"],
                        pc3[group == "Male"])
female_scores_3d <- cbind(pc1[group == "Female"],
                          pc2[group == "Female"],
                          pc3[group == "Female"])
# Plot PCA points
plot3d(pc1, pc2, pc3,
       col = colors,
       type = "s",
       radius = 0.002,
       xlab = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
       ylab = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
       zlab = paste0("PC3 (", round(variance_explained[3], 1), "%)"))
# Add convex hulls
add_hull_3d(male_scores_3d, male_color)
add_hull_3d(female_scores_3d, female_color)
# Add legend
legend3d("topright",
         legend = levels(group),
         pch = 16,
         col = c(male_color, female_color),
         cex = 1.2,
         inset = c(0.02))
### === size (CS) vs Shape (PD) allometry by sex groups separately ============= ###
proc_result_male <- procD.lm(male_coords ~ male_cs, iter = 999)                # Repeat for "proc_result_female" using "female_coords"
### === Procrustes ANOVA ======================================================= ###
proc_result <- procD.lm(gpa_result$coords ~ group, iter = 999)
### === PROCRUSTES ANOVA WITH AGE - DOES SD CHANGE WITH AGE? =================== ### 
n_specimens <- dim(gpa_result$coords)[3]
### === Trim merged_df to match the coordinates array === ###
merged_trim <- merged_df[1:n_specimens, ]
## ===# Create age group bins (young, middle, old) === ###
merged_trim$age_years <- as.numeric(as.character(merged_trim$age_years))
merged_trim$age_group <- cut(
  merged_trim$age_years,
  breaks = c(18, 40, 65, 100),
  labels = c("Young", "Middle", "Old"),
  right = TRUE, include.lowest = TRUE
)
valid_idx <- !is.na(merged_trim$age_group)
coords_sub <- gpa_result$coords[ , , valid_idx]
age_data <- merged_trim[valid_idx, ]
proc_result_age <- procD.lm(coords_sub ~ age_group, data = age_data, iter = 999)
### === CALCULATE MAHALANOBIS DISTANCE BETWEEN SEXES USING FIRST 17 PCs ======== ###
scores_17pc <- pca_scores[, 1:17]
scores_male <- scores_17pc[group == "Male", , drop = FALSE]                    # Repeat "scores_female"
### === Calculate group centroids === ###
centroid_male <- colMeans(scores_male)                                         # Repeat for "centroid_female" using "scores_female"
pooled_cov <- cov(rbind(scores_male, scores_female))
### === Mahalanobis distance between the two group centroids === ###
mahal_dist <- mahalanobis(centroid_male, center = centroid_female, cov = pooled_cov)
cat("Mahalanobis distance (Male vs Female, PC1–PC17):", round(mahal_dist, 3), "\n")
### === CALCULATE SEXUAL DIMORPHISM IN PROCRUSTES DISTANCES (PD) =============== ### 
male_coords   <- gpa_result$coords[, , group == "Male"] # Repeat for "female_coords"
mean_male     <- apply(male_coords, c(1, 2), mean) # Repeat for "mean_female" 
procrustes_distance <- sqrt(sum((mean_male - mean_female)^2))
cat("Procrustes Distance (Mean Male vs Female):", round(procrustes_distance, 4), "\n")
### === Calculate Euclidean distance for each landmark === ###
landmark_differences <- sqrt(rowSums((mean_male - mean_female)^2))
### === Create a dataframe and sort by largest difference === ###
landmark_df <- data.frame(
  Landmark = rownames(mean_male),
  Distance = landmark_differences
)
landmark_df <- landmark_df %>%
  arrange(desc(Distance))
### === TOP 20 MOST SEXUALLY DIMORPHIC LANDMARKS =============================== ###
top_n <- 20
top_landmarks <- landmark_df %>% slice(1:top_n)
ggplot(top_landmarks, aes(x = reorder(as.factor(Landmark), -Distance), y = Distance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = paste("Top", top_n, "Landmarks by Distance"),
       x = "Landmark",
       y = "Euclidean Distance") +
  theme_minimal()
### === Bootstrapped LDA (Use for first 5 or 17 PCs) =========================== ###
set.seed(123)
n_iter <- 1000
n_per_group <- 60
accuracy_vec <- numeric(n_iter)
pb <- progress_bar$new(
  format = "Bootstrapping [:bar] :percent ETA: :eta",
  total = n_iter, clear = FALSE, width = 60
)
for (i in 1:n_iter) {
  male_indices <- sample(which(group == "Male"), n_per_group, replace = TRUE)
  female_indices <- sample(which(group == "Female"), n_per_group, replace = TRUE)
  sampled_indices <- c(male_indices, female_indices)
  
  coords_sampled <- landmark_array[, , sampled_indices]
  group_sampled <- group[sampled_indices]
  
  gpa_result <- gpagen(coords_sampled, print.progress = FALSE)
  pca_result <- gm.prcomp(gpa_result$coords)
  pca_scores <- pca_result$x[, 1:5]  # Using 5 or 17 PCs
  
  lda_data <- data.frame(Sex = as.factor(group_sampled), pca_scores)
  lda_model <- lda(Sex ~ ., data = lda_data)
  lda_pred <- predict(lda_model)
  
  accuracy_vec[i] <- mean(lda_pred$class == lda_data$Sex)
  
  pb$tick()  
}
# === Summarize results === 
mean_acc <- mean(accuracy_vec)
ci <- quantile(accuracy_vec, probs = c(0.025, 0.975))
cat("\nBootstrapped LDA (5 PCs):\n")
cat("Mean Accuracy:", round(mean_acc * 100, 2), "%\n")
cat("95% CI:", round(ci[1] * 100, 2), "% -", round(ci[2] * 100, 2), "%\n")
# === visualise LDA (histogram) ===
accuracy_percent <- accuracy_vec * 100
ggplot(data.frame(Accuracy = accuracy_percent), aes(x = Accuracy)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black", alpha = 0.7) +
  geom_vline(aes(xintercept = mean(accuracy_percent)), color = "red", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = ci[1] * 100), color = "green", linetype = "dotted", size = 1) +
  geom_vline(aes(xintercept = ci[2] * 100), color = "green", linetype = "dotted", size = 1) +
  labs(title = "Title",
       x = "Classification Accuracy (%)",
       y = "Frequency") +
  theme_minimal()

### === LEFT-RIGHT BASICRANIUM LATERALITY ANALYSIS ============================= ### 
### === LOAD LEFT-RIGHT LANDMARKING PROTOCOL === ###
laterality_df <- read.csv("YOUR_FILE.csv", stringsAsFactors = FALSE)
for (id in names(landmarks)) {
  # reorder rows by label
  landmarks[[id]] <- landmarks[[id]][expected_labels, , drop=FALSE]
}
### === Convert to 3D array === ###
landmark_array <- abind(lapply(landmarks, as.matrix), along = 3)
dimnames(landmark_array)[[1]] <- expected_labels
dimnames(landmark_array)[[3]] <- names(landmarks)
### === Assign side vector for all landmarks (same for all specimens) === ###
lm_sides <- laterality_df$Side
names(lm_sides) <- expected_labels
### === Subset landmarks by left-right side === ###
left_indices <- which(lm_sides == "Left")                                      # Repeat for "right_indices" using lm_sides == "Right"
left_landmarks  <- landmark_array[left_indices, , , drop=FALSE]                # Repeat for "right_landmarks" using "right_indices"
### === GPA on left-right === ###
gpa_left  <- gpagen(left_landmarks, print.progress = TRUE)                     # Repeat for "gpa_right"
### === Run left-right PCA === ###
pca_left <- gm.prcomp(gpa_left$coords)                                         # Repeat for "pca_right"
### === Extract PC scores (e.g. PC1-PC10) === ###
max_pcs <- min(ncol(pca_left$x), ncol(pca_right$x), 10)
left_scores <- pca_left$x[, 1:max_pcs]                                         # Repeat for "right_scores" using "pca_right"
### === Calculate total variance (sum of variances per PC) === ###
left_variance <- sum(apply(left_scores, 2, var))                               # Repeat for "right_variance" using "right_scores"
### === Statistical testing === ###
combined_scores <- rbind(left_scores, right_scores)
group_labels <- factor(c(rep("Left", nrow(left_scores)), rep("Right", nrow(right_scores))))
n_perm <- 999
perm_diffs <- numeric(n_perm)
obs_diff <- abs(left_variance - right_variance)
set.seed(123) 
for (i in 1:n_perm) {
  shuffled <- sample(group_labels)
  left_grp <- combined_scores[shuffled == "Left", , drop = FALSE]
  right_grp <- combined_scores[shuffled == "Right", , drop = FALSE]
  
  v_left <- sum(apply(left_grp, 2, var))
  v_right <- sum(apply(right_grp, 2, var))
  
  perm_diffs[i] <- abs(v_left - v_right)
}
### === Laterality by sex === ###
specimen_sex <- meta$sex_code[match(names(landmarks), meta$deidentified_record_number)]
get_scores <- function(side = c("Left","Right"), sex = c("Male","Female")){
  s   <- match.arg(side)
  sx  <- match.arg(sex)
  pcs <- if(s == "Left") left_scores else right_scores
  pcs[specimen_sex == sx, , drop = FALSE]
}
left_male   <- get_scores("Left",  "Male")                                     # Repeat for "left_female" 
right_male  <- get_scores("Right", "Male")                                     # Repeat for "right_female"
group_var <- function(mat) sum(apply(mat, 2, var))
var_L_M <- group_var(left_male)                                                # Repeat for "var_L_F" using "left_female"
var_R_M <- group_var(right_male)                                               # Repeat for "var_R_F" using "right_female"
print(rbind(
  Left  = c(Male = var_L_M, Female = var_L_F),
  Right = c(Male = var_R_M, Female = var_R_F)
))
perm_var_test <- function(mat, sexes, n_perm = 999){
  obs <- abs(group_var(mat[sexes=="Male", ]) -
               group_var(mat[sexes=="Female", ]))
  perm <- replicate(n_perm, {
    shuf <- sample(sexes)
    abs(group_var(mat[shuf=="Male", ]) -
          group_var(mat[shuf=="Female", ]))
  })
  mean(perm >= obs)
}
p_left  <- perm_var_test(left_scores,  specimen_sex)                           # Repeat for "p_right" using "right_scores"
cat("Left‑side Male–Female variance p‑value :", round(p_left,  4), "\n")
cat("Right‑side Male–Female variance p‑value:", round(p_right, 4), "\n")
# Test shape difference by sex for LEFT side
proc_sex_left <- procD.lm(gpa_left$coords ~ specimen_sex, iter = 999)          # Repeat for "proc_sex_right" using "gpa_right"
# Combine scores and labels
combined_scores <- rbind(left_scores, right_scores)
combined_sex    <- rep(specimen_sex, 2)
combined_side   <- rep(c("Left", "Right"), each = length(specimen_sex))
# Build data frame
df <- data.frame(
  sex = factor(combined_sex),
  side = factor(combined_side)
)
### === Plot variance following ggplot2 template in this script === ###
### === ANTERIOR AND POSTERIOR BASICRANIAL SHAPE ANALYSIS ====================== ###
### === LOAD ANTERIOR-POSTERIOR LANDMARKING PROTOCOL === ###
landmark_regions <- read.csv("YOUR_FILE.csv", stringsAsFactors = FALSE)
anterior_labels <- landmark_regions$Label[landmark_regions$Region == "Anterior"]  
posterior_labels <- landmark_regions$Label[landmark_regions$Region == "Posterior"]
### === Extract arrays === ###
anterior_array <- abind(lapply(landmarks, \(x) x[anterior_labels, , drop=FALSE]), along=3)
posterior_array <- abind(lapply(landmarks, \(x) x[posterior_labels, , drop=FALSE]), along=3)
### === Anterior-posterior GPA and PCA === ###
gpa_anterior <- gpagen(anterior_array)                                                 # Repeat for "gpa_posterior" using "posterior_array"
pca_anterior <- gm.prcomp(gpa_anterior$coords)                                         # Repeat for "pca_posterior" using "gpa_posterior"
### === Estimate variance === ###
max_pcs <- min(ncol(pca_anterior$x), ncol(pca_posterior$x), 10)
var_anterior <- sum(apply(pca_anterior$x[, 1:max_pcs], 2, var))                        # Repeat for "var_posterior" using "pca_posterior"
### === Perform permutation tests === ### 
obs_diff <- abs(var_anterior - var_posterior)
perm_diffs <- replicate(999, {
  idx <- sample(c(rep("A", nrow(pca_anterior$x)), rep("P", nrow(pca_posterior$x))))
  v1 <- sum(apply(rbind(pca_anterior$x, pca_posterior$x)[idx == "A", 1:max_pcs], 2, var))
  v2 <- sum(apply(rbind(pca_anterior$x, pca_posterior$x)[idx == "P", 1:max_pcs], 2, var))
  abs(v1 - v2)
})
p_value <- mean(perm_diffs >= obs_diff)
### === Compute sex-based variance === ###
get_scores_by_sex <- function(scores, sex, label) sum(apply(scores[sex == label, ], 2, var))
var_ant_male <- get_scores_by_sex(pca_anterior$x[, 1:max_pcs], specimen_sex, "Male")   # Repeat for "var_ant_female"
var_post_male <- get_scores_by_sex(pca_posterior$x[, 1:max_pcs], specimen_sex, "Male") # Repeat for "var_post_female"
### === Procrustes ANOVA for anterior-posterior regions by sex === ###
proc_sex_ant <- procD.lm(gpa_anterior$coords ~ specimen_sex, iter=999)                 # Repeat for "proc_sex_post" using "gpa_posterior"
### === Per-specimen Procrustes variance === ###
calc_var <- \(coords) apply(gm.prcomp(coords)$x[, 1:max_pcs], 1, \(x) sum(x^2))
var_ant_specimen <- calc_var(gpa_anterior$coords)                                      # Repeat for "var_post_specimen" using "gpa_posterior"
### === Shape regression with ECV === ###
ecv_vector <- ecv_df$ECV[match(dimnames(gpa_anterior$coords)[[3]], ecv_df$deidentified_record_number)]
proc_reg_ant <- procD.lm(gpa_anterior$coords ~ ecv_vector, iter=999)                   # Repeat for "proc_reg_post" using "gpa_posterior"
### === Regression scores === ###
reg_score_ant <- two.d.array(gpa_anterior$coords) %*% proc_reg_ant$lmcoef[-1]          # Repeat for "reg_score_post" using "gpa_posterior"
# === Linear regression models: PC1 and PC2 ~ ECV, by region and sex === ###
# === requires extracting PC1 and PC2 scores for anterior and posterior landmarks === ### 
# === Example regression models for PC1 ant-post are shown below by sex === ###
lm_ant_male   <- lm(PC1_Anterior ~ ECV, data = subset(df_pc1, sex_code == "Male"))     # Repeat for "lm_ant_female"
lm_post_male  <- lm(PC1_Posterior ~ ECV, data = subset(df_pc1, sex_code == "Male"))    # Repeat for "lm_post_female"
### === LOCAL LANDMARK DIMENSIONS (USING JUST FIXED SINGLE LANDMARKS) ========== ###
### === Paths to .txt landmark folders of single landmarks and read them === ###
male_txt_path <- "YOUR_FILE"                                                           # Repeat for "female_txt_path"
read_txt_landmarks <- function(file) {
  as.matrix(read.table(file, sep = ",", header = FALSE))
}
### === Load and read .txt files === ###
male_files <- list.files(male_txt_path, pattern = ".txt$", full.names = TRUE)  # Repeat for "female_files"
male_landmarks <- lapply(male_files, read_txt_landmarks)                       # Repeat for "female_landmarks"
### === Convert to 3D arrays: (landmarks x 3 x specimens) === ###
male_array <- abind(male_landmarks, along = 3)                                 # Repeat for "female_array"
all_landmarks <- abind(male_array, female_array, along = 3)
group <- factor(c(rep("Male", length(male_files)), rep("Female", length(female_files))))
### === FIXED SINGLE LANDMARK GPA AND PCA === ###
gpa_result <- gpagen(all_landmarks, print.progress = TRUE)
pca_result <- gm.prcomp(gpa_result$coords)
pca_scores <- pca_result$x
male_scores <- pca_scores[group == "Male", 1:2]                                # Repeat for "female_scores"
### === Plot PCA with convex hull === ###
plot(pca_result, 
     pch = 21, 
     bg = c(female_color, male_color)[group], 
     col = "black", 
     cex = 1.4, 
     lwd = 1.2, 
     asp = 1, 
     main = "PCA of Basicranial Shape: Single Landmarks")
male_hull <- chull(male_scores)                                                # Repeat for "female_hull"
polygon(male_scores[male_hull, ], col = adjustcolor(male_color, alpha.f = 0.2), border = male_color, lwd = 1)
polygon(female_scores[female_hull, ], col = adjustcolor(female_color, alpha.f = 0.2), border = female_color, lwd = 1)
legend("topright", legend = levels(group), pch = 21, pt.bg = c(female_color, male_color), pt.cex = 1.3)
### === single landmarks Procrustes ANOVA === ###
procD.lm(gpa_result$coords ~ group, iter = 999)
### === Canonical Variates Analysis (CVA) ====================================== ###
cva_result <- CVA(gpa_result$coords, group)
cva_scores <- cva_result$CVscores
plot(cva_scores, col = c(female_color, male_color)[cva_result$groups], pch = 21,
     main = "CVA of Basicranial Shape", xlab = "CV1", ylab = "Score", asp = 1)
legend("topright", legend = levels(cva_result$groups), pch = 21, pt.bg = c(female_color, male_color))
t_test_result <- t.test(cva_result$CVscores ~ cva_result$groups)
boxplot(cva_result$CVscores ~ cva_result$groups, col = c(female_color, male_color),
        main = "Distribution of CV1 Scores by Group", xlab = "Group", ylab = "CV1 Score")
### === Extract shapes by group === ###
female_shapes <- gpa_result$coords[,,group == "Female"]                        # Repeat for "male_shapes"
### === Compute mean shape per group === ###
mean_female <- apply(female_shapes, c(1,2), mean)                              # Repeat for "mean_male"
### === Mean shape across both groups (i.e. pooled mean) === ###
mean_shape <- apply(gpa_result$coords, c(1,2), mean)
### === Mean shape differences and deformation visualization === ###
male_landmarks <- gpa_result$coords[, , group == "Male"]                       # Repeat for "female_landmarks"
mean_male <- apply(male_landmarks, c(1, 2), mean)                              # Repeat for "mean_female"
plotRefToTarget(mean_male, mean_female, method = "TPS", mag = 1.5, main = "Warps of Male to Female Shape")
plotRefToTarget(mean_male, mean_female, method = "TPS", mag = 1.5, main = "Deformation Arrows: Male to Female Shape")
procrustes_distance <- procD.lm(gpa_result$coords ~ group)
dist_male_female <- sqrt(sum((mean_male - mean_female)^2))
cat("Procrustes Distance (Mean Male vs Female):", dist_male_female, "\n")
### === Foramen Magnum (FM) Area Calculations ================================== ###
calc_fm_area <- function(coords_array) {
  apply(coords_array, 3, function(specimen) {
    length <- sqrt(sum((specimen[3, ] - specimen[4, ])^2))      # Basion (3) to Opisthion (4)
    breadth <- sqrt(sum((specimen[23, ] - specimen[24, ])^2))   # Left FM (23) to Right FM (24)
    area <- pi * (length / 2) * (breadth / 2)
    
    cat("Specimen ID: ", dimnames(coords_array)[[3]], "\n")
    cat("FM Length (Basion to Opisthion): ", length, " mm\n")
    cat("FM Breadth (Left FM to Right FM): ", breadth, " mm\n")
    cat("Calculated FM Area: ", area, " mm²\n")
    
    if (area < 1) {
      cat("Warning: Abnormally small FM area detected! Specimen ID: ", dimnames(coords_array)[[3]], "\n")
    }
    return(area)
  })
}
fm_area <- calc_fm_area(all_landmarks)
### === Specimen IDs cleanup and FM area dataframe === ###
specimen_ids <- c(basename(tools::file_path_sans_ext(male_files)), basename(tools::file_path_sans_ext(female_files)))
specimen_ids <- gsub("_skull|\\.mrk\\.json|\\.txt", "", specimen_ids)
fm_df <- data.frame(Specimen = specimen_ids, FM_Area = fm_area)
### === Merge FM Area with ECV data and run regressions === ###
ecv_df <- read.csv("YOUR FILE")
ecv_df$deidentified_record_number <- gsub("_skull|\\.mrk\\.json|\\.txt", "", ecv_df$deidentified_record_number)
fm_df$deidentified_record_number <- gsub("_skull|\\.mrk\\.json|\\.txt", "", fm_df$Specimen)
merged_df <- merge(ecv_df, fm_df[, c("deidentified_record_number", "FM_Area")], by = "deidentified_record_number", all.x = TRUE)
lm_fm <- lm(FM_Area ~ ECV, data = merged_df)
lm_fm_log <- lm(log(FM_Area) ~ log(ECV), data = merged_df)
summary(lm_fm)
summary(lm_fm_log)
### === Statistics, compare FM between sexes === ###
merged_df %>%
  group_by(sex_code) %>%
  summarise(
    mean_FM = mean(FM_Area, na.rm = TRUE),
    sd_FM = sd(FM_Area, na.rm = TRUE),
    n = n()
  )
t_test_result <- t.test(FM_Area ~ sex_code, data = merged_df)
### === Coefficient of Variation (CV) by sex === ###
merged_df %>%
  group_by(sex_code) %>%
  summarise(
    mean_FM = mean(FM_Area, na.rm = TRUE),
    sd_FM = sd(FM_Area, na.rm = TRUE),
    CV = sd_FM / mean_FM * 100
  )
### === Plot density violin plot for FM area by sex === ###
### === Normality Tests (Shapiro-Wilk) === ###
shapiro.test(merged_df$FM_Area[merged_df$sex_code == "male"])                  # Repeat for "female"
### === Wilcoxon Test (non-parametric alternative) === ###
wilcox.test(FM_Area ~ sex_code, data = merged_df)
### === WIREFRAME SHAPE FIGURES ================================================ ###
### === Define wireframe connections (created using MorphoJ initially) ========= ###
wireframe_connections <- matrix(c(
  4, 23, 3, 23, 3, 24, 4, 24, 4, 25, 3, 25, 3, 26, 4, 26,
  3, 4, 3, 21, 21, 27, 9, 27, 1, 9, 1, 10, 10, 28, 22, 28,
  3, 22, 2, 3, 1, 2, 30, 32, 18, 30, 8, 18, 6, 8, 6, 14,
  14, 16, 12, 16, 12, 20, 29, 31, 17, 29, 7, 17, 5, 7,
  5, 13, 13, 15, 11, 15, 11, 19
), ncol = 2, byrow = TRUE)
### === Extract PCA scores for males and females === ###
male_pca_scores <- pca_result$x[group == "Male", ]                # Repeat for female_pca_scores
                                                                  # Calculate male and female PC1 extremes using "[, 1]"
pc1_pos_male_idx <- which.max(male_pca_scores[, 1])               # Repeat for "pc1_pos_female_idx" using "female_pca_scores"
pc1_neg_male_idx <- which.min(male_pca_scores[, 1])               # Repeat for "pc1_neg_female_idx" using "female_pca_scores"
                                                                  # Repeat for male and female PC2 extremes using "[, 2]"
                                                                  # Retrieve coordinates for the extreme PC1 and PC2 values (both positive and negative)
pc1_pos_male_shape <- gpa_result$coords[, , pc1_pos_male_idx]     # Repeat for "pc1_neg_male_shape"
pc1_pos_female_shape <- gpa_result$coords[, , pc1_pos_female_idx] # Repeat for "pc1_neg_female_shape"
pc2_pos_male_shape <- gpa_result$coords[, , pc2_pos_male_idx]     # Repeat for "pc2_neg_male_shape"
pc2_pos_female_shape <- gpa_result$coords[, , pc2_pos_female_idx] # Repeat for "pc2_neg_female_shape"
### === Calculate mean shape === ###
mean_shape <- apply(gpa_result$coords, c(1, 2), mean)
### === Code for XY and YZ view === ###
plot_wireframe_custom_yz <- function(shape_data, mean_shape, wireframe, group_color, title) {
  shape_data <- as.data.frame(shape_data)
  mean_shape <- as.data.frame(mean_shape)
  # Set up plot window for YZ view (coordinates 2 and 3. XY coordinates 1 and 2)
  plot(1, type = "n", xlim = range(shape_data$Y), ylim = range(shape_data$Z),
       xlab = "", ylab = "", axes = FALSE, main = title)  # Remove axes and labels
  # Plot the wireframe for the shape
  for (i in 1:nrow(wireframe)) {
    segments(shape_data[wireframe[i, 1], 2], shape_data[wireframe[i, 1], 3], 
             shape_data[wireframe[i, 2], 2], shape_data[wireframe[i, 2], 3], 
             col = group_color, lwd = 2)
  }
  # Plot the mean shape wireframe
  for (i in 1:nrow(wireframe)) {
    segments(mean_shape[wireframe[i, 1], 2], mean_shape[wireframe[i, 1], 3], 
             mean_shape[wireframe[i, 2], 2], mean_shape[wireframe[i, 2], 3], 
             col = "black", lwd = 2, lty = 2)  # Dashed line for mean shape
  }
  # Add landmark dots (points) for both shape data and mean shape
  points(shape_data$Y, shape_data$Z, pch = 16, col = group_color, cex = 1.5)  # Landmark dots for the shape
  points(mean_shape$Y, mean_shape$Z, pch = 16, col = "black", cex = 1.5)      # Landmark dots for the mean shape
  # Add landmark numbers to the plot
  for (i in 1:nrow(shape_data)) {
    text(shape_data[i, 2], shape_data[i, 3], labels = i, pos = 4, cex = 0.7, col = "black")
  }
  
  
}
# Plot for PC1 Negative Male vs Mean (YZ View)
plot_wireframe_custom_yz(pc1_neg_male_shape, mean_shape, wireframe_connections, 
                         group_color = "#4FC3F7", 
                         title = "PC1 Negative Male vs Mean (YZ View)")        # Repeat for "pc1_pos_male_shape"
# Plot for PC1 Negative Female vs Mean (YZ View)
plot_wireframe_custom_yz(pc1_neg_female_shape, mean_shape, wireframe_connections, 
                         group_color = "#FF8A9A", 
                         title = "PC1 Negative Female vs Mean (YZ View)")      # Repeat for "pc1_pos_female_shape"
# Plot for PC2 Negative Male vs Mean (YZ View)
plot_wireframe_custom_yz(pc2_neg_male_shape, mean_shape, wireframe_connections, 
                         group_color = "#4FC3F7", 
                         title = "PC2 Negative Male vs Mean (YZ View)")        # Repeat for "pc2_pos_male_shape"
# Plot for PC2 Negative Female vs Mean (YZ View)
plot_wireframe_custom_yz(pc2_neg_female_shape, mean_shape, wireframe_connections, 
                         group_color = "#FF8A9A", 
                         title = "PC2 Negative Female vs Mean (YZ View)")      # Repeat for "pc2_pos_female_shape"

# repeat wireframes for XY (basal/inferior) view
# Plot male vs female wireframes together 
plot_wireframe_mf_yz <- function(male_shape, female_shape, wireframe, male_color, female_color, title) {
  male_shape <- as.data.frame(male_shape)
  female_shape <- as.data.frame(female_shape)
  
  # Set up YZ view (Y = column 2, Z = column 3)
  plot(1, type = "n", 
       xlim = range(c(male_shape$Y, female_shape$Y)), 
       ylim = range(c(male_shape$Z, female_shape$Z)),
       xlab = "", ylab = "", axes = FALSE, main = title)
  
  # Plot wireframes
  for (i in 1:nrow(wireframe)) {
    segments(male_shape[wireframe[i, 1], 2], male_shape[wireframe[i, 1], 3], 
             male_shape[wireframe[i, 2], 2], male_shape[wireframe[i, 2], 3], 
             col = male_color, lwd = 2)
    
    segments(female_shape[wireframe[i, 1], 2], female_shape[wireframe[i, 1], 3], 
             female_shape[wireframe[i, 2], 2], female_shape[wireframe[i, 2], 3], 
             col = female_color, lwd = 2)
  }
  
  # Plot landmarks
  points(male_shape$Y, male_shape$Z, pch = 16, col = male_color, cex = 1.5)
  points(female_shape$Y, female_shape$Z, pch = 16, col = female_color, cex = 1.5)
}
# Plot for PC1 Negative
plot_wireframe_mf_yz(pc1_neg_male_shape, pc1_neg_female_shape, wireframe_connections,
                     male_color = "#4FC3F7", female_color = "#FF8A9A",
                     title = "PC1 Negative: Male vs Female (YZ View)") 
# Repeat for "pc1_pos_male_shape, pc1_pos_female_shape"
# Plot for PC2 Negative
plot_wireframe_mf_yz(pc2_neg_male_shape, pc2_neg_female_shape, wireframe_connections,
                     male_color = "#4FC3F7", female_color = "#FF8A9A",
                     title = "PC2 Negative: Male vs Female (YZ View)") 
# Repeat for "pc2_pos_male_shape, pc2_pos_female_shape"


plot_wireframe_mf_xy <- function(male_shape, female_shape, wireframe, male_color, female_color, title) {
  male_shape <- as.data.frame(male_shape)
  female_shape <- as.data.frame(female_shape)
  
  # Set up XY view (X = column 1, Y = column 2)
  plot(1, type = "n", 
       xlim = range(c(male_shape$X, female_shape$X)), 
       ylim = range(c(male_shape$Y, female_shape$Y)),
       xlab = "", ylab = "", axes = FALSE, main = title)
  # Plot wireframes
  for (i in 1:nrow(wireframe)) {
    segments(male_shape[wireframe[i, 1], 1], male_shape[wireframe[i, 1], 2], 
             male_shape[wireframe[i, 2], 1], male_shape[wireframe[i, 2], 2], 
             col = male_color, lwd = 2)
    
    segments(female_shape[wireframe[i, 1], 1], female_shape[wireframe[i, 1], 2], 
             female_shape[wireframe[i, 2], 1], female_shape[wireframe[i, 2], 2], 
             col = female_color, lwd = 2)
  }
  
  # Plot landmarks
  points(male_shape$X, male_shape$Y, pch = 16, col = male_color, cex = 1.5)
  points(female_shape$X, female_shape$Y, pch = 16, col = female_color, cex = 1.5)
}
# Plot for PC1 Negative
plot_wireframe_mf_xy(pc1_neg_male_shape, pc1_neg_female_shape, wireframe_connections,
                     male_color = "#4FC3F7", female_color = "#FF8A9A",
                     title = "PC1 Negative: Male vs Female (XY View)")         # Repeat for "pc1_pos_male_shape, pc1_pos_female_shape"
# Plot for PC2 Negative
plot_wireframe_mf_xy(pc2_neg_male_shape, pc2_neg_female_shape, wireframe_connections,
                     male_color = "#4FC3F7", female_color = "#FF8A9A",
                     title = "PC2 Negative: Male vs Female (XY View)")         # Repeat for "pc2_pos_male_shape, pc2_pos_female_shape"
### === Linear model results split by sex ====================================== ###
lm_results <- function(df, var) {
  males <- lm(Stature_cm ~ get(var), data = df %>% filter(sex_code == "Male"))
  females <- lm(Stature_cm ~ get(var), data = df %>% filter(sex_code == "Female"))
  list(male = summary(males), female = summary(females))
}
### === Calculate clivus length (ho-ba) and anterior intercondylar distances (Lfob-Rfob) === ###
calc_clivus_length <- function(coords) {
  apply(coords, 3, function(sp) sqrt(sum((sp[1, ] - sp[3, ])^2)))
}
calc_intercondylar <- function(coords) {
  apply(coords, 3, function(sp) sqrt(sum((sp[21, ] - sp[22, ])^2)))
}

### === Summary statistics by sex === ###
summary_stats <- function(df, var) {
  df %>%
    group_by(Sex) %>%
    summarise(
      Mean = mean(.data[[var]], na.rm = TRUE),
      SD   = sd(.data[[var]], na.rm = TRUE),
      N    = n(),
      .groups = "drop"
    )
}
### === Clivus anf AIC dataframes === ###
clivus_df <- data.frame(
  Specimen = dimnames(landmarks)[[3]],
  Sex = rep(c("Male", "Female"), c(length(male_files), length(female_files))),
  Clivus_Length_mm = calc_clivus_length(landmarks)
)
intercondylar_df <- data.frame(
  Specimen = dimnames(landmarks)[[3]],
  Sex = rep(c("Male", "Female"), c(length(male_files), length(female_files))),
  Intercondylar_Distance = calc_intercondylar(landmarks)
)
### === Calculate CBA2 external cranial base flexion angle (ho-sphba-ba) ======= ###
calc_flexion_angle <- function(coords) {
  apply(coords, 3, function(sp) {
    v1 <- sp[1, ] - sp[2, ]
    v2 <- sp[3, ] - sp[2, ]
    acos(sum(v1 * v2) / (sqrt(sum(v1^2)) * sqrt(sum(v2^2)))) * 180 / pi
  })
}
### === Calculate CBA1 (na-sphba-ba external cranial base flexion angle (na-sphba-ba) === ###
calc_angle <- function(A, B, C) {
  BA <- A - B
  BC <- C - B
  cos_theta <- sum(BA * BC) / (sqrt(sum(BA^2)) * sqrt(sum(BC^2)))
  angle_rad <- acos(pmin(pmax(cos_theta, -1), 1)) # clamp to [-1,1]
  angle_deg <- angle_rad * 180 / pi
  return(angle_deg)
}
### === In this workflow - nasion was collected separately, and subsequently must
### === be incorporated into the other landmark files to calculate CBA1 === ###
alpaca_dir <- "YOUR_FILES"
new_dir <- "YOUR_FILES"
alpaca_files <- list.files(alpaca_dir, pattern = "\\.txt$", full.names = TRUE)
new_files <- list.files(new_dir, pattern = "\\.txt$", full.names = TRUE)
get_id <- function(filepath) tools::file_path_sans_ext(basename(filepath))
alpaca_ids <- sapply(alpaca_files, get_id)
new_ids <- sapply(new_files, get_id)
common_ids <- intersect(alpaca_ids, new_ids)
cba_results <- data.frame(ID = character(), CBA = numeric(), stringsAsFactors = FALSE)
for(id in common_ids) {
  alpaca_file <- alpaca_files[which(alpaca_ids == id)]
  new_file <- new_files[which(new_ids == id)]
  
  sphenobasion <- read_coords_from_line(alpaca_file, 2)  # Line 2
  basion <- read_coords_from_line(alpaca_file, 3)        # Line 3
  nasion <- read_coords_from_line(new_file, 3)           # Line 3
  
  if(any(is.na(sphenobasion), is.na(basion), is.na(nasion))) {
    warning(paste("Missing coords for", id))
    next
  }
  
  # Calculate CBA (angle at sphenobasion formed by nasion-sphenobasion-basion)
  cba <- calc_angle(nasion, sphenobasion, basion)
  
  cba_results <- rbind(cba_results, data.frame(ID = id, CBA = cba))
}
### === Example of CBA dataframe and statistical text - suitable for plotting === ###
flexion_angles <- calc_flexion_angle(landmarks)
results_df <- data.frame(Specimen = dimnames(landmarks)[[3]],
                         Sex = rep(c("Male", "Female"), c(length(male_files), length(female_files))),
                         Flexion_Angle = flexion_angles)
t_test_res <- t.test(Flexion_Angle ~ Sex, data = results_df)
### === ANTERIOR AND POSTERIOR BASICRANIAL SHAPE ANALYSIS ====================== ###
### === LOAD LANDMARKING PROTOCOL TO DEFINE ANERIOR-POSTERIOR REGIONS AND EXTRACT ARRAYS === ###
landmark_regions <- read.csv("YOUR_FILE.csv", stringsAsFactors = FALSE)
expected_labels <- as.character(landmark_regions$Label)
anterior_labels <- landmark_regions$Label[landmark_regions$Region == "Anterior"]           # Repeat for "posterior_labels"
anterior_array <- abind(lapply(landmarks, \(x) x[anterior_labels, , drop=FALSE]), along=3) # Repeat for "posterior_array" using "posterior_labels"
### === ANTERIOR-POSTERIOR GPA and PCA === ###
gpa_anterior <- gpagen(anterior_array)                                                     # Repeat for "gpa_posterior" using "posterior_array"
pca_anterior <- gm.prcomp(gpa_anterior$coords)                                             # Repeat for "pca_posterior" using "gpa_posterior"
### === CALCULATE VARIANCE === ###
max_pcs <- min(ncol(pca_anterior$x), ncol(pca_posterior$x), 10)
var_anterior <- sum(apply(pca_anterior$x[, 1:max_pcs], 2, var))                            # Repeat using "var_posterior" using "pca_posterior"
### === Permutation test === ###
obs_diff <- abs(var_anterior - var_posterior)
perm_diffs <- replicate(999, {
  idx <- sample(c(rep("A", nrow(pca_anterior$x)), rep("P", nrow(pca_posterior$x))))
  v1 <- sum(apply(rbind(pca_anterior$x, pca_posterior$x)[idx == "A", 1:max_pcs], 2, var))
  v2 <- sum(apply(rbind(pca_anterior$x, pca_posterior$x)[idx == "P", 1:max_pcs], 2, var))
  abs(v1 - v2)
})
p_value <- mean(perm_diffs >= obs_diff)
### === Sex-based variance === ###
get_scores_by_sex <- function(scores, sex, label) sum(apply(scores[sex == label, ], 2, var))
var_ant_male <- get_scores_by_sex(pca_anterior$x[, 1:max_pcs], specimen_sex, "Male")       # Repeat for "var_ant_female"
var_post_male <- get_scores_by_sex(pca_posterior$x[, 1:max_pcs], specimen_sex, "Male")     # Repeat for "var_post_female"
### === Procrustes ANOVA by sex === ###
proc_sex_ant <- procD.lm(gpa_anterior$coords ~ specimen_sex, iter=999)                     # Repeat "proc_sex_post" using "gpa_posterior"
### === MANOVA: sex with region === ###
combined_scores <- rbind(pca_anterior$x[, 1:max_pcs], pca_posterior$x[, 1:max_pcs])
group_df <- data.frame(sex = rep(specimen_sex, 2),
                       region = rep(c("Anterior", "Posterior"), each=length(specimen_sex)))
manova(combined_scores ~ sex * region, data = group_df)
### === Shape regression on ECV === ###
ecv_vector <- ecv_df$ECV[match(dimnames(gpa_anterior$coords)[[3]], ecv_df$deidentified_record_number)]
proc_reg_ant <- procD.lm(gpa_anterior$coords ~ ecv_vector, iter=999)                       # Repeat for "proc_reg_post" using "gpa_posterior"
### === Regression scores === ###
reg_score_ant <- two.d.array(gpa_anterior$coords) %*% proc_reg_ant$lmcoef[-1]              # Repeat for "reg_score_post" 
### === Linear models: PC1 ~ ECV, by region and sex === ###
lm_ant_male   <- lm(PC1_Anterior ~ ECV, data = subset(df_pc1, sex_code == "Male"))         # Repeat for "lm_ant_female"
lm_post_male  <- lm(PC1_Posterior ~ ECV, data = subset(df_pc1, sex_code == "Male"))        # Repeat for "lm_post_female"
#### === Combined sex anterior regression === ###
ant_df <- subset(df_long, Region == "Anterior")
lm_ant <- lm(PC1 ~ log_ECV, data = ant_df)
sum_ant <- summary(lm_ant)
pval_ant <- format.pval(sum_ant$coefficients[2,4], digits = 3)
r2_ant <- round(sum_ant$r.squared, 3)
subtitle_ant <- paste0("p = ", pval_ant, ", R² = ", r2_ant)
### === Combined sex posterior regression === ###
post_df <- subset(df_long, Region == "Posterior")
lm_post <- lm(PC1 ~ log_ECV, data = post_df)
sum_post <- summary(lm_post)
pval_post <- format.pval(sum_post$coefficients[2,4], digits = 3)
r2_post <- round(sum_post$r.squared, 3)
subtitle_post <- paste0("p = ", pval_post, ", R² = ", r2_post)
### === Anterior plot === ###
p_ant <- ggscatter(ant_df, x = "log_ECV", y = "PC1",
                   add = "reg.line",
                   conf.int = TRUE,
                   color = "#CD8162",
                   alpha = 0.7,
                   xlab = "Log(ECV)",
                   ylab = "PC1 (Anterior)") +
  ggtitle("Anterior PC1 vs Log(ECV)",
          subtitle = subtitle_ant) +
  theme(plot.title = element_text(size = 14, face = "plain"),
        plot.subtitle = element_text(size = 12))
### === Posterior plot === ###
p_post <- ggscatter(post_df, x = "log_ECV", y = "PC1",
                    add = "reg.line",
                    conf.int = TRUE,
                    color = "paleturquoise3",
                    alpha = 0.7,
                    xlab = "Log(ECV)",
                    ylab = "PC1 (Posterior)") +
  ggtitle("Posterior PC1 vs Log(ECV)",
          subtitle = subtitle_post) +
  theme(plot.title = element_text(size = 14, face = "plain"),
        plot.subtitle = element_text(size = 12))
grid.arrange(p_ant, p_post, ncol = 2)
### === WARPED MEAN HEATMAP MODEL ============================================== ### 
Male <- vcgImport("e.g. warpedmalemesh.ply")
Female <- vcgImport("e.g. warpedfemalemesh.ply")
diff_result <- meshDist(Male, Female, symmetric = TRUE)
### === CALCULATE MESH Z-LENGTH (ANTERO-POSTERIOR LENGTH, max-min vertex) ====== ### 
### === Load all meshes === ###
male_dir <- "MALE_FILE_TO_SKULL_MESHES"
female_dir <- "FEMALE_FILE_TO_SKULL_MESHES"
load_mesh_safe <- function(id, sex) {
  folder <- ifelse(sex == "Male", male_dir, female_dir)
  filename <- file.path(folder, paste0(id, "_skull.ply"))
  
  if (!file.exists(filename)) {
    warning(paste("Mesh file not found for ID:", id, "Expected path:", filename))
    return(NULL)
  }
  
  mesh <- tryCatch({
    vcgImport(filename)
  }, error = function(e) {
    warning(paste("Failed to import mesh for ID:", id, "Error:", e$message))
    return(NULL)
  })
  
  return(mesh)
}
meshes <- list()
for (i in seq_along(confirmed_ids)) {
  id <- confirmed_ids[i]
  sex <- as.character(group[i])
  mesh <- load_mesh_safe(id, sex)
  if (!is.null(mesh)) {
    meshes[[id]] <- mesh
  }
}
#### === Separate by sex === ###
male_meshes <- meshes[which(group == "Male" & names(meshes) %in% confirmed_ids[group == "Male"])] # Repeat for "female_meshes"
### === FUNCTION TO COMPUTE MESH LENGTH (Z AXIS) =============================== ###
meshLength <- function(mesh) {
  z_coords <- mesh$vb[3, ]
  total_length <- max(z_coords) - min(z_coords)
  return(total_length)
}
mesh_lengths <- sapply(confirmed_ids, function(id) {
  mesh <- meshes[[id]]
  if (!is.null(mesh)) {
    return(meshLength(mesh))
  } else {
    return(NA)  # If mesh is missing, assign NA
  }
})
### === COLLATE INTO ONE DATAFRAME === ###
mesh_lengths_df <- data.frame(
  Specimen = confirmed_ids,
  Mesh_Z_Length = mesh_lengths,
  Sex = group
)
### === TEST SEX DIFFERENCE === ###
t_test_result <- t.test(Mesh_Z_Length ~ Sex, data = mesh_lengths_df)
### === REPEAT CODE FOR MESH X-WIDTH (MEDIO-LATERAL DIMENSION) ================= ###
meshWidth <- function(mesh) {
  x_coords <- mesh$vb[1, ]
  total_width <- max(x_coords) - min(x_coords)
  return(total_width)
}
### === CALCULATE HOWELLS (1973) MEASUREMENTS USING ADDITIONAL CRANIOFACIAL LANDMARKS === ###
### === Load craniofacial landmark files for Howells measurements ===###
input_dir <- "YOUR_FILE_DIRECTORY"
output_dir <- file.path(dirname(input_dir), "txt_converted")
### === Create output folder if it doesn't exist === ###
if (!dir.exists(output_dir)) dir.create(output_dir)
### === Convert these 3DSlicer .json files to .txt files using the function used earlier === ###
landmark_array <- abind(lapply(landmarks, as.matrix), along = 3)
gpa_result <- gpagen(landmark_array, print.progress = TRUE)
### === List all landmark files === ###
files <- list.files(input_dir, pattern = "\\.mrk\\.json$", full.names = TRUE)
ids <- tools::file_path_sans_ext(basename(files))
### === Function to calculate Euclidean distance between two landmarks by label === ###
calc_distance <- function(coords, lm1, lm2) {
  p1 <- coords[lm1, ]
  p2 <- coords[lm2, ]
  sqrt(sum((p1 - p2)^2))
}
### === Combine all desired Howells (1973) measurements in one list === ###
landmark_pairs <- list(
  GOL = c("4", "6"),    # Glabella – Opisthocranion
  XCB = c("7", "8"),    # Left Euryon – Right Euryon
  NOL = c("3", "6"),    # Nasion – Opisthocranion
  AUB = c("9", "10"),    # Porion – Porion
  ZMB = c("21", "22")   # Zygomaxillare (left-right)
)
# === For bilateral measurements, use the average of both sides ===
calc_bilateral_avg <- function(coords, left_pair, right_pair) {
  if (all(c(left_pair, right_pair) %in% rownames(coords))) {
    left <- calc_distance(coords, left_pair[1], left_pair[2])
    right <- calc_distance(coords, right_pair[1], right_pair[2])
    return(mean(c(left, right)))
  } else {
    return(NA)
  }
}
### === Combine measurements into  results dataframe === ###
results <- data.frame(ID = character(),
                      GOL = numeric(), XCB = numeric(), NOL = numeric(),
                      AUB = numeric(), OBH = numeric(), OBB = numeric(),
                      ZMB = numeric(),
                      stringsAsFactors = FALSE)
### === Loop function to perform Howells measurements on all landmark configurations === ###
for (id in names(landmarks)) {
  coords <- landmarks[[id]]
  row <- list(ID = id)
  
  for (measure in names(landmark_pairs)) {
    pair <- landmark_pairs[[measure]]
    if (all(pair %in% rownames(coords))) {
      row[[measure]] <- calc_distance(coords, pair[1], pair[2])
    } else {
      warning(paste("Missing landmark(s) for", measure, "in ID", id))
      row[[measure]] <- NA
    }
  }
  # === Bilateral measures with averaging ===
  row$OBH <- calc_bilateral_avg(coords, c("11", "17"), c("12", "18"))
  row$OBB <- calc_bilateral_avg(coords, c("13", "15"), c("14", "16"))
  
  results <- rbind(results, row, stringsAsFactors = FALSE)
}
results <- data.frame(ID = character(), stringsAsFactors = FALSE)
### === CALCULATE BBH (BASION-BREGMA HEIGHT) === 
### === Basion was in a separate landmark directory from the Howells measurements of the NMDID meshes === ###
### === Basion and Bregma must first be merged into one dataframe to calculate BBH === ###
### === Set path to the folder with .txt landmark files === ###
bregma_txt_dir <- "YOUR_FILE_DIRECTORY"
txt_files <- list.files(bregma_txt_dir, pattern = "\\.txt$", full.names = TRUE)
### === Initialize list to hold Bregma coordinates (landmark no. 5, 5th line in .txt coordinates) === ###
bregma_list <- list()
for (file in txt_files) {
  lines <- readLines(file)
  
  if (length(lines) >= 5) {
    coords_str <- strsplit(lines[5], split = ",")[[1]]
    coords_num <- as.numeric(coords_str)
    
    if (length(coords_num) == 3 && !any(is.na(coords_num))) {
      specimen_name <- tools::file_path_sans_ext(basename(file))
      bregma_list[[specimen_name]] <- coords_num
    } else {
      warning(paste("Line 5 in", file, "does not contain exactly 3 numeric values"))
    }
  } else {
    warning(paste("File", file, "does not have at least 5 lines"))
  }
}
common_ids <- intersect(names(bregma_list), names(basion_list))
bbh_results <- data.frame(
  ID = character(),
  BBH = numeric(),
  stringsAsFactors = FALSE
)
### === EXTRACT BASION FROM THE OTHER FILE DIRECTORY AND COMBINE WITH BREGMA COORDINATES === ###
### === Repeat the code above for the following directory: Basion is coordinate line 3 === ###
basion_txt_dir <- "YOUR_FILE_DIRECTORY"
txt_files <- list.files(basion_txt_dir, pattern = "\\.txt$", full.names = TRUE)
### === Now we can calculate BBH for each specimen === ###
for (id in common_ids) {
  bregma <- bregma_list[[id]]
  basion <- basion_list[[id]]
  
  if (!any(is.na(bregma)) && !any(is.na(basion))) {
    bbh <- sqrt(sum((bregma - basion)^2))
    bbh_results <- rbind(bbh_results, data.frame(ID = id, BBH = bbh))
  } else {
    warning(paste("Missing or invalid coordinates for", id))
  }
}
write.csv(bbh_results, "YOUR_FILE_PATH.csv", row.names = FALSE)
### === Calculate means and t-test for GOL (repeat code for XCB etc) === ###
means <- data %>%
  group_by(sex_code) %>%
  summarise(mean_GOL = mean(GOL, na.rm = TRUE))
t_test_res <- t.test(GOL ~ sex_code, data = data)
p_val <- t_test_res$p.value
boxplot(GOL ~ sex_code, data = data, main = "GOL by Sex", ylab = "GOL (mm)", xlab = "Sex")
### === Utility Functions === ###
normalize <- function(v) v / sqrt(sum(v^2))

calc_angle <- function(p1, p2, p3) {
  v1 <- p1 - p2
  v2 <- p3 - p2
  cos_theta <- sum(v1 * v2) / (sqrt(sum(v1^2)) * sqrt(sum(v2^2)))
  angle_rad <- acos(pmin(pmax(cos_theta, -1), 1))
  angle_deg <- angle_rad * 180 / pi
  return(angle_deg)
}
read_coords_from_line <- function(file, line) {
  coords <- strsplit(readLines(file)[line], ",")[[1]]
  as.numeric(coords)
}
get_specimen_name <- function(path) {
  tools::file_path_sans_ext(basename(path))
}
### === ORIENTATE LANDMARK CONFIGURATIONS TO FRANKFURT HORIZONTAL (FH) PLANE === ###
### === define landmark labels for FH, it is standard to use the left cranial side === ###
or_l <- "11"  # Orbitale (Left)
po_l <- "9"   # Porion (Left)
po_r <- "10"  # Porion (Right)
### === Function to orient landmark configurations to FH === ###
orient_to_fh <- function(coords, or_l = "11", po_l = "9", po_r = "10") {
  OrL <- coords[or_l, ]
  PoL <- coords[po_l, ]
  PoR <- coords[po_r, ]
  
  x_vec <- normalize(PoR - PoL)
  z_vec <- normalize(cross(PoL - OrL, x_vec))
  y_vec <- cross(z_vec, x_vec)
  
  R <- t(cbind(x_vec, y_vec, z_vec))
  coords_rot <- t(R %*% t(coords))
  rownames(coords_rot) <- rownames(coords)
  return(coords_rot)
}
### === Remeasured flexion angle on FH-oriented array (CBA2: ho-sphba-ba) === ###
calc_flexion_angle_fh <- function(coords_array) {
  apply(coords_array, 3, function(specimen) {
    H  <- specimen["1", ]  # Hormion
    SB <- specimen["2", ]  # Sphenobasion
    Ba <- specimen["3", ]  # Basion
    calc_angle(H, SB, Ba)
  })
}
### === Load and Orient Landmarks === ###
oriented_landmarks <- lapply(landmarks, orient_to_fh)
oriented_array <- abind::abind(lapply(oriented_landmarks, as.matrix), along = 3)
### === CBA2 results === ###
male_specimens <- sapply(male_files, get_specimen_name)                        # Repeat for "female_specimens"
group <- rep(c("Male", "Female"), c(length(male_specimens), length(female_specimens)))
flexion_angles_fh <- calc_flexion_angle_fh(oriented_array)

results_df <- data.frame(
  Specimen = c(male_specimens, female_specimens),
  Sex = factor(group),
  Flexion_Angle_FH = flexion_angles_fh
)

summary_flexion <- results_df %>%
  group_by(Sex) %>%
  summarise(
    Mean = mean(Flexion_Angle_FH, na.rm = TRUE),
    SD = sd(Flexion_Angle_FH, na.rm = TRUE),
    .groups = "drop"
  )
### === FH-Aligned CBA (Per Specimen Files) === ###
alpaca_dir <- "E:/NMDID/ALPACA_single_landmarks/Combined_txt_landmarks"
new_dir <- "E:/NMDID/Craniofacial_landmarks/txt_converted"
get_id <- function(f) tools::file_path_sans_ext(basename(f))

common_ids <- intersect(
  sapply(list.files(alpaca_dir, "\\.txt$", full.names = TRUE), get_id),
  sapply(list.files(new_dir, "\\.txt$", full.names = TRUE), get_id)
)

cba_results <- data.frame(ID = character(), CBA = numeric())
for (id in common_ids) {
  alpaca_file <- file.path(alpaca_dir, paste0(id, ".txt"))
  new_file <- file.path(new_dir, paste0(id, ".txt"))
  
  coords <- do.call(rbind, lapply(readLines(new_file), \(l) as.numeric(strsplit(l, ",")[[1]])))
  if (nrow(coords) < 11 || anyNA(coords)) next
  rownames(coords) <- as.character(seq_len(nrow(coords)))
  
  if (!all(c("9", "10", "11") %in% rownames(coords))) next
  aligned <- orient_to_fh(coords)
  
  nasion <- aligned["3", ]
  sphenobasion <- read_coords_from_line(alpaca_file, 2)
  basion <- read_coords_from_line(alpaca_file, 3)
  if (anyNA(c(nasion, sphenobasion, basion))) next
  
  angle <- calc_angle(nasion, sphenobasion, basion)
  cba_results <- rbind(cba_results, data.frame(ID = id, CBA = angle))
}
write.csv(cba_results, "E:/NMDID/CBA_fh_aligned.csv", row.names = FALSE)
### === Load Metadata and Compute Scaled Metrics === ###
metadata <- read_csv("E:/NMDID/Spreadsheets/METADATA.csv")
metadata$Sex <- as.factor(metadata$Sex)
metadata <- metadata %>%
  mutate(Scaled_Intercondylar = Intercondylar_Distance / XCB)
### === Violin plots function for CBAs, AIC, etc =============================== ###
plot_by_sex <- function(var, title, ylab, file) {
  stats <- metadata %>%
    group_by(Sex) %>%
    summarise(Mean = mean(.data[[var]], na.rm = TRUE),
              SD = sd(.data[[var]], na.rm = TRUE),
              .groups = "drop")
  t_result <- t.test(reformulate("Sex", var), data = metadata)
  subtitle <- sprintf(
    "Male: %.2f ± %.2f | Female: %.2f ± %.2f | %s",
    stats$Mean[stats$Sex == "Male"], stats$SD[stats$Sex == "Male"],
    stats$Mean[stats$Sex == "Female"], stats$SD[stats$Sex == "Female"],
    if (t_result$p.value < 0.001) "p < 0.001" else paste0("p = ", signif(t_result$p.value, 3))
  )
  p <- ggplot(metadata, aes(x = Sex, y = .data[[var]], fill = Sex)) +
    geom_violin(trim = FALSE, alpha = 0.6) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.4) +
    scale_fill_manual(values = c("Male" = "#4FC3F7", "Female" = "#FF8A9A")) +
    labs(title = title, subtitle = subtitle, x = "Sex", y = ylab) +
    theme_minimal() + theme(legend.position = "none")
  ggsave(file, p, width = 7, height = 5, dpi = 300)
}
### === Generate Plots === ###
plot_by_sex("Scaled_Intercondylar", "Scaled Intercondylar Distance", "Ratio", "")
plot_by_sex("CBA..Na.Sphba.Ba.", "CBA1 (na–sphba–ba)", "CBA1 (degrees)", "")
plot_by_sex("CBA1_FH", "CBA1 aligned to FH (na–sphba–ba)", "CBA1 (degrees)", "")
plot_by_sex("CBA..Ho.Sphba.Ba.", "CBA2 (ho–sphba–ba)", "CBA2 (degrees)", "")
### === FURTHER LINEAR REGRESSIONS PLOTTING FUNCTION ========================== ###
plot_lm_by_sex <- function(data, response, predictor, resp_label, pred_label, save_path) {
  male <- subset(data, Sex == "Male")
  female <- subset(data, Sex == "Female")
  
  lm_m <- lm(as.formula(paste(response, "~", predictor)), data = male)
  lm_f <- lm(as.formula(paste(response, "~", predictor)), data = female)
  
  r2_m <- summary(lm_m)$r.squared
  r2_f <- summary(lm_f)$r.squared
  p_m <- summary(lm_m)$coefficients[2, 4]
  p_f <- summary(lm_f)$coefficients[2, 4]
  
  p <- ggplot(data, aes_string(x = predictor, y = response, color = "Sex")) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", aes(group = Sex), se = TRUE, linewidth = 1.1) +
    scale_color_manual(values = c("Male" = "#4FC3F7", "Female" = "#FF8A9A")) +
    labs(
      title = paste(resp_label, "vs", pred_label),
      subtitle = paste0("Male: R²=", round(r2_m,3), ", p=", signif(p_m,3),
                        " | Female: R²=", round(r2_f,3), ", p=", signif(p_f,3)),
      x = pred_label, y = resp_label, color = "Sex"
    ) +
    theme_minimal(base_size = 14)
  
  ggsave(save_path, p, width = 7, height = 5, dpi = 300)
}
### === OTHER ALLOMETRTY PLOTS ================================================= ###
var_pairs <- list(
  c("FM_Area", "GOL"),
  c("FM_Area", "XCB"),
  c("FM_Breadth", "GOL"),
  c("FM_Breadth", "XCB"),
  c("FM_Length", "GOL"),
  c("FM_Length", "XCB"),
  c("FM_Area", "ECV"),
  c("FM_Area", "Stature"),
  c("FM_Length", "ECV"),
  c("FM_Length", "Stature"),
  c("FM_Breadth", "ECV"),
  c("FM_Breadth", "Stature")
)

for (pair in var_pairs) {
  response <- pair[1]
  predictor <- pair[2]
  save_path <- paste0("E:/NMDID/", tolower(response), "_vs_", tolower(predictor), ".png")
  plot_lm_by_sex(metadata, response, predictor, response, predictor, save_path)
}
### === ENDOMAKER: CODE TO CALCULATE ECV AND PRODUCE VIRTUAL ENDOCASTS ========= ###
### === Set file path for cranial mesh === ###
file_path <- "YOUR_MODEL.ply"
### === import the skull model === ###
skull.model <- tryCatch({
  vcgImport(file_path)
}, error = function(e) {
  cat("Error importing:", file_path, "\n", e$message, "\n")
  return(NULL)
})
### === GENERATE ENDOCAST AND ECV === ###
if (!is.null(skull.model)) {
  endocast <- tryCatch({
    endomaker(
      mesh = skull.model,
      param1_endo = 1.5,  # Adjusted parameter
      alpha_vol = 70,     # Lowered alpha for better segmentation
      npovs = 100,        # Increased viewpoints for accuracy
      scalendo = 0.6,     # Rescale endocast
      decmesh = 50000,    # Mesh decimation for efficiency
      volume = TRUE,      # Extract ECV values
      save = TRUE         # Do not save endocast model
    )
  }, error = function(e) {
    cat("Error processing:", file_path, "\n", e$message, "\n")
    return(NULL)
  })
  
  if (!is.null(endocast) && !is.null(endocast$volume)) {
    ecv_value <- endocast$volume
    cat("ECV value for", basename(file_path), ":", ecv_value, "\n")
  } else {
    cat("No ECV value for:", file_path, "\n")
  }
}
### === CONVERT AND DECIMATE AN .OBJ CRANIAL MESH TO .PLY ====================== ###
# === This speeds up processing time and reduces storage requirements ===
folder_path <- "YOUR_FILE_PATH"
obj_files <- list.files(folder_path, pattern = "\\.obj$", full.names = TRUE,recursive=TRUE)
# === Loop through each .obj file and convert it to .ply ===
for (obj_file in obj_files) {
  # Read the .obj file and convert it to a mesh object
  mesh <- file2mesh(obj_file)
  decim_mesh<-vcgQEdecim(mesh,tarface = 1000000)
  # Define the output .ply file path
  ply_file <- sub("\\.obj$", "decimated_1M.ply", obj_file)
  
  # Write the mesh object to a .ply file
  vcgPlyWrite(decim_mesh, ply_file)
}


### === CLIVUS/BASILAR PROCESS LENGTH (ho-ba) === ###
### === Collect summary stats from 'metadata' === ###
mean_male <- round(mean(metadata$Clivus_Length[metadata$Sex == "Male"], na.rm = TRUE), 2) # Repeat for "mean_female"
ttest_res <- t.test(Clivus_Length ~ Sex, data = metadata)
### === SCALED CLIVUS/BASILAR PROCESS LENGTH (ho-ba) TO GOL === ###
### === Compute scaled clivus length === ###
metadata <- metadata %>%
  mutate(
    Clivus_to_GOL = Clivus_Length / GOL,
    Sex = as.factor(Sex)
  )
### === Collect summary stats === ###
mean_male <- round(mean(metadata$Clivus_to_GOL[metadata$Sex == "Male"], na.rm = TRUE), 3) # Repeat for "mean_female"
ttest_res <- t.test(Clivus_to_GOL ~ Sex, data = metadata)
### === SCALED AIC DISTANCE (Loca-Roca) TO XCB === ###
### === Compute scaled anterior intercondylar (AIC) distance === ###
metadata <- metadata %>%
  mutate(Scaled_Intercondylar = Intercondylar_Distance / XCB)
### === Collect summary statistics === ###
mean_male <- round(mean(metadata$Scaled_Intercondylar[metadata$Sex == "Male"], na.rm = TRUE), 2) # Repeat for "mean_female"
ttest_ic <- t.test(Scaled_Intercondylar ~ Sex, data = metadata)
### === blank template violin plot for unscaled/scaled clivus length and AIC distance === ###
ggplot(metadata, aes(x = , y = , fill = Sex)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
  geom_jitter(width = 0.1, alpha = 0.4) +
  scale_fill_manual(values = c("Male" = "#4FC3F7", "Female" = "#FF8A9A")) +
  labs(
    title = "",
    subtitle = ,
    x = "",
    y = ""
  ) +
  theme_minimal() +
  theme(legend.position = "")
### === PARTIAL LEAST SQUARES (PLS) ANALYSIS =================================== ###
### === Load basicranial landmark directory ===
skull_base_dir <- "YOUR_FILES.txt"
skull_files <- list.files(path = skull_base_dir, pattern = "\\.txt$", full.names = TRUE)
read_landmarks <- function(filepath) {
  lines <- readLines(filepath)
  coords_list <- strsplit(lines, split = ",")
  coords_mat <- do.call(rbind, lapply(coords_list, function(x) as.numeric(x)))
  return(coords_mat)
}
skull_base_landmarks <- lapply(skull_files, read_landmarks)
# Convert to array (assuming all have same landmarks count)
n_landmarks_skull <- nrow(skull_base_landmarks[[1]])
n_specimens_skull <- length(skull_base_landmarks)
skull_base_array <- array(unlist(skull_base_landmarks), dim = c(n_landmarks_skull, 3, n_specimens_skull))
gpa_skull_base <- gpagen(skull_base_array)
pls_result <- two.b.pls(gpa_result$coords, gpa_skull_base$coords, iter = 999)
### === PLS OF BASICRANIAL SHAPE VS ECV BY SEX ================================= ###
### === Logical indices === ###
male_index <- which(sex_vector == "Male")                                      # Repeat for "female_index"
### === Subset coordinates and ECV === ###
coords_base_male <- gpa_skull_base$coords[,,male_index]                        # Repeat for "coords_female" using "female_index"
ecv_male <- metadata$ECV[male_index]                                           # Repeat for "ecv_female" using "female_index"
pls_male <- two.b.pls(coords_base_male, ecv_male, iter = 999)                  # Repeat for "pls_female" using "coords_female, ecv_female"
### === Extract PLS1 scores for Skull Base vs ECV regression === ###
pls1_shape_male <- pls_male$XScores[, 1]                                       # Repeat for "pls1_shape_female"
### === Combine scores and ECV into one data frame === ###
pls_ecv_df <- data.frame(
  PLS1_Shape = c(pls1_shape_male, pls1_shape_female),
  log_ECV = log(c(ecv_male, ecv_female)),  # log-transform ECV
  Sex = factor(c(rep("Male", length(pls1_shape_male)),
                 rep("Female", length(pls1_shape_female))))
)
### === Linear regressions using PLS1 === ###
lm_m_ecv <- lm(PLS1_Shape ~ log_ECV, data = subset(pls_ecv_df, Sex == "Male")) # Repeat for "lm_f_ecv"
### === PLS OF BASICRANIAL SHAPE VS FACIAL SHAPE BY SEX === ### 
coords_face_male <- gpa_result$coords[,,male_index]                            # Repeat for "coords_face_female"
pls_face_base_male <- two.b.pls(coords_face_male, coords_base_male, iter = 999)# Repeat for "pls_face_base_female"  
### === Extract PLS1 scores basicranial shape vs facial shape regression === ###
pls1_face_male   <- pls_face_base_male$XScores[, 1]                            # Repeat for "pls1_face_female" 
pls1_base_male   <- pls_face_base_male$YScores[, 1]                            # Repeat for "pls1_base_female"
### === Combine scores into one data frame === ###
pls_fb_df <- data.frame(
  Face_PLS1 = c(pls1_face_male, pls1_face_female),
  Base_PLS1 = c(pls1_base_male, pls1_base_female),
  Sex = factor(c(rep("Male", length(pls1_face_male)),
                 rep("Female", length(pls1_face_female))))
)
### === Linear regressions using PLS1 ========================================== ###
lm_m <- lm(Base_PLS1 ~ Face_PLS1, data = subset(pls_fb_df, Sex == "Male"))     # Repeat for "lm_f" using "Female"
### === Calculate face centroid size (CS) from GPA and subset them === ###
facial_cs <- gpa_result$Csize
face_cs_male       <- facial_cs[group == "Male"]                               # Repeat for "face_cs_female"

# Make sure specimen names exist
specimen_names_male   <- dimnames(coords_base_male)[[3]]
specimen_names_female <- dimnames(coords_base_female)[[3]]
# Create CS matrices with names to match shape arrays
face_cs_male_mat <- matrix(face_cs_male, ncol = 1,
                           dimnames = list(specimen_names_male, "CS"))
face_cs_female_mat <- matrix(face_cs_female, ncol = 1,
                             dimnames = list(specimen_names_female, "CS"))
### === Run PLS analysis of basicranial shape vs facial CS (size) === ###
pls_base_vs_cs_male <- two.b.pls(coords_base_male, face_cs_male_mat, iter = 999)# Repeat for "pls_base_vs_cs_female"
### === Log-transform facial CS for plotting === ###
log_cs_male   <- log(as.numeric(face_cs_male))                                 # Repeat for "log_cs_female"
# Combine into a data frame
pls_cs_df <- data.frame(
  Sex       = c(rep("Male", length(pls1_base_male)), rep("Female", length(pls1_base_female))),
  log_CS    = c(log_cs_male, log_cs_female),
  PLS1_Shape = c(pls1_base_male, pls1_base_female)
)
### === Linear regressions using PLS1 === ###
lm_male   <- lm(PLS1_Shape ~ log_CS, data = subset(pls_cs_df, Sex == "Male"))  # Repeat for "lm_female"
### === Plot linear regression scatterplot models with ggplot2 and ggscatter as above === ###




