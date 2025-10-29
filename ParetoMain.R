source("C:/Users/Ankur/Downloads/ParetoUtils.R")
remove.packages("ggplot2")
remove.packages("dplyr")
remove.packages("tidyr")
install.packages("ggplot2", repos = "https://cloud.r-project.org")
install.packages("dplyr", repos = "https://cloud.r-project.org")
install.packages("tidyr", repos = "https://cloud.r-project.org")
library(ggplot2)
library(dplyr)
library(tidyr)

# setup
base <- "https://raw.githubusercontent.com/AnkurBishayee/Pareto-Front-GBM-Project/main"

P2  <- readRDS(url(paste0(base, "/landscape.Rds")))
P6  <- readRDS(url(paste0(base, "/landscape%20%282%29.Rds")))
P11 <- readRDS(url(paste0(base, "/landscape%20%283%29.Rds")))

PTemporal <- rbind(P2, P6, P11)
PTemporal <- PTemporal[PTemporal$fq | PTemporal$nn, ]
PTemporal$Location <- "Temporal"

P3  <- readRDS(url(paste0(base, "/landscape%20%281%29.Rds")))
P7  <- readRDS(url(paste0(base, "/landscape%20%285%29.Rds")))
P8  <- readRDS(url(paste0(base, "/landscape%20%286%29.Rds")))
P17 <- readRDS(url(paste0(base, "/landscape%20%287%29.Rds")))
P23 <- readRDS(url(paste0(base, "/landscape%20%288%29.Rds")))
P27 <- readRDS(url(paste0(base, "/landscape%20%289%29.Rds")))
P32 <- readRDS(url(paste0(base, "/landscape%20%2810%29.Rds")))

PFrontal <- rbind(P3, P7, P8, P17, P23, P27, P32)
PFrontal <- PFrontal[PFrontal$fq | PFrontal$nn, ]
PFrontal$Location <- "Frontal"

frontal_means <- aggregate(mean ~ k, data = PFrontal, FUN = mean)
temporal_means <- aggregate(mean ~ k, data = PTemporal, FUN = mean)

# frontal-origin intersection
shared_karys_F <- intersect(frontal_means$k, temporal_means$k)
F_origin_df <- merge(
  frontal_means[frontal_means$k %in% shared_karys_F, ],
  temporal_means,
  by = "k",
  suffixes = c("_Glycolysis", "_OxPhos")
)

# temporal-origin intersection
shared_karys_T <- intersect(temporal_means$k, frontal_means$k)
T_origin_df <- merge(
  temporal_means[temporal_means$k %in% shared_karys_T, ],
  frontal_means,
  by = "k",
  suffixes = c("_OxPhos", "_Glycolysis")
)

fitness_mat <- as.matrix(F_origin_df[, c("mean_Glycolysis", "mean_OxPhos")])
pareto_mask <- is_pareto_efficient(fitness_mat)
F_origin_df$pareto_label <- ifelse(pareto_mask, "Pareto", "NonPareto")
pareto_points <- F_origin_df[pareto_mask, ]

################################################################################
# Figure 1A
################################################################################
ggplot(F_origin_df, aes(x = mean_Glycolysis, y = mean_OxPhos)) +
  geom_point(alpha = 0.5) +
  geom_point(data = pareto_points, color = "blue", size = 2.5) +
  geom_path(data = pareto_points[order(pareto_points$mean_Glycolysis), ],
            aes(x = mean_Glycolysis, y = mean_OxPhos),
            color = "blue", linewidth = 1) +
  labs(
    title = "Frontal (Glycolytic) Fitness vs. Temporal (OxPhos) Fitness",
    x = "Glycolytic Fitness",
    y = "OxPhos Fitness"
  ) +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.3))

################################################################################
# Figure 1B
################################################################################

ggplot(F_origin_df, aes(x = mean_Glycolysis, y = mean_OxPhos)) +
  geom_point(alpha = 0.5) +
  geom_point(data = pareto_points, color = "blue", size = 2.5) +
  geom_smooth(data = pareto_points, aes(x = mean_Glycolysis, y = mean_OxPhos),
              method = "loess", span=2, se = FALSE, color = "blue", linewidth = 1.2) +
  labs(
    title = "Frontal (Glycolytic) Fitness vs. Temporal (OxPhos) Fitness with Smoothing",
    x = "Glycolytic Fitness",
    y = "OxPhos Fitness"
  ) +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))

################################################################################
# Figure 1C
################################################################################

F_origin_df$epsilon <- mapply(
  compute_epsilon,
  F_origin_df$mean_Glycolysis,
  F_origin_df$mean_OxPhos,
  MoreArgs = list(pareto_df = pareto_points)
)

# use 75th percentile to identify generalist/specialists
epsilon_threshold <- quantile(F_origin_df$epsilon, 0.75)
F_origin_df$specialist <- F_origin_df$epsilon > epsilon_threshold
F_origin_df$specialist_label <- ifelse(F_origin_df$specialist, "Specialist", "Generalist")
F_origin_df$specialist_label <- factor(F_origin_df$specialist_label, levels = c("Generalist", "Specialist"))

ggplot(F_origin_df, aes(x = mean_Glycolysis, y = mean_OxPhos, color = specialist_label)) +
  geom_point(alpha = 0.5) +
  geom_point(data = pareto_points, color = "blue", size = 2.5) +
  geom_smooth(data = pareto_points, aes(x = mean_Glycolysis, y = mean_OxPhos),
              method = "loess", span=2, se = FALSE, color = "blue", linewidth = 1.2) +
  labs(
    title = "Specialist vs Generalist Classification with Smoothing",
    # subtitle = paste("75th percentile threshold =", round(epsilon_threshold, 3)),
    x = "Glycolytic Fitness",
    y = "OxPhos Fitness",
    color = "Classification"
  ) +
  scale_color_manual(values = c("Generalist" = "black", "Specialist" = "purple2")) +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.4))

specialist_karyotypes <- F_origin_df[F_origin_df$specialist, ]
print(data.frame(Specialist_Karyotype = specialist_karyotypes$k))

################################################################################
# Table 1
View(pareto_points)
################################################################################

################################################################################
# Figures 2A & 2B
################################################################################
F_origin_df$k <- as.character(F_origin_df$k)

# find number of deletion events (copy number <2)
loss_by_chr_df <- F_origin_df %>%
  rowwise() %>%
  mutate(
    loss_vector = list(as.integer(unlist(strsplit(k, "\\."))) < 2)
  ) %>%
  ungroup() %>%
  mutate(sample_id = row_number()) %>%
  unnest_longer(loss_vector, indices_include = TRUE) %>%
  rename(chr = loss_vector_id, loss = loss_vector) %>%
  filter(loss) %>%
  group_by(specialist_label, chr) %>%
  summarise(count = n(), .groups = "drop")

ggplot(loss_by_chr_df, aes(x = as.factor(chr), y = count, fill = specialist_label)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Specialist" = "purple", "Generalist" = "black")) +
  labs(
    x = "Chromosome",
    y = "Number of Samples with Loss",
    title = "Copy Number Losses by Chromosome: Specialist vs Generalist"
  ) +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.55))


# find number of gain events (copy number >2)
gain_by_chr_df <- F_origin_df %>%
  rowwise() %>%
  mutate(
    gain_vector = list(as.integer(unlist(strsplit(k, "\\."))) > 2)
  ) %>%
  ungroup() %>%
  mutate(sample_id = row_number()) %>%
  unnest_longer(gain_vector, indices_include = TRUE) %>%
  rename(chr = gain_vector_id, gain = gain_vector) %>%
  filter(gain) %>%
  group_by(specialist_label, chr) %>%
  summarise(count = n(), .groups = "drop")

ggplot(gain_by_chr_df, aes(x = as.factor(chr), y = count, fill = specialist_label)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Specialist" = "purple", "Generalist" = "black")) +
  labs(
    x = "Chromosome",
    y = "Number of Samples with Gain",
    title = "Copy Number Gains by Chromosome: Specialist vs Generalist"
  ) +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.55))

################################################################################
# Figure 3A
################################################################################
P2$Pt  <- "P2"
P6$Pt  <- "P6"
P11$Pt <- "P11"
PTemporal <- rbind(P2, P6, P11)
PTemporal <- PTemporal[PTemporal$fq | PTemporal$nn, ]
PTemporal$Location <- "Temporal"

P3$Pt  <- "P3"
P7$Pt  <- "P7"
P8$Pt  <- "P8"
P17$Pt <- "P17"
P23$Pt <- "P23"
P27$Pt <- "P27"
P32$Pt <- "P32"
PFrontal <- rbind(P3, P7, P8, P17, P23, P27, P32)
PFrontal <- PFrontal[PFrontal$fq | PFrontal$nn, ]
PFrontal$Location <- "Frontal"

temporal_list <- split(PTemporal, PTemporal$Pt)
frontal_list  <- split(PFrontal, PFrontal$Pt)


all_lines_df <- data.frame()

# loop to create all interpatient fronts
for(temp_pt in names(temporal_list)) {
  for(front_pt in names(frontal_list)) {
    temp_df <- temporal_list[[temp_pt]]
    front_df <- frontal_list[[front_pt]]
    
    temp_agg <- aggregate(mean ~ k, data = temp_df, FUN = mean)
    front_agg <- aggregate(mean ~ k, data = front_df, FUN = mean)
    
    shared_k <- intersect(temp_agg$k, front_agg$k)
    if(length(shared_k) < 2) next  # skip if not enough shared points
    
    temp_shared <- temp_agg[temp_agg$k %in% shared_k, ]
    front_shared <- front_agg[front_agg$k %in% shared_k, ]
    
    comp_df <- merge(temp_shared, front_shared, by = "k", suffixes = c("_Temporal", "_Frontal"))
    
    fit_mat <- as.matrix(comp_df[, c("mean_Frontal", "mean_Temporal")])
    pareto_mask <- is_pareto_efficient(fit_mat)
    comp_df <- comp_df[pareto_mask, ]
    
    
    comp_df$pair <- paste0(temp_pt, "_vs_", front_pt)
    comp_df <- comp_df[order(comp_df$mean_Frontal), ]
    
    all_lines_df <- rbind(all_lines_df, comp_df)
  }
}

# plotting interpatient fronts
ggplot(all_lines_df, aes(x = mean_Frontal, y = mean_Temporal, color = pair, group = pair)) +
  geom_line(alpha = 0.7) +
  geom_point(alpha = 0.5) +
  labs(
    title = "Frontal vs Temporal Patient Pair Pareto Fronts",
    x = "Frontal Fitness",
    y = "Temporal Fitness",
    color = "Patient Pair"
  ) + coord_cartesian() +
  theme(plot.title = element_text(hjust = 0.5))

################################################################################
#Figure 3B
################################################################################
regression_results <- data.frame(
  Pair = character(),
  Slope = numeric(),
  R2 = numeric(),
  Group = character(),
  stringsAsFactors = FALSE
)

results_FF <- run_pairwise_regressions(frontal_list, frontal_list, "Frontal/Frontal")
results_TT <- run_pairwise_regressions(temporal_list, temporal_list, "Frontal/Temporal")
results_FT <- run_pairwise_regressions(temporal_list, frontal_list, "Temporal/Temporal")
regression_results <- bind_rows(results_FF, results_TT, results_FT)

ggplot(regression_results %>% filter(R2 > 0.7), aes(x = Group, y = Slope, fill = Group)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 2, alpha = 0.7) +
  theme_minimal(base_size = 14) +
  labs(title = "Distribution of Pareto Front Slopes by Lobe Group",
       x = "Comparison Group",
       y = "Slope of Pareto Front") +
  scale_fill_manual(values = c("Frontal/Frontal" = "blue4",
                               "Frontal/Temporal" = "blue",
                               "Temporal/Temporal" = "lightblue1")) +
  theme(legend.position = "none") + theme(plot.title = element_text(hjust = 0.5))
anova_result <- aov(Slope ~ Group, data = regression_results)
summary(anova_result)

TukeyHSD(anova_result)
