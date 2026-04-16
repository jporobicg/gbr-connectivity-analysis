# ============================================================================
# SPAWNING HOUR ANALYSIS
# Analysis of spawning time of day (start and end times)
# ============================================================================

library(tidyverse)
library(ggplot2)

# Read and process data
data <- read.csv("SpawningDatabase.csv", header = TRUE)

# Family classification function
classify_family <- function(genus) {
  genus <- as.character(genus)
  if (genus %in% c("Acropora", "Montipora", "Astreopora", "Anacropora", "Alveopora")) return("Acroporidae")
  if (genus %in% c("Astrea", "Caulastraea", "Coelastrea", "Cyphastrea", "Dipsastraea",
                   "Echinopora", "Favites", "Goniastrea", "Leptoria", "Merulina",
                   "Oulophyllia", "Paragoniastrea", "Pectinia", "Platygyra")) return("Merulinidae")
  if (genus %in% c("Porites", "Goniopora")) return("Poritidae")
  if (genus %in% c("Diploastrea")) return("Diploastreidae")
  if (genus %in% c("Agaricia", "Pavona", "Leptoseris", "Gardineroseris")) return("Agariciidae")
  if (genus %in% c("Lobophyllia", "Symphyllia", "Acanthastrea", "Cynarina",
                   "Echinophyllia", "Homophyllia", "Micromussa", "Oxypora", "Scolymia")) return("Lobophylliidae")
  if (genus %in% c("Euphyllia", "Catalaphyllia", "Nemenzophyllia", "Plerogyra")) return("Euphylliidae")
  return("Other")
}

data$Family <- sapply(data$Genus, classify_family)
randall_families <- c("Acroporidae", "Merulinidae", "Poritidae", 
                     "Diploastreidae", "Agariciidae", "Lobophylliidae", "Euphylliidae")
data_filtered <- data %>% filter(Family %in% randall_families)

cbPalette <- c('#f98b65', '#cae0fe', '#1a3c6e', '#894c38', '#f1b435', 
               "#56B4E9", '#516091', '#F67280', '#D6F8B8')

main_families <- c("Acroporidae", "Merulinidae", "Poritidae", "Lobophylliidae", 
                   "Diploastreidae", "Agariciidae", "Euphylliidae")

cat("================================================================================\n")
cat("SPAWNING HOUR ANALYSIS\n")
cat("================================================================================\n\n")

# ============================================================================
# PART 1: START TIME ANALYSIS
# ============================================================================

cat("PART 1: START TIME ANALYSIS\n")
cat(rep("-", 80), "\n\n", sep = "")

start_data <- data_filtered %>%
  filter(!is.na(Start_decimal)) %>%
  filter(Family %in% main_families) %>%
  mutate(
    Hour_round = round(Start_decimal),
    # Convert hour 0 to 24, hour 1 to 25
    Hour_round = ifelse(Hour_round == 0, 24, 
                 ifelse(Hour_round == 1, 25, Hour_round))
  )

# Calculate proportion for each family
family_proportions_start <- start_data %>%
  group_by(Family, Hour_round) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(Family) %>%
  mutate(
    total = sum(n),
    proportion = n / total
  ) %>%
  ungroup()

# Total proportion (all data)
total_proportion_start <- start_data %>%
  group_by(Hour_round) %>%
  summarise(n = n(), .groups = 'drop') %>%
  arrange(Hour_round) %>%
  mutate(total = sum(n), proportion = n / total) %>% 
  select(-n, -total)

# Weighted average (equal family weight)
weighted_avg_start <- family_proportions_start %>%
  group_by(Hour_round) %>%
  summarise(avg_proportion = mean(proportion), .groups = 'drop') %>%
  arrange(Hour_round) %>%
  mutate(
    normalized_proportion = avg_proportion / sum(avg_proportion),
    cumulative = cumsum(normalized_proportion)
  )

# Calculate percentiles
weighted_q10_start <- weighted_avg_start %>% filter(cumulative >= 0.10) %>% slice(1) %>% pull(Hour_round)
weighted_q50_start <- weighted_avg_start %>% filter(cumulative >= 0.50) %>% slice(1) %>% pull(Hour_round)
total_q10_start <- total_proportion_start %>% mutate(cumulative = cumsum(proportion)) %>% filter(cumulative >= 0.10) %>% slice(1) %>% pull(Hour_round)
total_q50_start <- total_proportion_start %>% mutate(cumulative = cumsum(proportion)) %>% filter(cumulative >= 0.50) %>% slice(1) %>% pull(Hour_round)

cat(sprintf("Total 10th percentile: %d:00\n", total_q10_start))
cat(sprintf("Total 50th percentile: %d:00\n", total_q50_start))
cat(sprintf("Weighted 10th percentile: %d:00\n", weighted_q10_start))
cat(sprintf("Weighted 50th percentile: %d:00\n\n", weighted_q50_start))

# Create plot
max_density <- max(density(start_data$Start_decimal)$y)

p_start <- ggplot() +
  geom_density(data = start_data, aes(x = Start_decimal, fill = Family), alpha = 0.3) + 
  scale_fill_manual(values = cbPalette) +
  geom_smooth(data = total_proportion_start, aes(x = Hour_round, y = proportion),
              method = "loess", span = 0.3, se = FALSE, linetype = "dashed",
              color = "#050505", linewidth = 1.5, alpha = 0.9) +
  geom_smooth(data = weighted_avg_start, aes(x = Hour_round, y = normalized_proportion),
              method = "loess", span = 0.3, se = FALSE, linetype = "dashed",
              color = "#1067ca", linewidth = 1.5, alpha = 0.8) +
  geom_vline(xintercept = weighted_q50_start, color = "#e95014", linetype = "dashed", 
             linewidth = 1.3, alpha = 0.8) +
  geom_vline(xintercept = total_q10_start, color = "#050505", linetype = "dashed", 
             linewidth = 1.3, alpha = 0.8) +
  annotate("text", x = total_q10_start, y = 0.48, 
           label = sprintf("Total 10%%\n(%dh)", total_q10_start), 
           color = "#050505", fontface = "bold", size = 3, hjust = -0.1) +
  annotate("text", x = weighted_q50_start, y = 0.48, 
           label = sprintf("Weighted 50%%\n(%dh)", weighted_q50_start), 
           color = "#e95014", fontface = "bold", size = 3, hjust = -0.1) +
  annotate("text", x = 23, y = 0.35, 
           label = "Total Proportion\n(all data)", 
           color = "#050505", fontface = "bold", size = 3.5, lineheight = 0.9) +
  annotate("text", x = 23, y = 0.30, 
           label = "Weighted Proportion\n(equal family weight)", 
           color = "#1067ca", fontface = "bold", size = 3.5, lineheight = 0.9) +
  ggtitle("Spawning Start Time - Family Densities + Proportion Curves") +
  xlab("Hour of Day") +
  scale_x_continuous(breaks = 17:25, limits = c(17, 25)) +
  scale_y_continuous(name = "Density / Proportion") +
  theme_bw() +
  theme(text = element_text(size = 12),
        legend.position = c(0.12, 0.68),
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5))

ggsave("figures/spawning_start_time.png", p_start, width = 14, height = 8, dpi = 300)

# ============================================================================
# PART 2: END TIME ANALYSIS
# ============================================================================

cat("PART 2: END TIME ANALYSIS\n")
cat(rep("-", 80), "\n\n", sep = "")

end_data <- data_filtered %>%
  filter(!is.na(End_decimal)) %>%
  filter(Family %in% main_families) %>%
  mutate(
    Hour_round = round(End_decimal),
    Hour_round = ifelse(Hour_round == 0, 24, 
                 ifelse(Hour_round == 1, 25, Hour_round))
  )

# Calculate proportions
family_proportions_end <- end_data %>%
  group_by(Family, Hour_round) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(Family) %>%
  mutate(total = sum(n), proportion = n / total) %>%
  ungroup()

total_proportion_end <- end_data %>%
  group_by(Hour_round) %>%
  summarise(n = n(), .groups = 'drop') %>%
  arrange(Hour_round) %>%
  mutate(total = sum(n), proportion = n / total) %>% 
  select(-n, -total)

weighted_avg_end <- family_proportions_end %>%
  group_by(Hour_round) %>%
  summarise(avg_proportion = mean(proportion), .groups = 'drop') %>%
  arrange(Hour_round) %>%
  mutate(
    normalized_proportion = avg_proportion / sum(avg_proportion),
    cumulative = cumsum(normalized_proportion)
  )

# Calculate percentiles
weighted_q50_end <- weighted_avg_end %>% filter(cumulative >= 0.50) %>% slice(1) %>% pull(Hour_round)
total_q50_end <- total_proportion_end %>% mutate(cumulative = cumsum(proportion)) %>% filter(cumulative >= 0.50) %>% slice(1) %>% pull(Hour_round)
weighted_q90_end <- weighted_avg_end %>% filter(cumulative >= 0.90) %>% slice(1) %>% pull(Hour_round)
total_q90_end <- total_proportion_end %>% mutate(cumulative = cumsum(proportion)) %>% filter(cumulative >= 0.90) %>% slice(1) %>% pull(Hour_round)

cat(sprintf("Total 50th percentile: %d:00\n", total_q50_end))
cat(sprintf("Total 90th percentile: %d:00\n", total_q90_end))
cat(sprintf("Weighted 50th percentile: %d:00\n", weighted_q50_end))
cat(sprintf("Weighted 90th percentile: %d:00\n\n", weighted_q90_end))

# Create plot
max_density_end <- max(density(end_data$End_decimal)$y)

p_end <- ggplot() +
  geom_density(data = end_data, aes(x = End_decimal, fill = Family), alpha = 0.3) + 
  scale_fill_manual(values = cbPalette) +
  geom_smooth(data = total_proportion_end, aes(x = Hour_round, y = proportion),
              method = "loess", span = 0.3, se = FALSE, linetype = "dashed",
              color = "#050505", linewidth = 1.5, alpha = 0.9) +
  geom_smooth(data = weighted_avg_end, aes(x = Hour_round, y = normalized_proportion),
              method = "loess", span = 0.3, se = FALSE, linetype = "dashed",
              color = "#1067ca", linewidth = 1.5, alpha = 0.8) +
  geom_vline(xintercept = weighted_q90_end, color = "#e95014", linetype = "dashed", 
             linewidth = 1.3, alpha = 0.8) +
  geom_vline(xintercept = total_q50_end, color = "#050505", linetype = "dashed", 
             linewidth = 1.3, alpha = 0.8) +
  annotate("text", x = total_q50_end, y = 0.58, 
           label = sprintf("Total 50%%\n(%dh)", total_q50_end), 
           color = "#050505", fontface = "bold", size = 3, hjust = -0.1) +
  annotate("text", x = weighted_q90_end, y = 0.58, 
           label = sprintf("Weighted 90%%\n(%dh)", weighted_q90_end), 
           color = "#e95014", fontface = "bold", size = 3, hjust = -0.1) +
  annotate("text", x = 23, y = 0.45, 
           label = "Total Proportion\n(all data)", 
           color = "#050505", fontface = "bold", size = 3.5, lineheight = 0.9) +
  annotate("text", x = 23, y = 0.40, 
           label = "Weighted Proportion\n(equal family weight)", 
           color = "#1067ca", fontface = "bold", size = 3.5, lineheight = 0.9) +
  ggtitle("Spawning End Time - Family Densities + Proportion Curves") +
  xlab("Hour of Day") +
  scale_x_continuous(breaks = 17:25, limits = c(17, 25)) +
  scale_y_continuous(name = "Density / Proportion") +
  theme_bw() +
  theme(text = element_text(size = 12),
        legend.position = c(0.12, 0.68),
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5))

ggsave("figures/spawning_end_time.png", p_end, width = 14, height = 8, dpi = 300)

# ============================================================================
# PART 3: COMBINED START-END COMPARISON
# ============================================================================

cat("PART 3: SPAWNING DURATION ANALYSIS\n")
cat(rep("-", 80), "\n\n", sep = "")

data_both <- data_filtered %>%
  filter(!is.na(Start_decimal) & !is.na(End_decimal)) %>%
  mutate(
    Duration = End_decimal - Start_decimal,
    Duration_adj = ifelse(Duration < 0, Duration + 24, Duration)
  )

cat(sprintf("Observations with both start and end: %d\n", nrow(data_both)))
cat(sprintf("Median spawning duration: %.2f hours (%.0f minutes)\n", 
            median(data_both$Duration_adj), median(data_both$Duration_adj) * 60))
cat(sprintf("Mean spawning duration: %.2f hours (%.0f minutes)\n\n", 
            mean(data_both$Duration_adj), mean(data_both$Duration_adj) * 60))

# ============================================================================
# PART 4: PROPORTION SPAWNING BY HOUR (20, 21, 22, 23, 24)
# ============================================================================

cat("PART 4: PROPORTION SPAWNING BY HOUR (20:00, 21:00, 22:00, 23:00, 24:00)\n")
cat(rep("-", 80), "\n\n", sep = "")

# Calculate proportion at each hour (20, 21, 22, 23, 24) for start times
hours_20_24 <- c(20, 21, 22, 23, 24)

# By family - start times (proportion at each hour)
# First, ensure we have all families and all hours
all_families_hours <- expand_grid(Family = main_families, Hour_round = hours_20_24)

family_hourly_start <- family_proportions_start %>%
  filter(Hour_round %in% hours_20_24) %>%
  select(Family, Hour_round, proportion) %>%
  # Join with all combinations to fill missing hours with 0
  right_join(all_families_hours, by = c("Family", "Hour_round")) %>%
  mutate(proportion = replace_na(proportion, 0)) %>%
  # Pivot to wide format
  pivot_wider(names_from = Hour_round, values_from = proportion, names_prefix = "Hour_") %>%
  # Reorder columns
  select(Family, Hour_20, Hour_21, Hour_22, Hour_23, Hour_24) %>%
  arrange(Family)

# Total proportion at each hour (all data combined)
all_hours_df <- data.frame(Hour_round = hours_20_24)

total_hourly_start <- total_proportion_start %>%
  filter(Hour_round %in% hours_20_24) %>%
  select(Hour_round, proportion) %>%
  # Join with all hours to fill missing with 0
  right_join(all_hours_df, by = "Hour_round") %>%
  mutate(proportion = replace_na(proportion, 0)) %>%
  # Pivot to wide format
  pivot_wider(names_from = Hour_round, values_from = proportion, names_prefix = "Hour_")

# Weighted proportion at each hour (equal family weight)
weighted_hourly_start <- weighted_avg_start %>%
  filter(Hour_round %in% hours_20_24) %>%
  select(Hour_round, normalized_proportion) %>%
  # Join with all hours to fill missing with 0
  right_join(all_hours_df, by = "Hour_round") %>%
  mutate(normalized_proportion = replace_na(normalized_proportion, 0)) %>%
  # Pivot to wide format
  pivot_wider(names_from = Hour_round, values_from = normalized_proportion, names_prefix = "Hour_")

# Print results
cat("START TIME - Proportion by Hour:\n\n")
cat("Total (all data):\n")
if (nrow(total_hourly_start) > 0) {
  cat(sprintf("  Hour 20: %.4f (%.2f%%)\n", total_hourly_start$Hour_20[1], total_hourly_start$Hour_20[1] * 100))
  cat(sprintf("  Hour 21: %.4f (%.2f%%)\n", total_hourly_start$Hour_21[1], total_hourly_start$Hour_21[1] * 100))
  cat(sprintf("  Hour 22: %.4f (%.2f%%)\n", total_hourly_start$Hour_22[1], total_hourly_start$Hour_22[1] * 100))
  cat(sprintf("  Hour 23: %.4f (%.2f%%)\n", total_hourly_start$Hour_23[1], total_hourly_start$Hour_23[1] * 100))
  cat(sprintf("  Hour 24: %.4f (%.2f%%)\n", total_hourly_start$Hour_24[1], total_hourly_start$Hour_24[1] * 100))
}
cat("\n")

cat("Weighted (equal family):\n")
if (nrow(weighted_hourly_start) > 0) {
  cat(sprintf("  Hour 20: %.4f (%.2f%%)\n", weighted_hourly_start$Hour_20[1], weighted_hourly_start$Hour_20[1] * 100))
  cat(sprintf("  Hour 21: %.4f (%.2f%%)\n", weighted_hourly_start$Hour_21[1], weighted_hourly_start$Hour_21[1] * 100))
  cat(sprintf("  Hour 22: %.4f (%.2f%%)\n", weighted_hourly_start$Hour_22[1], weighted_hourly_start$Hour_22[1] * 100))
  cat(sprintf("  Hour 23: %.4f (%.2f%%)\n", weighted_hourly_start$Hour_23[1], weighted_hourly_start$Hour_23[1] * 100))
  cat(sprintf("  Hour 24: %.4f (%.2f%%)\n", weighted_hourly_start$Hour_24[1], weighted_hourly_start$Hour_24[1] * 100))
}
cat("\n")

cat("By Family:\n")
for (i in 1:nrow(family_hourly_start)) {
  fam <- family_hourly_start$Family[i]
  cat(sprintf("  %s:\n", fam))
  cat(sprintf("    Hour 20: %.4f (%.2f%%)\n", family_hourly_start$Hour_20[i], family_hourly_start$Hour_20[i] * 100))
  cat(sprintf("    Hour 21: %.4f (%.2f%%)\n", family_hourly_start$Hour_21[i], family_hourly_start$Hour_21[i] * 100))
  cat(sprintf("    Hour 22: %.4f (%.2f%%)\n", family_hourly_start$Hour_22[i], family_hourly_start$Hour_22[i] * 100))
  cat(sprintf("    Hour 23: %.4f (%.2f%%)\n", family_hourly_start$Hour_23[i], family_hourly_start$Hour_23[i] * 100))
  cat(sprintf("    Hour 24: %.4f (%.2f%%)\n", family_hourly_start$Hour_24[i], family_hourly_start$Hour_24[i] * 100))
  cat("\n")
}

# Export results
write.csv(family_hourly_start, "figures/family_hourly_proportions_start.csv", row.names = FALSE)
write.csv(total_hourly_start, "figures/total_hourly_proportions_start.csv", row.names = FALSE)
write.csv(weighted_hourly_start, "figures/weighted_hourly_proportions_start.csv", row.names = FALSE)

# ============================================================================
# EXPORT RESULTS
# ============================================================================

write.csv(total_proportion_start, "figures/start_time_total_proportions.csv", row.names = FALSE)
write.csv(weighted_avg_start, "figures/start_time_weighted_proportions.csv", row.names = FALSE)
write.csv(total_proportion_end, "figures/end_time_total_proportions.csv", row.names = FALSE)
write.csv(weighted_avg_end, "figures/end_time_weighted_proportions.csv", row.names = FALSE)

cat("✓ Figures saved to figures/ directory:\n")
cat("  - spawning_start_time.png\n")
cat("  - spawning_end_time.png\n\n")
cat("✓ Data exported to figures/ directory:\n")
cat("  - proportion_20_24_hours.csv\n")
cat("  - family_proportion_20_24_start.csv\n")
cat("  - family_proportion_20_24_end.csv\n\n")
cat("✓ Analysis complete!\n")

