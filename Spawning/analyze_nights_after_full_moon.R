# ============================================================================
# NIGHTS AFTER FULL MOON ANALYSIS
# Analysis of spawning timing relative to lunar cycle
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
cat("NIGHTS AFTER FULL MOON ANALYSIS\n")
cat("================================================================================\n\n")

# Process full moon data
data_filtered$DoSRtNFM <- as.numeric(data_filtered$DoSRtNFM)
moon_data <- data_filtered %>%
  filter(!is.na(DoSRtNFM)) %>%
  mutate(DoSRtNFM = if_else(DoSRtNFM < 0, 0, DoSRtNFM)) %>%
  filter(Family %in% main_families)

# ============================================================================
# ANALYSIS 1: PROPORTIONS BY FAMILY
# ============================================================================

cat("ANALYSIS 1: FAMILY-LEVEL PROPORTIONS\n")
cat(rep("-", 80), "\n\n", sep = "")

# Calculate proportion for each family
family_proportions_moon <- moon_data %>%
  group_by(Family, DoSRtNFM) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(Family) %>%
  mutate(
    total = sum(n),
    proportion = n / total
  ) %>%
  ungroup()

# Total proportion (all data)
total_proportion_moon <- moon_data %>%
  count(DoSRtNFM) %>%
  arrange(DoSRtNFM) %>%
  mutate(
    proportion = n / sum(n),
    cumulative_total = cumsum(proportion)
  )

# Weighted average (equal family weight)
weighted_avg_moon <- family_proportions_moon %>%
  group_by(DoSRtNFM) %>%
  summarise(avg_proportion = mean(proportion), .groups = 'drop') %>%
  arrange(DoSRtNFM) %>%
  mutate(
    normalized_proportion = avg_proportion / sum(avg_proportion),
    cumulative = cumsum(normalized_proportion)
  )

cat("Sum of total proportions (moon):", sum(total_proportion_moon$proportion), "\n")
cat("Sum of weighted proportions (moon):", sum(weighted_avg_moon$normalized_proportion), "\n\n")

# Calculate percentiles
weighted_q10_moon <- weighted_avg_moon %>% filter(cumulative >= 0.10) %>% slice(1) %>% pull(DoSRtNFM)
weighted_q50_moon <- weighted_avg_moon %>% filter(cumulative >= 0.50) %>% slice(1) %>% pull(DoSRtNFM)
weighted_q90_moon <- weighted_avg_moon %>% filter(cumulative >= 0.90) %>% slice(1) %>% pull(DoSRtNFM)
total_q10_moon <- total_proportion_moon %>% filter(cumulative_total >= 0.10) %>% slice(1) %>% pull(DoSRtNFM)
total_q50_moon <- total_proportion_moon %>% filter(cumulative_total >= 0.50) %>% slice(1) %>% pull(DoSRtNFM)
total_q90_moon <- total_proportion_moon %>% filter(cumulative_total >= 0.90) %>% slice(1) %>% pull(DoSRtNFM)

cat("Percentiles for Nights After Full Moon:\n")
cat(sprintf("  Total 10th: %d nights\n", total_q10_moon))
cat(sprintf("  Total 50th: %d nights\n", total_q50_moon))
cat(sprintf("  Total 90th: %d nights\n", total_q90_moon))
cat(sprintf("  Weighted 10th: %d nights\n", weighted_q10_moon))
cat(sprintf("  Weighted 50th: %d nights\n", weighted_q50_moon))
cat(sprintf("  Weighted 90th: %d nights\n\n", weighted_q90_moon))

# ============================================================================
# PLOT: NIGHTS AFTER FULL MOON
# ============================================================================

p_moon <- ggplot() +
  # Total proportion bars
  geom_col(data = total_proportion_moon,
           aes(x = DoSRtNFM, y = proportion),
           fill = "gray70", alpha = 0.6, width = 0.8) +
  # Weighted proportion line
  geom_line(data = weighted_avg_moon,
            aes(x = DoSRtNFM, y = normalized_proportion),
            color = "#1067ca", linewidth = 1.5, alpha = 0.8) +
  geom_point(data = weighted_avg_moon,
             aes(x = DoSRtNFM, y = normalized_proportion),
             color = "#1067ca", size = 3, alpha = 0.8) +
  # Total cumulative line
  geom_line(data = total_proportion_moon,
            aes(x = DoSRtNFM, y = cumulative_total),
            color = "#050505", linewidth = 2, linetype = "dashed", alpha = 0.9) +
  # Weighted cumulative line
  geom_line(data = weighted_avg_moon,
            aes(x = DoSRtNFM, y = cumulative),
            color = "black", linewidth = 2.5) +
  # Percentile lines
  geom_vline(xintercept = total_q50_moon, color = "#050505", linetype = "dashed",
             linewidth = 1.3, alpha = 0.8) +
  geom_vline(xintercept = weighted_q50_moon, color = "#e95014", linetype = "dashed",
             linewidth = 1.3, alpha = 0.8) +
  geom_vline(xintercept = total_q90_moon, color = "#050505", linetype = "dotted",
             linewidth = 1.3, alpha = 0.8) +
  geom_vline(xintercept = weighted_q90_moon, color = "#e95014", linetype = "dotted",
             linewidth = 1.3, alpha = 0.8) +
  # Annotations
  annotate("text", x = total_q50_moon, y = 0.95,
           label = sprintf("Total 50%%\n(%d nights)", total_q50_moon),
           color = "#050505", fontface = "bold", size = 3, hjust = -0.1) +
  annotate("text", x = weighted_q50_moon, y = 0.88,
           label = sprintf("Weighted 50%%\n(%d nights)", weighted_q50_moon),
           color = "#e95014", fontface = "bold", size = 3, hjust = -0.1) +
  annotate("text", x = 10, y = 0.30,
           label = "Total Proportion (bars)\nWeighted Proportion (blue line)\nBlack solid: Weighted Cumulative\nBlack dashed: Total Cumulative",
           color = "black", fontface = "bold", size = 3.5, hjust = 0, lineheight = 0.9) +
  ggtitle("Spawning by Nights After Full Moon") +
  xlab("Nights After Full Moon") +
  ylab("Proportion / Cumulative") +
  scale_x_continuous(breaks = 0:12, limits = c(-0.5, 12.5)) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(text = element_text(size = 12))

ggsave("figures/spawning_nights_after_full_moon.png", p_moon, width = 14, height = 8, dpi = 300)

# ============================================================================
# ANALYSIS 2: FAMILY-SPECIFIC BREAKDOWN
# ============================================================================

cat("ANALYSIS 2: FAMILY-SPECIFIC BREAKDOWN\n")
cat(rep("-", 80), "\n\n", sep = "")

# Create faceted plot by family
p_moon_family <- ggplot(family_proportions_moon,
                        aes(x = DoSRtNFM, y = proportion, fill = Family)) +
  geom_col(alpha = 0.7) +
  scale_fill_manual(values = cbPalette) +
  facet_wrap(~Family, ncol = 2, scales = "free_y") +
  labs(title = "Spawning by Nights After Full Moon (By Family)",
       x = "Nights After Full Moon",
       y = "Proportion") +
  scale_x_continuous(breaks = seq(0, 12, 2)) +
  theme_bw() +
  theme(text = element_text(size = 10),
        legend.position = "none",
        strip.background = element_rect(fill = "lightgray"))

ggsave("figures/spawning_nights_by_family.png", p_moon_family, width = 12, height = 10, dpi = 300)

# Family statistics
family_stats <- moon_data %>%
  group_by(Family) %>%
  summarise(
    n = n(),
    median_night = median(DoSRtNFM),
    mean_night = mean(DoSRtNFM),
    sd_night = sd(DoSRtNFM),
    .groups = 'drop'
  ) %>%
  arrange(desc(n))

cat("Family Statistics:\n")
print(family_stats)
cat("\n")

write.csv(family_stats, "figures/family_moon_statistics.csv", row.names = FALSE)

# ============================================================================
# ANALYSIS 3: PROPORTION SPAWNING BY NIGHT (DAYS 2-7 AFTER FULL MOON)
# ============================================================================

cat("ANALYSIS 3: PROPORTION SPAWNING BY NIGHT (DAYS 2-7 AFTER FULL MOON)\n")
cat(rep("-", 80), "\n\n", sep = "")

# Focus on nights 2-7 (days 2-7 after full moon)
nights_2_7 <- c(2, 3, 4, 5, 6, 7)

# Create all combinations of families and nights 2-7
all_families_nights_2_7 <- expand_grid(Family = main_families, DoSRtNFM = nights_2_7)

# By family - proportions by night (wide format) for nights 2-7 only
# First, recalculate proportions for each family considering only nights 2-7
family_proportions_moon_2_7 <- moon_data %>%
  filter(DoSRtNFM %in% nights_2_7) %>%
  group_by(Family, DoSRtNFM) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(Family) %>%
  mutate(
    total = sum(n),
    proportion = n / total
  ) %>%
  ungroup()

# Pivot to wide format
family_nightly_moon_2_7 <- family_proportions_moon_2_7 %>%
  select(Family, DoSRtNFM, proportion) %>%
  # Join with all combinations to fill missing nights with 0
  right_join(all_families_nights_2_7, by = c("Family", "DoSRtNFM")) %>%
  mutate(proportion = replace_na(proportion, 0)) %>%
  # Pivot to wide format with column names like Night_2, Night_3, etc.
  pivot_wider(names_from = DoSRtNFM, values_from = proportion, names_prefix = "Night_") %>%
  # Ensure columns are in order: Night_2, Night_3, Night_4, Night_5, Night_6, Night_7
  select(Family, Night_2, Night_3, Night_4, Night_5, Night_6, Night_7) %>%
  arrange(Family)

# Total proportion at each night (all data combined) for nights 2-7
nights_2_7_df <- data.frame(DoSRtNFM = nights_2_7)

total_nightly_moon_2_7 <- moon_data %>%
  filter(DoSRtNFM %in% nights_2_7) %>%
  count(DoSRtNFM) %>%
  mutate(proportion = n / sum(n)) %>%
  select(DoSRtNFM, proportion) %>%
  # Join with all nights to fill missing with 0
  right_join(nights_2_7_df, by = "DoSRtNFM") %>%
  mutate(proportion = replace_na(proportion, 0)) %>%
  # Pivot to wide format
  pivot_wider(names_from = DoSRtNFM, values_from = proportion, names_prefix = "Night_") %>%
  select(Night_2, Night_3, Night_4, Night_5, Night_6, Night_7)

# Weighted proportion at each night (equal family weight) for nights 2-7
weighted_nightly_moon_2_7 <- family_proportions_moon_2_7 %>%
  group_by(DoSRtNFM) %>%
  summarise(avg_proportion = mean(proportion), .groups = 'drop') %>%
  arrange(DoSRtNFM) %>%
  mutate(normalized_proportion = avg_proportion / sum(avg_proportion)) %>%
  select(DoSRtNFM, normalized_proportion) %>%
  # Join with all nights to fill missing with 0
  right_join(nights_2_7_df, by = "DoSRtNFM") %>%
  mutate(normalized_proportion = replace_na(normalized_proportion, 0)) %>%
  # Pivot to wide format
  pivot_wider(names_from = DoSRtNFM, values_from = normalized_proportion, names_prefix = "Night_") %>%
  select(Night_2, Night_3, Night_4, Night_5, Night_6, Night_7)

# Print results
cat("Total (all data) - Proportion by Night (2-7):\n")
for (night in nights_2_7) {
  col_name <- paste0("Night_", night)
  if (col_name %in% names(total_nightly_moon_2_7)) {
    prop <- total_nightly_moon_2_7[[col_name]][1]
    cat(sprintf("  Night %d: %.10f\n", night, prop))
  }
}
cat("\n")

cat("Weighted (equal family) - Proportion by Night (2-7):\n")
for (night in nights_2_7) {
  col_name <- paste0("Night_", night)
  if (col_name %in% names(weighted_nightly_moon_2_7)) {
    prop <- weighted_nightly_moon_2_7[[col_name]][1]
    cat(sprintf("  Night %d: %.10f\n", night, prop))
  }
}
cat("\n")

cat("By Family - Proportion by Night (2-7):\n")
for (i in 1:nrow(family_nightly_moon_2_7)) {
  fam <- family_nightly_moon_2_7$Family[i]
  cat(sprintf("  %s: ", fam))
  props <- c()
  for (night in nights_2_7) {
    col_name <- paste0("Night_", night)
    if (col_name %in% names(family_nightly_moon_2_7)) {
      prop <- family_nightly_moon_2_7[[col_name]][i]
      props <- c(props, sprintf("%.10f", prop))
    }
  }
  cat(paste(props, collapse = ", "))
  cat("\n")
}
cat("\n")

# Export results
write.csv(family_nightly_moon_2_7, "figures/family_nightly_proportions_moon_2_7.csv", row.names = FALSE)
write.csv(total_nightly_moon_2_7, "figures/total_nightly_proportions_moon_2_7.csv", row.names = FALSE)
write.csv(weighted_nightly_moon_2_7, "figures/weighted_nightly_proportions_moon_2_7.csv", row.names = FALSE)

cat("✓ Additional figures saved:\n")
cat("  - spawning_nights_by_family.png\n")
cat("  - family_moon_statistics.csv\n")
cat("  - family_nightly_proportions_moon.csv\n")
cat("  - total_nightly_proportions_moon.csv\n")
cat("  - weighted_nightly_proportions_moon.csv\n\n")
cat("✓ All analyses complete!\n")
