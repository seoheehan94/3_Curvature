# ---------*----------*---------*---------*---------*---------*---------*---------#
# ------------------------------------------------------------------------------- #
#                               R2                          #
#
#            By Seohee Han
# ------------------------------------------------------------------------------- #
# Date: 2024-09-14
# Environment: R Studio Cloud, Windows 10 / macOS Big Sur
# ------------------------------------------------------------------------------- #
# ---------*----------*---------*---------*---------*---------*---------*---------#

# get ready
rm(list=ls())
set.seed(4228) # for replication

# load packages 
pacman::p_load(tidyverse, readr, emmeans, tidyr, dplyr, knitr, ggpubr, rstatix, lmerTest, ggdist, purrr)
options(knitr.kable.NA = '') # hide NA with knitr function

## ---------------------------------------------------------- #
## 1. load data ####
## ---------------------------------------------------------- #

# R2 for MLV (to be called "contour")
all_mlv <- read_csv("allROI_MLV.csv") %>%
  mutate(method = "contour")

# R2 for filter
all_filter <- read_csv("allROI_filter.csv") %>%
  mutate(method = "filter")

# Combine
all_data <- bind_rows(all_mlv, all_filter) %>%
  mutate(
    ROI = factor(ROI, levels = c("V1", "V2", "V3", "hV4", "OPA", "PPA", "RSC")),
    method = factor(method, levels = c("contour", "filter"))
  )

all_data_pos <- all_data %>%
  filter(R2 > 0)

## ---------------------------------------------------------- #
## 2a. Descriptive - R2 ####
## ---------------------------------------------------------- #
## Descriptives for MLV only ---------------------------------------
all_data_mlv <- all_data %>% filter(method == "contour")

min(all_data_mlv$R2)
max(all_data_mlv$R2)
hist(all_data_mlv$R2, main = "R2 histogram (contour)")

summary_all_mlv <- all_data_mlv %>%
  summarise(
    n          = n(),
    n_R2_pos   = sum(R2 > 0, na.rm = TRUE),
    prop_R2_pos = mean(R2 > 0, na.rm = TRUE),
    mean_R2_pos = ifelse(n_R2_pos > 0,
                         mean(R2[R2 > 0], na.rm = TRUE),
                         NA_real_)
  )
summary_all_mlv

summary_by_roi_mlv <- all_data %>%
  group_by(ROI) %>%
  summarise(
    n          = n(),
    n_R2_pos   = sum(R2 > 0, na.rm = TRUE),
    prop_R2_pos = mean(R2 > 0, na.rm = TRUE),
    mean_R2_pos = ifelse(n_R2_pos > 0,
                         mean(R2[R2 > 0], na.rm = TRUE),
                         NA_real_)
  )
summary_by_roi_mlv

## Descriptives for filter only ------------------------------------
all_data_filter <- all_data %>% filter(method == "filter")

min(all_data_filter$R2)
max(all_data_filter$R2)
hist(all_data_filter$R2, main = "R2 histogram (filter)")

summary_all_filter <- all_data_filter %>%
  summarise(
    n          = n(),
    n_R2_pos   = sum(R2 > 0, na.rm = TRUE),
    prop_R2_pos = mean(R2 > 0, na.rm = TRUE),
    mean_R2_pos = ifelse(n_R2_pos > 0,
                         mean(R2[R2 > 0], na.rm = TRUE),
                         NA_real_)
  )
summary_all_filter

summary_by_roi_filter <- all_data_filter %>%
  group_by(ROI) %>%
  summarise(
    n          = n(),
    n_R2_pos   = sum(R2 > 0, na.rm = TRUE),
    prop_R2_pos = mean(R2 > 0, na.rm = TRUE),
    mean_R2_pos = ifelse(n_R2_pos > 0,
                         mean(R2[R2 > 0], na.rm = TRUE),
                         NA_real_)
  )
summary_by_roi_filter

## ---------------------------------------------------------- #
## 2a. Descriptive - curvPref ####
## ---------------------------------------------------------- #
curv_summary <- all_data_pos %>%
  group_by(ROI, method) %>%
  summarise(
    n_vox    = n(),
    mean_curv = mean(curvPref, na.rm = TRUE),
    median_curv = median(curvPref, na.rm = TRUE),
    .groups = "drop"
  )
curv_summary

curv_hist <- all_data_pos %>%
  group_by(ROI, method, curvPref) %>%
  summarise(n_vox = n(), .groups = "drop")

curv_hist_prop <- curv_hist %>%
  group_by(ROI, method) %>%
  mutate(prop = n_vox / sum(n_vox)) %>%
  ungroup()


## ---------------------------------------------------------- #
## 3a. Statistical tests - R2####
## ---------------------------------------------------------- #
# Per subject × ROI × method
subj_summary <- all_data %>%   
  group_by(subj, ROI, method) %>%
  summarise(
    n_vox        = n(),
    n_R2_pos     = sum(R2 > 0, na.rm = TRUE),
    prop_R2_pos  = mean(R2 > 0, na.rm = TRUE),
    mean_R2_pos  = ifelse(n_R2_pos > 0,
                          mean(R2[R2 > 0], na.rm = TRUE),
                          NA_real_),
    .groups = "drop"
  )

# "Proportion of positive R2 values by ROI"
# Wide format: contour vs filter per subj × ROI
subj_prop_wide <- subj_summary %>%
  select(subj, ROI, method, prop_R2_pos) %>%
  tidyr::pivot_wider(
    id_cols = c(subj, ROI),
    names_from = method,
    values_from = prop_R2_pos
  )

# Paired test per ROI
prop_tests_subj <- subj_prop_wide %>%
  group_by(ROI) %>%
  summarise(
    n_subj = sum(!is.na(contour) & !is.na(filter)),
    mean_prop_contour = mean(contour, na.rm = TRUE),
    mean_prop_filter  = mean(filter,  na.rm = TRUE),
    p_prop_wilcox = if (n_subj > 1)
      wilcox.test(contour, filter, paired = TRUE, exact = FALSE)$p.value
    else NA_real_,
    .groups = "drop"
  ) %>%
  mutate(
    p_prop_wilcox_fdr = p.adjust(p_prop_wilcox, method = "BH"),
    sig_prop          = case_when(
      is.na(p_prop_wilcox_fdr)         ~ NA_character_,
      p_prop_wilcox_fdr < 0.001        ~ "***",
      p_prop_wilcox_fdr < 0.01         ~ "**",
      p_prop_wilcox_fdr < 0.05         ~ "*",
      TRUE                             ~ "ns"
    )
  )

prop_tests_subj



# "Mean positive R2 values by ROI"

# Wide format: contour vs filter mean R2 per subj × ROI
subj_mean_wide <- subj_summary %>%
  select(subj, ROI, method, mean_R2_pos) %>%
  tidyr::pivot_wider(
    id_cols = c(subj, ROI),
    names_from = method,
    values_from = mean_R2_pos
  )

wilcox_results_subj <- subj_mean_wide %>%
  group_by(ROI) %>%
  summarise(
    n_subj = sum(!is.na(contour) & !is.na(filter)),
    median_contour = median(contour, na.rm = TRUE),
    median_filter  = median(filter,  na.rm = TRUE),
    mean_contour   = mean(contour,   na.rm = TRUE),
    mean_filter    = mean(filter,    na.rm = TRUE),
    p_wilcox = if (n_subj > 1)
      wilcox.test(contour, filter,
                  paired = TRUE, exact = FALSE)$p.value
    else NA_real_,
    .groups = "drop"
  ) %>%
  mutate(
    p_wilcox_fdr = p.adjust(p_wilcox, method = "BH"),
    sig_wilcox   = case_when(
      is.na(p_wilcox_fdr)        ~ NA_character_,
      p_wilcox_fdr < 0.001       ~ "***",
      p_wilcox_fdr < 0.01        ~ "**",
      p_wilcox_fdr < 0.05        ~ "*",
      TRUE                       ~ "ns"
    )
  )

wilcox_results_subj

## ---------------------------------------------------------- #
## 3a. Statistical tests - curvPref####
## ---------------------------------------------------------- #

# Add method-specific voxel index
curv_indexed <- all_data_pos %>%
  group_by(subj, ROI, method) %>%
  mutate(voxel_idx = row_number()) %>%  # index within each subj × ROI × method
  ungroup()

# Split and rename
curv_contour <- curv_indexed %>%
  filter(method == "contour") %>%
  select(subj, ROI, voxel_idx, curv_contour = curvPref)

curv_filter <- curv_indexed %>%
  filter(method == "filter") %>%
  select(subj, ROI, voxel_idx, curv_filter = curvPref)

# Join on subj, ROI, voxel_idx: ensures both methods present for same "voxel"
curv_wide <- curv_contour %>%
  inner_join(curv_filter,
             by = c("subj", "ROI", "voxel_idx"))

# Correlation per subj × ROI
curv_corr_by_subj <- curv_wide %>%
  group_by(subj, ROI) %>%
  summarise(
    n_vox   = n(),
    r_spear = cor(curv_contour, curv_filter, method = "spearman"),
    .groups = "drop"
  ) %>%
  mutate(
    z_spear = atanh(r_spear)
  )

# One-sample t-tests on Fisher z, per ROI
corr_tests <- curv_corr_by_subj %>%
  group_by(ROI) %>%
  summarise(
    n_subj = sum(!is.na(z_spear)),
    
    # Spearman-based group test
    mean_r_spear = mean(r_spear, na.rm = TRUE),
    mean_z_spear = mean(z_spear, na.rm = TRUE),
    p_spear = if (n_subj > 1)
      t.test(z_spear, mu = 0)$p.value
    else NA_real_,
    
    .groups = "drop"
  ) %>%
  mutate(
    p_spear_fdr = p.adjust(p_spear, method = "BH")
  )

corr_tests <- corr_tests %>%
  mutate(
    sig_spear = case_when(
      is.na(p_spear_fdr)      ~ NA_character_,
      p_spear_fdr < 0.001     ~ "***",
      p_spear_fdr < 0.01      ~ "**",
      p_spear_fdr < 0.05      ~ "*",
      TRUE                ~ "ns"
    )
  )
corr_tests

## ---------------------------------------------------------- #
## 4a. Plots - R2 ####
## ---------------------------------------------------------- #
summary_by_roi_method <- all_data %>%
  group_by(ROI, method) %>%
  summarise(
    n            = n(),
    n_R2_pos     = sum(R2 > 0, na.rm = TRUE),
    prop_R2_pos  = mean(R2 > 0, na.rm = TRUE),
    se_prop_R2   = sqrt(prop_R2_pos * (1 - prop_R2_pos) / n),
    mean_R2_pos  = ifelse(n_R2_pos > 0,
                          mean(R2[R2 > 0], na.rm = TRUE),
                          NA_real_),
    se_mean_R2   = ifelse(n_R2_pos > 1,
                          sd(R2[R2 > 0], na.rm = TRUE) / sqrt(n_R2_pos),
                          NA_real_),
    .groups = "drop"
  )

# "Proportion of positive R2 values by ROI"
prop_labels <- prop_tests_subj %>%
  mutate(
    # place stars a bit above the voxel-level max bar
    y_pos = pmax(mean_prop_contour, mean_prop_filter) + 0.03
  ) %>%
  select(ROI, sig_prop, y_pos)

dodge <- position_dodge(width = 0.7)

ggplot(summary_by_roi_method,
       aes(x = ROI, y = prop_R2_pos, fill = method)) +
  geom_col(position = dodge, width = 0.6) +
  geom_errorbar(
    aes(ymin = prop_R2_pos - se_prop_R2,
        ymax = prop_R2_pos + se_prop_R2),
    position = dodge,
    width = 0.2
  ) +
  coord_cartesian(ylim = c(0.2, 0.95)) +
  geom_text(data = prop_labels,
            aes(x = ROI, y = y_pos, label = sig_prop),
            inherit.aes = FALSE,
            vjust = 0) +
  labs(
    x = "ROI",
    y = "Proportion of R2 > 0",
    fill = "Method"
  ) +
  theme_classic(base_size = 18)


# "Mean positive R2 values by ROI"
mean_labels <- wilcox_results_subj %>%
  mutate(
    # a bit above the larger subject-level mean
    y_pos = pmax(mean_contour, mean_filter) + 0.002
  ) %>%
  select(ROI, sig_wilcox, y_pos)

dodge <- position_dodge(width = 0.7)

ggplot(summary_by_roi_method,
       aes(x = ROI, y = mean_R2_pos, fill = method)) +
  geom_col(position = dodge, width = 0.6) +
  geom_errorbar(
    aes(ymin = mean_R2_pos - se_mean_R2,
        ymax = mean_R2_pos + se_mean_R2),
    position = dodge,
    width = 0.2
  ) +
  coord_cartesian(ylim = c(0.002, 0.06)) +
  geom_text(data = mean_labels,
            aes(x = ROI, y = y_pos, label = sig_wilcox),
            inherit.aes = FALSE,
            vjust = 0) +
  labs(
    x = "ROI",
    y = "Mean R2 (R2 > 0 only)",
    fill = "Method"
  ) +
  theme_classic(base_size = 18)

## ---------------------------------------------------------- #
## 4a. Plots - curvPref ####
## ---------------------------------------------------------- #
dodge <- position_dodge(width = 0.7)

# Early visual ROIs: V1–hV4
curv_hist_prop_ev <- curv_hist_prop %>%
  filter(ROI %in% c("V1", "V2", "V3", "hV4"))

ggplot(curv_hist_prop_ev,
       aes(x = factor(curvPref), y = prop, fill = method)) +
  geom_col(position = dodge, width = 0.6) +
  facet_wrap(~ ROI, nrow = 1) +   # no scales = "free_y" → shared y axis
  scale_y_continuous(
    limits = c(0, 0.45),           # adjust to your max
    breaks = seq(0, 0.45, by = 0.1),
    expand = expansion(mult = c(0, 0.05))  # no gap at 0
  ) +
  labs(
    x = "Preferred curvature level",
    y = "Proportion of voxels",
    fill = "Method"
  ) +
  theme_classic(base_size = 16) +
  theme(
    strip.background = element_blank(),
    panel.spacing.x = unit(0.8, "lines"),
    plot.margin = margin(t = 5, r = 5, b = 10, l = 5)
  )


# Scene ROIs: OPA, PPA, RSC
curv_hist_prop_scene <- curv_hist_prop %>%
  filter(ROI %in% c("OPA", "PPA", "RSC"))

ggplot(curv_hist_prop_scene,
       aes(x = factor(curvPref), y = prop, fill = method)) +
  geom_col(position = dodge, width = 0.6) +
  facet_wrap(~ ROI, nrow = 1) +   # no scales = "free_y" → shared y axis
  scale_y_continuous(
    limits = c(0, 0.45),           # adjust to your max
    breaks = seq(0, 0.45, by = 0.1),
    expand = expansion(mult = c(0, 0.05))  # no gap at 0
  ) +
  labs(
    x = "Preferred curvature bin",
    y = "Proportion of voxels",
    fill = "Method"
  ) +
  theme_classic(base_size = 16) +
  theme(
    strip.background = element_blank(),
    panel.spacing.x = unit(0.8, "lines"),
    plot.margin = margin(t = 5, r = 5, b = 10, l = 5)
  )

