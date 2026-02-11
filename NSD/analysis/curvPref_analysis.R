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

# Curvature preference for MLV (contour)
curv_mlv <- read_csv("allROI_curvPref_MLV.csv") %>%
  mutate(method = "contour")

# Curvature preference for filter
curv_filter <- read_csv("allROI_curvPref_filter.csv") %>%
  mutate(method = "filter")

all_curv <- bind_rows(curv_mlv, curv_filter) %>%
  mutate(
    ROI = factor(ROI, levels = c("V1", "V2", "V3", "hV4", "OPA", "PPA", "RSC")),
    method = factor(method, levels = c("contour", "filter"))
  )

## ---------------------------------------------------------- #
## 2. Descriptive ####
## ---------------------------------------------------------- #


curv_summary <- all_curv %>%
  group_by(ROI, method) %>%
  summarise(
    n_vox    = n(),
    mean_curv = mean(curvPref, na.rm = TRUE),
    median_curv = median(curvPref, na.rm = TRUE),
    .groups = "drop"
  )
curv_summary

curv_hist <- all_curv %>%
  group_by(ROI, method, curvPref) %>%
  summarise(n_vox = n(), .groups = "drop")

curv_hist_prop <- curv_hist %>%
  group_by(ROI, method) %>%
  mutate(prop = n_vox / sum(n_vox)) %>%
  ungroup()

## ---------------------------------------------------------- #
## 3. Statistical tests ####
## ---------------------------------------------------------- #
# Add method-specific voxel index
curv_indexed <- all_curv %>%
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
      is.na(p_spear)      ~ NA_character_,
      p_spear < 0.001     ~ "***",
      p_spear < 0.01      ~ "**",
      p_spear < 0.05      ~ "*",
      TRUE                ~ "ns"
    )
  )
corr_tests
## ---------------------------------------------------------- #
## 4. Plots ####
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






