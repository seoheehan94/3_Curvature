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
R2_mlv <- read_csv("allROI_R2_MLV.csv") %>%
  mutate(method = "contour")

# R2 for filter
R2_filter <- read_csv("allROI_R2_filter.csv") %>%
  mutate(method = "filter")

# Combine
all_R2 <- bind_rows(R2_mlv, R2_filter) %>%
  mutate(
    ROI = factor(ROI, levels = c("V1", "V2", "V3", "hV4", "OPA", "PPA", "RSC")),
    method = factor(method, levels = c("contour", "filter"))
  )

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
## Descriptives for MLV only ---------------------------------------
all_data_mlv <- all_R2 %>% filter(method == "contour")

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

summary_by_roi_mlv <- all_data_mlv %>%
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
all_data_filter <- all_R2 %>% filter(method == "filter")

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
## 3. Statistical tests ####
## ---------------------------------------------------------- #
# Per subject × ROI × method
subj_summary <- all_R2 %>%   
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
      is.na(p_prop_wilcox)         ~ NA_character_,
      p_prop_wilcox < 0.001        ~ "***",
      p_prop_wilcox < 0.01         ~ "**",
      p_prop_wilcox < 0.05         ~ "*",
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
      is.na(p_wilcox)        ~ NA_character_,
      p_wilcox < 0.001       ~ "***",
      p_wilcox < 0.01        ~ "**",
      p_wilcox < 0.05        ~ "*",
      TRUE                       ~ "ns"
    )
  )

wilcox_results_subj


## ---------------------------------------------------------- #
## 4. Plots ####
## ---------------------------------------------------------- #
summary_by_roi_method <- all_R2 %>%
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
  theme_classic(base_size = 14)


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
  coord_cartesian(ylim = c(0.002, 0.055)) +
  geom_text(data = mean_labels,
            aes(x = ROI, y = y_pos, label = sig_wilcox),
            inherit.aes = FALSE,
            vjust = 0) +
  labs(
    x = "ROI",
    y = "Mean R2 (R2 > 0 only)",
    fill = "Method"
  ) +
  theme_classic(base_size = 14)



