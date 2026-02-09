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
pacman::p_load(tidyverse, emmeans, tidyr, dplyr, knitr, ggpubr, rstatix, lmerTest, ggdist)
options(knitr.kable.NA = '') # hide NA with knitr function

## ---------------------------------------------------------- #
## 1. load data ####
## ---------------------------------------------------------- #
getwd()

files <- c(
  "V1R2_MLV.csv", "V2R2_MLV.csv", "V3R2_MLV.csv",
  "hV4R2_MLV.csv", "OPAR2_MLV.csv", "PPAR2_MLV.csv",
  "RSCR2_MLV.csv"
)

get_roi <- function(fname) {
  # drop extension
  base <- sub("\\.csv$", "", fname)
  # take everything before "R2_"
  roi  <- sub("R2_.*$", "", base)
  roi
}

read_with_roi <- function(fname) {
  df <- read.csv(fname, stringsAsFactors = FALSE)
  names(df)[1] <- "R2"
  df$ROI <- get_roi(fname)
  df
}

dat_list <- lapply(files, read_with_roi)
names(dat_list) <- get_roi(files)

all_data <- do.call(rbind, dat_list)
all_data$ROI <- factor(
  all_data$ROI,
  levels = c("V1", "V2", "V3", "hV4", "OPA", "PPA", "RSC")
)

## ---------------------------------------------------------- #
## 2. Descriptive ####
## ---------------------------------------------------------- #
min(all_data$R2)
max(all_data$R2)
hist(all_data$R2)

summary_all <- all_data %>%
  summarise(
    n        = n(),                         # total number of rows
    n_R2_pos = sum(R2 > 0, na.rm = TRUE),   # count of R2 > 0
    prop_R2_pos = mean(R2 > 0, na.rm = TRUE), # proportion R2 > 0
    mean_R2_pos = ifelse(n_R2_pos > 0,
                         mean(R2[R2 > 0], na.rm = TRUE),
                         NA_real_)          # mean of R2 where R2>0
  )
summary_all

summary_by_roi <- all_data %>%
  group_by(ROI) %>%
  summarise(
    n        = n(),                         # total number of rows
    n_R2_pos = sum(R2 > 0, na.rm = TRUE),   # count of R2 > 0
    prop_R2_pos = mean(R2 > 0, na.rm = TRUE), # proportion R2 > 0
    mean_R2_pos = ifelse(n_R2_pos > 0,
                         mean(R2[R2 > 0], na.rm = TRUE),
                         NA_real_)          # mean of R2 where R2>0
  )
summary_by_roi

## ---------------------------------------------------------- #
## 3. Plots ####
## ---------------------------------------------------------- #
# Bar plot with custom labels and colors

# "Proportion of positive R2 values by ROI"

ggplot(summary_by_roi,
       aes(x = ROI, y = prop_R2_pos)) +
  geom_col() +
  coord_cartesian(ylim = c(0.2, 0.9)) +
  labs(
    x = "ROI",
    y = "Proportion of R2 > 0"
  ) +
  theme_classic()+
  theme(
    legend.position="none")

# "Mean positive R2 values by ROI"
ggplot(summary_by_roi,
       aes(x = ROI, y = mean_R2_pos)) +
  geom_col() +
  coord_cartesian(ylim = c(0.002, 0.045)) +
  labs(
    x = "ROI",
    y = "Mean R2 (R2 > 0 only)"
  ) +
  theme_classic()+
  theme(
    legend.position="none")


