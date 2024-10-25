# Load packages 
library(stats)
library(rstatix)
library(TOSTER)
library(MOTE)
library(tidyverse)

set.seed(21)

# Replication data loading and prep ---------------------
rep_data <- read_csv("knee_valgus_data.csv")
head(rep_data)

# Convert to long dataset

rep_data_long <- rep_data %>%
  gather(key = "condition", value = "valgus_angle", compressive, noncompressive)
head(rep_data_long, 3)

# add differences column to wide dataset

rep_data <- rep_data %>% 
  mutate(differences =  noncompressive - compressive) 

## Replication descriptives  ---------------

rep_desc <- rep_data_long %>%
  group_by(condition) %>%
  summarise(count = n(),
            mean = mean(valgus_angle),
            sd = sd(valgus_angle)) %>%
  mutate(mean_diff = mean(rep_data$differences), 
         sd_diff = sd(rep_data$differences))

## Resolving assumptions  ---------------------------------------

## Distribution check 

ggplot(rep_data_long, aes(valgus_angle)) +
  geom_histogram(color="black", fill="white", 
                 bins = 10)

ggplot(rep_data_long, aes(condition, valgus_angle, color = condition)) +
  geom_boxplot(show.legend = FALSE) +
  theme_minimal()

### Outliers check 

rep_data %>%
  identify_outliers(differences)

### Normality check  

rep_data %>% shapiro_test(differences) 

## Paired t-test  -----------------------------

rep_data_long$condition <- as.factor(rep_data_long$condition)

# R compares conditions alphabetically, I am reordering here to match the original study

rep_data_long$condition <- forcats::fct_relevel(rep_data_long$condition, "noncompressive", "compressive")

#replication_ttest <- t.test(valgus_angle ~ condition, rep_data_long, 
#                            alternative = "two.sided", paired = TRUE, conf.level = 0.95) %>%
#  tidy()
#replication_ttest

replication_ttest <- t.test(rep_data$noncompressive, rep_data$compressive, paired = TRUE)
replication_ttest

### Replication effect size calculation ------

rep_dz <- d.dep.t.diff.t(t = replication_ttest$statistic, n = rep_desc$count[1], a = 0.05)
rep_dz

rep_dav <- d.dep.t.avg(m1 = rep_desc$mean[1], m2 = rep_desc$mean[2], 
                       sd1 = rep_desc$sd[1], sd2 = rep_desc$sd[2],
                       n = rep_desc$count[1], a = 0.05)
rep_dav


# Original study data ------

orig_data <- read_csv("original_data.csv")
head(orig_data)

# Convert to long dataset

orig_data_long <- orig_data %>%
  gather(key = "condition", value = "valgus_angle", compressive, noncompressive)
head(orig_data_long, 3)

# add differences column to wide dataset

orig_data <- orig_data %>% 
  mutate(differences =  compressive - noncompressive) 

## Original descriptives  ---------------

orig_desc <- orig_data_long %>%
  group_by(condition) %>%
  summarise(count = n(),
            mean = mean(valgus_angle),
            sd = sd(valgus_angle)) %>%
  mutate(mean_diff = mean(orig_data$differences), 
         sd_diff = sd(orig_data$differences))
orig_desc

## Resolving assumptions  ---------------------------------------

## Distribution check 

ggplot(orig_data_long, aes(valgus_angle)) +
  geom_histogram(color="black", fill="white", 
                 bins = 10)

ggplot(orig_data_long, aes(condition, valgus_angle, color = condition)) +
  geom_boxplot(show.legend = FALSE) +
  theme_minimal()

### Outliers check 

orig_data %>%
  identify_outliers(differences)

### Normality check  

orig_data %>% shapiro_test(differences) 

## Paired t-test  -----------------------------

#original_ttest <- t.test(valgus_angle ~ condition, orig_data_long, 
#                        alternative = "two.sided", paired = TRUE, conf.level = 0.95) %>%
#  tidy()
#original_ttest

original_ttest <- t.test(orig_data$noncompressive, orig_data$compressive, paired = TRUE)
original_ttest

### Original effect size calculation ------

orig_dz <- d.dep.t.diff.t(t = original_ttest$statistic, n = orig_desc$count[1], a = 0.05)
orig_dz

# Note - the calculated effect size is slightly different to the original (d = 0.58) so will check using Cohen's dav

orig_dav <- d.dep.t.avg(m1 = orig_desc$mean[1], m2 = orig_desc$mean[2], 
                        sd1 = orig_desc$sd[1], sd2 = orig_desc$sd[2],
                        n = orig_desc$count[1], a = 0.05)
orig_dav

# This calculation matches the original to the decimal place


# Replication analyses - z-test with appropriate effect size --------

rep_test <- compare_smd(
  smd1 = abs(orig_dz$d),
  n1 = orig_desc$count[1],
  smd2 = rep_dz$d,
  n2 = rep_desc$count[1],
  paired = TRUE,
  alternative = "greater")
rep_test

# Replication analyses - z-test with reported effect size --------

rep_test <- compare_smd(
  smd1 = abs(orig_dav$d),
  n1 = orig_desc$count[1],
  smd2 = -rep_dav$d,
  n2 = rep_desc$count[1],
  paired = TRUE,
  alternative = "greater")
rep_test
