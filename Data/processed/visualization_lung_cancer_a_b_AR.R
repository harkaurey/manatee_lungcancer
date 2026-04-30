# Validation of the diagnostic lung cancer score in an independent prospectively enrolled cohort.
# Made by Alessandra Rossetti 30.04.26

# Preparation:
library(pastecs)
library(ggpubr); theme_set(theme_bw())
library(GGally)
library(performance)
library(DHARMa)
library(car)
library(ggeffects)
library(ggplot2)
library(knitr) 
library(dplyr)
library(ggeffects)
library(MuMIn)


dd<- read.csv("../../Desktop/Project/Data/processed/41698_2025_1043_MOESM3_ESM.csv",stringsAsFactors= T)
View(dd)

# Visualize the Boxplot:

dd$group<- factor(dd$group,levels= levels(dd$group)[c(2, 1, 3)])

ggplot(dd, aes(x = group, y = lungcancerscore)) +
  
  # Violinplot
  geom_violin(fill = "grey90",
              color = "grey60",
              alpha = 0.8,
              trim = FALSE) +
  
  # Points
  geom_point(aes(color = group),
             size = 2,
             alpha = 0.25,
             position = position_jitter(width = 0.15)) +
  
  # Boxplot
  geom_boxplot(width = 0.12,
               fill = "white",
               outlier.shape = NA) +
  
  # Achsenbeschriftung
  labs(x = "",
       y = "lung cancer score",
       title = "Lung cancer score, validation (UK)") +
  
  # Theme
  theme_minimal() +
  
  # Do not put any legend
  theme(legend.position = "none")
