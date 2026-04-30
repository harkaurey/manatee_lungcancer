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

# Change the x-axis to contrl, benign and lung cancer.
dd$group<- factor(dd$group,levels= levels(dd$group)[c(2, 1, 3)])

# Calculate the p values between the groops:
t.test(lungcancerscore ~ group,
       data = subset(dd,
                     group %in% c("control", "benign")))


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
  
  stat_compare_means(
    comparisons = list(
      c("control", "benign"),
      c("benign", "lung cancer"),
      c("control", "lung cancer")
    ),
    method = "t.test", # Berechnet den t-Test
    
  ) +
  
  # Theme
  theme_minimal() +
  
  # Do not put any legend
  theme(legend.position = "none")

# Find the p values:

group1 <- dd$lungcancerscore[dd$group == "control"]
group2 <- dd$lungcancerscore[dd$group == "lung cancer"]

# Berechne den t-Test
test_result <- t.test(group1, group2)

# Zeige den p-Wert
print(test_result$p.value)


install.packages("pROC")

library(pROC)

# ROC 1
dd$lc_vs_all <- ifelse(dd$group == "lung cancer", 1, 0)

roc1 <- roc(dd$lc_vs_all,
            dd$lungcancerscore)

# ROC 2
dd2 <- subset(dd,
              group %in% c("benign", "lung cancer"))

dd2$lc_vs_benign <- ifelse(dd2$group == "lung cancer", 1, 0)

roc2 <- roc(dd2$lc_vs_benign,
            dd2$lungcancerscore)

# Plot
plot(roc1,
     col = "darkgrey",
     lwd = 3,
     legacy.axes = TRUE,
     main = "Lung cancer score, validation (UK)",
     xlab = "false positive rate (1-specificity)",
     ylab = "true positive rate (sensitivity)")

plot(roc2,
     add = TRUE,
     col = "lightgrey",
     lwd = 3)

legend("bottomright",
       legend = c(
         
         paste0(
           "lung cancer vs control/benign\n",
           "AUC = ", round(auc(roc1), 3),
           " (95% CI ",
           round(ci.auc(roc1)[1], 3),
           "-",
           round(ci.auc(roc1)[3], 3),
           ")"
         ),
         
         paste0(
           "lung cancer vs benign\n",
           "AUC = ", round(auc(roc2), 3),
           " (95% CI ",
           round(ci.auc(roc2)[1], 3),
           "-",
           round(ci.auc(roc2)[3], 3),
           ")"
         )
         
       ),
       col = c("darkgrey", "lightgrey"),
       lwd = 3,
       bty = "n")

