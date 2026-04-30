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

# Visualize the Boxplot of Figure 4A:

# Change the x-axis to contrl, benign and lung cancer.
dd$group<- factor(dd$group,levels= levels(dd$group)[c(2, 1, 3)])

# Calculate the p values between the groups:

t.test(lungcancerscore ~ group,
       data = subset(dd,
                     group %in% c("control", "benign")))

t.test(lungcancerscore ~ group,
       data = subset(dd,
                     group %in% c("benign", "lung cancer")))

t.test(lungcancerscore ~ group,
       data = subset(dd,
                     group %in% c("control", "lung cancer")))

# Find the p values in a different way, to see if the calculation was correct.

group1 <- dd$lungcancerscore[dd$group == "control"]
group2 <- dd$lungcancerscore[dd$group == "lung cancer"]

# Berechne den t-Test
test_result <- t.test(group1, group2)

# Zeige den p-Wert
print(test_result$p.value)



# Figure 4B:

ggplot(dd, aes(group, lungcancerscore)) +
  
  geom_violin(fill = "grey90") +
  
  geom_boxplot(width = 0.12,
               fill = "white",
               outlier.shape = NA) +
  
  geom_point(aes(color = group),
             alpha = 0.3,
             position = position_jitter(width = 0.15)) +
  labs(
    x = "",
    y = "Lung cancer score"
  ) +
  
  stat_compare_means(
    comparisons = list(
      c("control", "benign"),
      c("benign", "lung cancer"),
      c("control", "lung cancer")
    ),
    method = "wilcox.test",
    p.adjust.method = "bonferroni"
  ) +
  
  theme_minimal() +
  
  # Legende entfernen
  theme(legend.position = "none", axis.text.x = element_text(size = 20), axis.title.y = element_text(size = 20))



# Visualize Figure 4C:

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

