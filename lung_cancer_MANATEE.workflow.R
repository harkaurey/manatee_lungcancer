# Lung Cancer Signature Analysis with MANATEE

#### Loading packages and functions ----
library(COCONUT)
library(ROCR)
library(MetaIntegrator)
library(ggplot2)
library(parallel)
library(OptimalCutpoints)
library(data.table)
library(xgboost)
library(caret)

# The following R script files can be found in https://github.com/Khatri-Lab/manatee_pnas/tree/main/MANATEE_functions
source("MANATEE_functions.R")
source("MANATEE_functions_new.R")
source("runManatee.R")
source("filterManatee.R")
source("bss_functions.R")
source("search_functions.R")

# scripts for lung cancer signature
source("lung_cancer_scripts.R")

#### Load, filter and format data ----
LC.GSEs = readRDS("/labs/khatrilab/hongzheng/lung/obj/LC.GSEs.rds")
LC.GSEs.names = readRDS("/labs/khatrilab/hongzheng/lung/obj/LC.GSEs.names.rds")
load("/labs/khatrilab/adityamr/lungcancer/All_datasets/gse12771gpl6097_Cosmos.RData")
load("/labs/khatrilab/adityamr/lungcancer/All_datasets/gse12771gpl6097_EPIC.RData")

# Diseases to be filtered out of the analysis
remove.names = c("Remove", "Sepsis Post Admission (1 Day)", "Sepsis Post Admission (2 Days)", "pNTM (Post-Treatment)", "TB Post-Treatment (28 weeks)",
                 "Metastatic Melanoma (Subclass M1b)", "TB Post-Treatment (8 weeks)", "TB Post-Treatment (6 months)", "TB Post-Treatment (3 months)",
                 "Unknown Illness", "TB Post-Treatment (2 weeks)", "TB Post-Treatment (2 months)", "TB Post-Treatment (6 months)", "ARDS (Non-infectious)",
                 "Malaria", "Abscess (MSSA)", "TB Post Treatment (1 year)", "Colon Cancer/NSCLC", "CF")
# Add both Post and Post- versions)
remove.names = unique(c(remove.names, gsub("Post-", "Post ", grep("Post-", remove.names, value= T)), gsub("Post ", "Post-", grep("Post ", remove.names, value= T))))

# Set the healthy names (i.e. which samples should be treated as baseline samples for COCONUT) 
healthy.names = c("Latent TB", "Healthy (Latent TB)", "Healthy (Smoker)", "Healthy (Snorer)", "Healthy (Suspicious Mammogram)", "Healthy (Breast-Feeding)", "Healthy (Moist Snuff User)",
                  "Healthy (Motor Complications)", "Healthy (Former Smoker)", "Healthy (Minor Illness)", "Healthy (Smoker/Minor Illness)", "Healthy (Environmental Exposure)",
                  "Healthy (Former Smoker/Environmental Exposure)", "Healthy (Weight Loss)", "Healthy (Hypertension)", "Healthy (Mild Pain Disorder)", "Healthy (Bipolar Disorder/Mild Pain Disorder)",
                  "Healthy (MDD/Mild Pain Disorder)", "Healthy (PTSD)", "Healthy (Pooled)")
healthy.names.COCONUT = c("Healthy", paste0(healthy.names, " - Forced Healthy"))

# Set the dummy columns (i.e. which columns should be included even if they're not present in every dataset. 
my.dummy.cols = c("Age", "child_adult", "Sex", "Smoking.Status", "LC.Stage", "BC.Stage", "CRC.Stage", "PDAC.Stage", "RCC.Stage", "Melanoma.Stage",
                  "MM.Stage", "NLPHL.Stage", "COPD.Stage", "Ethnicity", "Survivor", "TST", "APACHE.at.admission", "APACHEII.Score", "FEV1.FVC",
                  "FEV1", "FEV1.L", "Eosinophils", "Glucose.mg.L")

LC.GSEs = lapply(LC.GSEs, function(gse){
  # first, convert to character
  gse$pheno = unfactor_df(gse$pheno)
  # get rid of extraneous samples
  remove = which(gse$pheno$group %in% remove.names)
  if(length(remove) > 0){
    gse = removeSamples(gse, remove, expr=("expr" %in% names(gse)))
  }
  # substitute mark pseudohealthy samples
  if(!"Healthy" %in% gse$pheno$group && any(healthy.names %in% gse$pheno$group)){
    gse$pheno$group[gse$pheno$group %in% healthy.names] = paste0(gse$pheno$group[gse$pheno$group %in% healthy.names], " - Forced Healthy")
    cat(gse$formattedName)
    cat(" Psuedohealthy\n")
  }
  # change healthy to smoker/former smoker if there are also at least 5 true healthies in the dataset
  if("Smoking.Status" %in% names(gse$pheno) && !all(is.na(gse$pheno$Smoking.Status))){
    if("Healthy" %in% gse$pheno$group && (sum(na.omit(gse$pheno$Smoking.Status == "No")) >= 5) ){
      gse$pheno$group[gse$pheno$group == "Healthy" & gse$pheno$Smoking.Status == "Quit"] = "Healthy (Former Smoker)"
      gse$pheno$group[gse$pheno$group == "Healthy" & gse$pheno$Smoking.Status == "Yes"] = "Healthy (Smoker)"
      cat(gse$formattedName)
      cat(" Smoking\n")
    }
  }
  # Standardize Post vs. Post- (for now just remove the dash)
  gse$pheno$group = gsub("Post-", "Post ", gse$pheno$group)
  #add dummy columns
  for(col in my.dummy.cols){
    if(!col %in% colnames(gse$pheno)){
      gse$pheno[[col]] = NA
    }
  }
  return(gse)
})

# formatManateeData
LC.GSEs = formatManateeData(LC.GSEs, AddDatasetName = TRUE, AddControl0Class = TRUE, ControlNames = healthy.names.COCONUT)

#### Split data ----

# Split data into four batches, each containing one lung cancer dataset and > 50 other datasets
LC.GSEs.1 = LC.GSEs[LC.GSEs.names$LC.GSEs.1]
LC.GSEs.2 = LC.GSEs[LC.GSEs.names$LC.GSEs.2]
LC.GSEs.3 = LC.GSEs[LC.GSEs.names$LC.GSEs.3]
LC.GSEs.4 = LC.GSEs[LC.GSEs.names$LC.GSEs.4]

#### Conormalize and further split data ----

# Firstly, randomly split each set of GSE data into two equally-sized groups. 
# Then, conormalize each group separately. Use COCOImputation with NAcutoff set to 0.1, which means that any gene with 
# more than 10% NAs across all the samples will be removed, and any gene with fewer than 10% NAs will be kept. 
# Also, set the split.min.healthy to 5, meaning that if a dataset has fewer than 5 samples, it will not be split 
# between the two groups and instead will be kept entirely within one group or the other)

LC.DV.1 = makeDiscoValid(LC.GSEs.1, use.pheno.class = TRUE, remove.samples = NULL, seed = 1337, conorm.type = "COCOImputation", Impute.NAcutoff = 0.1, healthy.names = healthy.names.COCONUT,dummy.cols = my.dummy.cols, split.prop.discovery = 0.5, split.min.healthy = 5)
LC.DV.2 = makeDiscoValid(LC.GSEs.2, use.pheno.class = TRUE, remove.samples = NULL, seed = 1337, conorm.type = "COCOImputation", Impute.NAcutoff = 0.1, healthy.names = healthy.names.COCONUT,dummy.cols = my.dummy.cols, split.prop.discovery = 0.5, split.min.healthy = 5)
LC.DV.3 = makeDiscoValid(LC.GSEs.3, use.pheno.class = TRUE, remove.samples = NULL, seed = 1337, conorm.type = "COCOImputation", Impute.NAcutoff = 0.1, healthy.names = healthy.names.COCONUT,dummy.cols = my.dummy.cols, split.prop.discovery = 0.5, split.min.healthy = 5)
LC.DV.4 = makeDiscoValid(LC.GSEs.4, use.pheno.class = TRUE, remove.samples = NULL, seed = 1337, conorm.type = "COCOImputation", Impute.NAcutoff = 0.1, healthy.names = healthy.names.COCONUT,dummy.cols = my.dummy.cols, split.prop.discovery = 0.5, split.min.healthy = 5)

# Organize post-COCONUT data
# Remove samples from children (since this analysis is focused on people at higher risk of lung cancer)
# Use the LC_create_labels function to annotate the samples at varying levels of organization

LC.DV.1$Discovery$pheno$group = gsub(" - Forced Healthy", "", LC.DV.1$Discovery$pheno$group)
LC.DV.1$Validation$pheno$group = gsub(" - Forced Healthy", "", LC.DV.1$Validation$pheno$group)
LC.DV.1$Discovery = removeSamples(LC.DV.1$Discovery,which(LC.DV.1$Discovery$pheno$child_adult == "Child"), expr = F)
LC.DV.1$Validation = removeSamples(LC.DV.1$Validation,which(LC.DV.1$Validation$pheno$child_adult == "Child"), expr = F)
LC.DV.1$Discovery$pheno$group[LC.DV.1$Discovery$pheno$group == "Healthy" & LC.DV.1$Discovery$pheno$Smoking.Status == "Quit" & !is.na(LC.DV.1$Discovery$pheno$Smoking.Status)] = "Healthy (Former Smoker)"
LC.DV.1$Discovery$pheno$group[LC.DV.1$Discovery$pheno$group == "Healthy" & LC.DV.1$Discovery$pheno$Smoking.Status == "Yes" & !is.na(LC.DV.1$Discovery$pheno$Smoking.Status)] = "Healthy (Smoker)"
LC.DV.1$Validation$pheno$group[LC.DV.1$Validation$pheno$group == "Healthy" & LC.DV.1$Validation$pheno$Smoking.Status == "Quit" & !is.na(LC.DV.1$Validation$pheno$Smoking.Status)] = "Healthy (Former Smoker)"
LC.DV.1$Validation$pheno$group[LC.DV.1$Validation$pheno$group == "Healthy" & LC.DV.1$Validation$pheno$Smoking.Status == "Yes" & !is.na(LC.DV.1$Validation$pheno$Smoking.Status)] = "Healthy (Smoker)"
LC.DV.1$Discovery = LC_create_labels(LC.DV.1$Discovery)
LC.DV.1$Validation = LC_create_labels(LC.DV.1$Validation)
LC.DV.1$Discovery$class = makeClassVector(LC.DV.1$Discovery$pheno$group, "Lung Cancer")
LC.DV.1$Validation$class = makeClassVector(LC.DV.1$Validation$pheno$group, "Lung Cancer")

LC.DV.2$Discovery$pheno$group = gsub(" - Forced Healthy", "", LC.DV.2$Discovery$pheno$group)
LC.DV.2$Validation$pheno$group = gsub(" - Forced Healthy", "", LC.DV.2$Validation$pheno$group)
LC.DV.2$Discovery = removeSamples(LC.DV.2$Discovery,which(LC.DV.2$Discovery$pheno$child_adult == "Child"), expr = F)
LC.DV.2$Validation = removeSamples(LC.DV.2$Validation,which(LC.DV.2$Validation$pheno$child_adult == "Child"), expr = F)
LC.DV.2$Discovery$pheno$group[LC.DV.2$Discovery$pheno$group == "Healthy" & LC.DV.2$Discovery$pheno$Smoking.Status == "Quit" & !is.na(LC.DV.2$Discovery$pheno$Smoking.Status)] = "Healthy (Former Smoker)"
LC.DV.2$Discovery$pheno$group[LC.DV.2$Discovery$pheno$group == "Healthy" & LC.DV.2$Discovery$pheno$Smoking.Status == "Yes" & !is.na(LC.DV.2$Discovery$pheno$Smoking.Status)] = "Healthy (Smoker)"
LC.DV.2$Validation$pheno$group[LC.DV.2$Validation$pheno$group == "Healthy" & LC.DV.2$Validation$pheno$Smoking.Status == "Quit" & !is.na(LC.DV.2$Validation$pheno$Smoking.Status)] = "Healthy (Former Smoker)"
LC.DV.2$Validation$pheno$group[LC.DV.2$Validation$pheno$group == "Healthy" & LC.DV.2$Validation$pheno$Smoking.Status == "Yes" & !is.na(LC.DV.2$Validation$pheno$Smoking.Status)] = "Healthy (Smoker)"
LC.DV.2$Discovery = LC_create_labels(LC.DV.2$Discovery)
LC.DV.2$Validation = LC_create_labels(LC.DV.2$Validation)
LC.DV.2$Discovery$class = makeClassVector(LC.DV.2$Discovery$pheno$group, "Lung Cancer")
LC.DV.2$Validation$class = makeClassVector(LC.DV.2$Validation$pheno$group, "Lung Cancer")

LC.DV.3$Discovery$pheno$group = gsub(" - Forced Healthy", "", LC.DV.3$Discovery$pheno$group)
LC.DV.3$Validation$pheno$group = gsub(" - Forced Healthy", "", LC.DV.3$Validation$pheno$group)
LC.DV.3$Discovery = removeSamples(LC.DV.3$Discovery,which(LC.DV.3$Discovery$pheno$child_adult == "Child"), expr = F)
LC.DV.3$Validation = removeSamples(LC.DV.3$Validation,which(LC.DV.3$Validation$pheno$child_adult == "Child"), expr = F)
LC.DV.3$Discovery$pheno$group[LC.DV.3$Discovery$pheno$group == "Healthy" & LC.DV.3$Discovery$pheno$Smoking.Status == "Quit" & !is.na(LC.DV.3$Discovery$pheno$Smoking.Status)] = "Healthy (Former Smoker)"
LC.DV.3$Discovery$pheno$group[LC.DV.3$Discovery$pheno$group == "Healthy" & LC.DV.3$Discovery$pheno$Smoking.Status == "Yes" & !is.na(LC.DV.3$Discovery$pheno$Smoking.Status)] = "Healthy (Smoker)"
LC.DV.3$Validation$pheno$group[LC.DV.3$Validation$pheno$group == "Healthy" & LC.DV.3$Validation$pheno$Smoking.Status == "Quit" & !is.na(LC.DV.3$Validation$pheno$Smoking.Status)] = "Healthy (Former Smoker)"
LC.DV.3$Validation$pheno$group[LC.DV.3$Validation$pheno$group == "Healthy" & LC.DV.3$Validation$pheno$Smoking.Status == "Yes" & !is.na(LC.DV.3$Validation$pheno$Smoking.Status)] = "Healthy (Smoker)"
LC.DV.3$Discovery = LC_create_labels(LC.DV.3$Discovery)
LC.DV.3$Validation = LC_create_labels(LC.DV.3$Validation)
LC.DV.3$Discovery$class = makeClassVector(LC.DV.3$Discovery$pheno$group, "Lung Cancer")
LC.DV.3$Validation$class = makeClassVector(LC.DV.3$Validation$pheno$group, "Lung Cancer")

LC.DV.4$Discovery$pheno$group = gsub(" - Forced Healthy", "", LC.DV.4$Discovery$pheno$group)
LC.DV.4$Validation$pheno$group = gsub(" - Forced Healthy", "", LC.DV.4$Validation$pheno$group)
LC.DV.4$Discovery = removeSamples(LC.DV.4$Discovery,which(LC.DV.4$Discovery$pheno$child_adult == "Child"), expr = F)
LC.DV.4$Validation = removeSamples(LC.DV.4$Validation,which(LC.DV.4$Validation$pheno$child_adult == "Child"), expr = F)
LC.DV.4$Discovery$pheno$group[LC.DV.4$Discovery$pheno$group == "Healthy" & LC.DV.4$Discovery$pheno$Smoking.Status == "Quit" & !is.na(LC.DV.4$Discovery$pheno$Smoking.Status)] = "Healthy (Former Smoker)"
LC.DV.4$Discovery$pheno$group[LC.DV.4$Discovery$pheno$group == "Healthy" & LC.DV.4$Discovery$pheno$Smoking.Status == "Yes" & !is.na(LC.DV.4$Discovery$pheno$Smoking.Status)] = "Healthy (Smoker)"
LC.DV.4$Validation$pheno$group[LC.DV.4$Validation$pheno$group == "Healthy" & LC.DV.4$Validation$pheno$Smoking.Status == "Quit" & !is.na(LC.DV.4$Validation$pheno$Smoking.Status)] = "Healthy (Former Smoker)"
LC.DV.4$Validation$pheno$group[LC.DV.4$Validation$pheno$group == "Healthy" & LC.DV.4$Validation$pheno$Smoking.Status == "Yes" & !is.na(LC.DV.4$Validation$pheno$Smoking.Status)] = "Healthy (Smoker)"
LC.DV.4$Discovery = LC_create_labels(LC.DV.4$Discovery)
LC.DV.4$Validation = LC_create_labels(LC.DV.4$Validation)
LC.DV.4$Discovery$class = makeClassVector(LC.DV.4$Discovery$pheno$group, "Lung Cancer")
LC.DV.4$Validation$class = makeClassVector(LC.DV.4$Validation$pheno$group, "Lung Cancer")

#### Run basic Manatee ----

LC.Basic.1 = runManatee(LC.DV.1$Discovery, manatee.type = "Basic", seed = 1337, runLeaveOneOutAnalysis = FALSE)
LC.Basic.2 = runManatee(LC.DV.1$Validation, manatee.type = "Basic", seed = 1337, runLeaveOneOutAnalysis = FALSE)
LC.Basic.3 = runManatee(LC.DV.2$Discovery, manatee.type = "Basic", seed = 1337, runLeaveOneOutAnalysis = FALSE)
LC.Basic.4 = runManatee(LC.DV.2$Validation, manatee.type = "Basic", seed = 1337, runLeaveOneOutAnalysis = FALSE)
LC.Basic.5 = runManatee(LC.DV.3$Discovery, manatee.type = "Basic", seed = 1337, runLeaveOneOutAnalysis = FALSE)
LC.Basic.6 = runManatee(LC.DV.3$Validation, manatee.type = "Basic", seed = 1337, runLeaveOneOutAnalysis = FALSE)
LC.Basic.7 = runManatee(LC.DV.4$Discovery, manatee.type = "Basic", seed = 1337, runLeaveOneOutAnalysis = FALSE)
LC.Basic.8 = runManatee(LC.DV.4$Validation, manatee.type = "Basic", seed = 1337, runLeaveOneOutAnalysis = FALSE)

# For each group, filter out the top 1000 genes based on their absolute effect
LC.Basic.1 = filterManatee(LC.Basic.1, NumTopGenesThresh = 1000, NumTopGenesMetric = "EffectSize", isLeaveOneOut = F)
LC.Basic.2 = filterManatee(LC.Basic.2, NumTopGenesThresh = 1000, NumTopGenesMetric = "EffectSize", isLeaveOneOut = F)
LC.Basic.3 = filterManatee(LC.Basic.3, NumTopGenesThresh = 1000, NumTopGenesMetric = "EffectSize", isLeaveOneOut = F)
LC.Basic.4 = filterManatee(LC.Basic.4, NumTopGenesThresh = 1000, NumTopGenesMetric = "EffectSize", isLeaveOneOut = F)
LC.Basic.5 = filterManatee(LC.Basic.5, NumTopGenesThresh = 1000, NumTopGenesMetric = "EffectSize", isLeaveOneOut = F)
LC.Basic.6 = filterManatee(LC.Basic.6, NumTopGenesThresh = 1000, NumTopGenesMetric = "EffectSize", isLeaveOneOut = F)
LC.Basic.7 = filterManatee(LC.Basic.7, NumTopGenesThresh = 1000, NumTopGenesMetric = "EffectSize", isLeaveOneOut = F)
LC.Basic.8 = filterManatee(LC.Basic.8, NumTopGenesThresh = 1000, NumTopGenesMetric = "EffectSize", isLeaveOneOut = F)

#### Select the top genes ----
upgenes1 = LC.Basic.1$filterResults$nTop.e.po1000_tgp_cf$upGeneNames
downgenes1 = LC.Basic.1$filterResults$nTop.e.po1000_tgp_cf$downGeneNames

upgenes2 = LC.Basic.2$filterResults$nTop.e.po1000_tgp_cf$upGeneNames
downgenes2 = LC.Basic.2$filterResults$nTop.e.po1000_tgp_cf$downGeneNames

upgenes3 = LC.Basic.3$filterResults$nTop.e.po1000_tgp_cf$upGeneNames
downgenes3 = LC.Basic.3$filterResults$nTop.e.po1000_tgp_cf$downGeneNames

upgenes4 = LC.Basic.4$filterResults$nTop.e.po1000_tgp_cf$upGeneNames
downgenes4 = LC.Basic.4$filterResults$nTop.e.po1000_tgp_cf$downGeneNames

upgenes5 = LC.Basic.5$filterResults$nTop.e.po1000_tgp_cf$upGeneNames
downgenes5 = LC.Basic.5$filterResults$nTop.e.po1000_tgp_cf$downGeneNames

upgenes6 = LC.Basic.6$filterResults$nTop.e.po1000_tgp_cf$upGeneNames
downgenes6 = LC.Basic.6$filterResults$nTop.e.po1000_tgp_cf$downGeneNames

upgenes7 = LC.Basic.7$filterResults$nTop.e.po1000_tgp_cf$upGeneNames
downgenes7 = LC.Basic.7$filterResults$nTop.e.po1000_tgp_cf$downGeneNames

upgenes8 = LC.Basic.8$filterResults$nTop.e.po1000_tgp_cf$upGeneNames
downgenes8 = LC.Basic.8$filterResults$nTop.e.po1000_tgp_cf$downGeneNames

all.upgenes = Reduce(union, list(upgenes1, upgenes2, upgenes3, upgenes4, upgenes5, upgenes6, upgenes7, upgenes8))
all.downgenes = Reduce(union, list(downgenes1, downgenes2, downgenes3, downgenes4, downgenes5, downgenes6, downgenes7, downgenes8))

# Make a note of how many times each gene shows up in the top 1000 genes across all 8 groups, 
# and also identify which genes showed up at least twice ("presence 2+") and which genes showed up at least three times ("presence 3+")
upgene.presence = sapply(all.upgenes, function(gene){
  sum(c(gene %in% upgenes1, gene %in% upgenes2, gene %in% upgenes3, gene %in% upgenes4, 
        gene %in% upgenes5, gene %in% upgenes6, gene %in% upgenes7, gene %in% upgenes8))
})
upgene.presence = sort(upgene.presence,decreasing = T)
upgene.presence_2plus = names(upgene.presence[upgene.presence >= 2])
upgene.presence_3plus = names(upgene.presence[upgene.presence >= 3])

downgene.presence = sapply(all.downgenes, function(gene){
  sum(c(gene %in% downgenes1, gene %in% downgenes2, gene %in% downgenes3, gene %in% downgenes4, 
        gene %in% downgenes5, gene %in% downgenes6, gene %in% downgenes7, gene %in% downgenes8))
})

downgene.presence = sort(downgene.presence,decreasing = T)
downgene.presence_2plus = names(downgene.presence[downgene.presence >= 2])
downgene.presence_3plus = names(downgene.presence[downgene.presence >= 3])

# Remove overlap between up and down genes
overlap.genes = intersect(all.upgenes, all.downgenes)
all.upgenes = all.upgenes[-which(all.upgenes %in% overlap.genes)]
all.downgenes = all.downgenes[-which(all.downgenes %in% overlap.genes)]

# Remove upgenes that were either a downgene or had less than 0.1 effect size as an upgene
# Remove downgenes that were either a upgenes or had less than 0.1 effect size as an downgenes
# Because gse42834 has so few lung cancer samples, we loosened the criteria for those two groups (3 & 4)
upgenes1.remove = c(rownames(LC.Basic.1$manateeResults$downgenes), rownames(LC.Basic.1$manateeResults$upgenes[LC.Basic.1$manateeResults$upgenes$effectSize < 0.1,]))
downgenes1.remove = c(rownames(LC.Basic.1$manateeResults$upgenes), rownames(LC.Basic.1$manateeResults$downgenes[LC.Basic.1$manateeResults$downgenes$effectSize > -0.1,]))

upgenes2.remove = c(rownames(LC.Basic.2$manateeResults$downgenes), rownames(LC.Basic.2$manateeResults$upgenes[LC.Basic.2$manateeResults$upgenes$effectSize < 0.1,]))
downgenes2.remove = c(rownames(LC.Basic.2$manateeResults$upgenes), rownames(LC.Basic.2$manateeResults$downgenes[LC.Basic.2$manateeResults$downgenes$effectSize > -0.1,]))

upgenes5.remove = c(rownames(LC.Basic.5$manateeResults$downgenes), rownames(LC.Basic.5$manateeResults$upgenes[LC.Basic.5$manateeResults$upgenes$effectSize < 0.1,]))
downgenes5.remove = c(rownames(LC.Basic.5$manateeResults$upgenes), rownames(LC.Basic.5$manateeResults$downgenes[LC.Basic.5$manateeResults$downgenes$effectSize > -0.1,]))

upgenes6.remove = c(rownames(LC.Basic.6$manateeResults$downgenes), rownames(LC.Basic.6$manateeResults$upgenes[LC.Basic.6$manateeResults$upgenes$effectSize < 0.1,]))
downgenes6.remove = c(rownames(LC.Basic.6$manateeResults$upgenes), rownames(LC.Basic.6$manateeResults$downgenes[LC.Basic.6$manateeResults$downgenes$effectSize > -0.1,]))

upgenes7.remove = c(rownames(LC.Basic.7$manateeResults$downgenes), rownames(LC.Basic.7$manateeResults$upgenes[LC.Basic.7$manateeResults$upgenes$effectSize < 0.1,]))
downgenes7.remove = c(rownames(LC.Basic.7$manateeResults$upgenes), rownames(LC.Basic.7$manateeResults$downgenes[LC.Basic.7$manateeResults$downgenes$effectSize > -0.1,]))

upgenes8.remove = c(rownames(LC.Basic.8$manateeResults$downgenes), rownames(LC.Basic.8$manateeResults$upgenes[LC.Basic.8$manateeResults$upgenes$effectSize < 0.1,]))
downgenes8.remove = c(rownames(LC.Basic.8$manateeResults$upgenes), rownames(LC.Basic.8$manateeResults$downgenes[LC.Basic.8$manateeResults$downgenes$effectSize > -0.1,]))

all.remove.upgenes = Reduce(union, list(upgenes1.remove, upgenes2.remove, upgenes5.remove, upgenes6.remove, upgenes7.remove, upgenes8.remove))
all.upgenes = all.upgenes[-which(all.upgenes %in% all.remove.upgenes)]

all.remove.downgenes = Reduce(union, list(downgenes1.remove, downgenes2.remove, downgenes5.remove, downgenes6.remove, 
                                          downgenes7.remove, downgenes8.remove))
all.downgenes = all.downgenes[-which(all.downgenes %in% all.remove.downgenes)]

#### Preparing data for multigroup forward search ----
g1 = unique(LC.Basic.1$pheno$group2)[unique(LC.Basic.1$pheno$group2) != "Lung Cancer"]
g2 = unique(LC.Basic.2$pheno$group2)[unique(LC.Basic.2$pheno$group2) != "Lung Cancer"]
g3 = unique(LC.Basic.3$pheno$group2)[unique(LC.Basic.3$pheno$group2) != "Lung Cancer"]
g4 = unique(LC.Basic.4$pheno$group2)[unique(LC.Basic.4$pheno$group2) != "Lung Cancer"]
g5 = unique(LC.Basic.5$pheno$group2)[unique(LC.Basic.5$pheno$group2) != "Lung Cancer"]
g6 = unique(LC.Basic.6$pheno$group2)[unique(LC.Basic.6$pheno$group2) != "Lung Cancer"]
g7 = unique(LC.Basic.7$pheno$group2)[unique(LC.Basic.7$pheno$group2) != "Lung Cancer"]
g8 = unique(LC.Basic.8$pheno$group2)[unique(LC.Basic.8$pheno$group2) != "Lung Cancer"]

l1 = lapply(g1, function(name){
  return(removeSamples(LC.Basic.1,which(!LC.Basic.1$pheno$group2 %in% c(name,"Lung Cancer")),expr = F))
})
l2 = lapply(g2,function(name){
  return(removeSamples(LC.Basic.2,which(!LC.Basic.2$pheno$group2 %in% c(name,"Lung Cancer")),expr = F))
})
l3 = lapply(g3,function(name){
  return(removeSamples(LC.Basic.3,which(!LC.Basic.3$pheno$group2 %in% c(name,"Lung Cancer")),expr = F))
})
l4 = lapply(g4,function(name){
  return(removeSamples(LC.Basic.4,which(!LC.Basic.4$pheno$group2 %in% c(name,"Lung Cancer")),expr = F))
})
l5 = lapply(g5,function(name){
  return(removeSamples(LC.Basic.5,which(!LC.Basic.5$pheno$group2 %in% c(name,"Lung Cancer")),expr = F))
})
l6 = lapply(g6,function(name){
  return(removeSamples(LC.Basic.6,which(!LC.Basic.6$pheno$group2 %in% c(name,"Lung Cancer")),expr = F))
})
l7 = lapply(g7,function(name){
  return(removeSamples(LC.Basic.7,which(!LC.Basic.7$pheno$group2 %in% c(name,"Lung Cancer")),expr = F))
})
l8 = lapply(g8,function(name){
  return(removeSamples(LC.Basic.8,which(!LC.Basic.8$pheno$group2 %in% c(name,"Lung Cancer")),expr = F))
})

# version with gse12771 (adjust to make sure it's weighted correctly - 25 copies of each)
new = list(gse12771gpl6102_Cosmos, gse12771gpl6102_EPIC, gse12771gpl6102_BC)
new.dup = c(rep(new, 25))
dataList2 = c(l1, l2, l3, l4, l5, l6, l7, l8, new.dup)

# add version with all 5 gse12771 datasets
new = list(gse12771gpl6102_Cosmos, gse12771gpl6102_EPIC, gse12771gpl6102_BC, gse12771gpl6097_Cosmos, gse12771gpl6097_EPIC)
new.dup = c(rep(new, 25))
dataList5 = c(l1, l2, l3, l4, l5, l6, l7, l8, new.dup)

#### Multigroup forward search ----
forward_multigroup_2_presence3plus = multigroupForwardSearch(pos.genes = upgene.presence_3plus, neg.genes = downgene.presence_3plus, pooledDataObjectList = dataList2, numCores = 12)
forward_multigroup_5_presence2plus = multigroupForwardSearch(pos.genes = upgene.presence_2plus, neg.genes = downgene.presence_2plus, pooledDataObjectList = dataList5, numCores = 12)

up22 = forward_multigroup_2_presence3plus$upgenes
down22 = forward_multigroup_2_presence3plus$downgenes

up25 = forward_multigroup_5_presence2plus$upgenes
down25 = forward_multigroup_5_presence2plus$downgenes

up_45 = sort(union(forward_multigroup_2_presence3plus$upgenes,
                   forward_multigroup_5_presence2plus$upgenes))
down_45 = sort(union(forward_multigroup_2_presence3plus$downgenes,
                     forward_multigroup_5_presence2plus$downgenes))
