# Additional scripts for lung Cancer Signature Analysis with MANATEE

multigroupForwardSearch <- function(pos.genes, neg.genes, forwardThresh = 0, pooledDataObjectList, init.pos = NULL, init.neg = NULL, force.posneg = TRUE, numCores = 1)
{
  perf.meas = "AUC"
  if(is.null(init.pos)){
    pos = NULL
  }else{
    pos = init.pos
  }
  if(is.null(init.neg)){
    neg = NULL
  }else{
    neg = init.neg
  }
  if(force.posneg){
    if(!is.null(init.pos) || !is.null(init.neg)){
      warning("force.posneg cannot be used if init.pos or init.neg is set, so it will be ignored")
      force.posneg = FALSE
    }else if(length(pos.genes) == 0 || length(neg.genes) == 0){
      warning("force.posneg cannot be used if pos.genes or neg.genes is empty, so it will be ignored")
      force.posneg = FALSE
    }
  }
  count = 0
  perf.max = 0
  continue = TRUE
  genelist = c(pos.genes, neg.genes)
  bestgenelist = list()
  max.genes = length(pos.genes) + length(neg.genes)
  
  if(any(!init.pos %in% pos.genes)){stop("There are some genes in init.pos that do not appear in pos.genes")}
  if(any(!init.neg %in% neg.genes)){stop("There are some genes in init.neg that do not appear in neg.genes")}
  
  while(continue && (count < max.genes)){
    auclist = mclapply(mc.cores = numCores, genelist, function(gene){
      currpos = pos
      currneg = neg
      if(gene %in% pos.genes){
        currpos = append(currpos, gene)
      }else{
        currneg = append(currneg, gene)
      }
      
      aucVec = rep(0, length(pooledDataObjectList))
      for(i in 1:length(pooledDataObjectList)){
        if(!any(c(currpos,currneg) %in% rownames(pooledDataObjectList[[i]]$genes))){
          scores = rep(0, ncol(pooledDataObjectList[[i]]$genes))
        }else{
          scores = getGeneScores(pooledDataObjectList[[i]]$genes, currpos,currneg,out.missing = F)
        }
        aucVec[i] = as.numeric(calculateROC(as.numeric(as.character(pooledDataObjectList[[i]]$class)), as.numeric(scores))$auc)
      }
      auc = mean(aucVec, na.rm = T)
      if(force.posneg && count == 1){
        if(is.null(pos) && gene %in% neg.genes){auc = 0}
        if(is.null(neg) && gene %in% pos.genes){auc = 0}
      }
      return(auc)

    })
    # check for NULL values
    if(any(sapply(auclist,function(x) is.null(x)))){
      warning("Some values of auclist are NULL - this is due to an mclapply error. Reduce numCores and try again.")
    }
    auclist = unlist(auclist)
    curr.max = max(auclist, na.rm = T)
    perf.diff = curr.max - perf.max
    cat(sprintf("next best: %s\n", perf.diff))
    count = count + 1
    if(perf.diff > forwardThresh){
      perf.max = curr.max
      topindex = which.max(auclist) 
      if(length(topindex) > 1){topindex = topindex[1]} 
      topgene = genelist[topindex]
      if(topgene %in% pos.genes){
        cat(sprintf("Adding %s (up)\n", topgene))
        pos = append(pos, topgene)
      }else{
        cat(sprintf("Adding %s (down)\n", topgene))
        neg = append(neg, topgene)
      }
      genelist = genelist[-topindex]
      
      bestauc.vec = bestauchi.vec = bestauclo.vec = rep(0, length(pooledDataObjectList))
      for(i in 1:length(pooledDataObjectList)){
        best.scores = getGeneScores(pooledDataObjectList[[i]]$genes, pos, neg,out.missing = F)
        auc.all = calculateROC(as.numeric(as.character(pooledDataObjectList[[i]]$class)), as.numeric(best.scores))
        bestauc.vec[i] = as.numeric(auc.all$auc)
        bestauchi.vec[i] = as.numeric(auc.all$auc.CI[1])
        bestauclo.vec[i] = as.numeric(auc.all$auc.CI[2])
      }
      best.perf = list(auc = mean(bestauc.vec), auc.CI = c(mean(bestauchi.vec), mean(bestauclo.vec)))
      
      best.auc = as.numeric(best.perf$auc)
      
      blurb = sprintf("%s Genes, %s =%s (95%% CI %s-%s): %s up and %s down", count, perf.meas, round(best.auc,3), round(best.perf$auc.CI[1], 3),
                      round(best.perf$auc.CI[2],3), paste(pos,collapse = ", "), paste(neg, collapse = ", "))
      bestgenelist[[as.character(count)]] = list(numGenes = count, upgenes = pos, downgenes = neg, perf.meas = perf.meas, perf = best.auc,
                                                 perf.ci.lower = best.perf$auc.CI[1], perf.ci.upper = best.perf$auc.CI[2], blurb = blurb)
    }else{
      cat(sprintf("\nFinal %s =%s\n%s upgenes and %s downgenes chosen\n",perf.meas,curr.max, length(pos), length(neg)))
      continue = FALSE
    }
    if(continue && count >= max.genes){
      cat(sprintf("\nFinal %s =%s\n%s upgenes and %s downgenes chosen\n", perf.meas,curr.max, length(pos), length(neg)))
    }
  }
  
  bestgeneDT = data.table(do.call(rbind, lapply(bestgenelist, function(x){
    data.frame(numGenes = x$numGenes, perf = x$perf, perf.ci.lower = x$perf.ci.lower, perf.ci.upper = x$perf.ci.upper,
               upgenes = paste(x$upgenes,collapse = " / "), downgenes = paste(x$downgenes,collapse = " / "), blurb=x$blurb,
               stringsAsFactors = FALSE)
  })))
  names(bestgeneDT)[2] = perf.meas
  names(bestgeneDT)[3] = paste0(perf.meas, ".ci.lower")
  names(bestgeneDT)[4] = paste0(perf.meas, ".ci.upper")
  
  return(list(upgenes = pos, downgenes = neg, bestgenelist = bestgenelist, bestgeneDT = bestgeneDT))
}

replace_names <- function(labelVec, currNames, replaceName)
{
  currNames = gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", currNames) #escape all special characters
  labelVec[grep(paste0(currNames,collapse ="|"), labelVec)] = gsub(paste0(currNames,collapse ="|"), replaceName, labelVec[grep(paste0(currNames,collapse ="|"), labelVec)])
  return(labelVec)
}

LC_create_labels <- function(GSElist)
{
  #start making group labels
  GSElist$pheno$group10 = GSElist$pheno$group
  
  
  #first level is collapsing really similar stuff and removing extraneous info
  GSElist$pheno$group9 = GSElist$pheno$group10
  
  benign.pulmonary.tumor.names = c("GGO Nodules", "Pulmonary Hamartoma", "Benign Pulmonary Nodule", "Benign Lung Tumor")
  post.vaccination.names = c("Advanced RCC (Post Vaccination)")
  unknown.COPD.names = c("Lung Adenocarcinoma/Unknown COPD", "Lung B-Adenocarcinoma/Unknown COPD", "Benign Pulmonary Tumor/Unknown COPD", "Lung Cirrhosis/Unknown COPD",
                         "Lung Psoriasis/Unknown COPD", "Lung SCC/Unknown COPD", "Granulomatous Lung Inflammation/Unknown COPD", "PAH/Unknown COPD")
  treatment.time.names = c("ALL/cGVHD Post Treatment (1 Month)", "ALL/cGVHD Post Treatment (2 Months)", "Aplastic Anemia/cGVHD Post Treatment (2 Months)", "CLL/cGVHD Post Treatment (1 Month)",
                           "CLL/cGVHD Post Treatment (2 Months)", "CML/cGVHD Post Treatment (1 Month)", "HIV Post Treatment (1 Week)", "HIV Post Treatment (12 Weeks)", "HIV Post Treatment (2 Weeks)",
                           "HIV Post Treatment (4 Weeks)", "HIV Post Treatment (24 Weeks)", "MDS/cGVHD Post Treatment (1 Month)", "MDS/cGVHD Post Treatment (2 Months)",
                           "Myelofibrosis/cGVHD Post Treatment (1 Month)", "Myelofibrosis/cGVHD Post Treatment (2 Months)", "Non-Hodgkins Lymphoma/cGVHD Post Treatment (1 Month)",
                           "Non-Hodgkins Lymphoma/cGVHD Post Treatment (2 Months)", "Post Kidney Transplant (1 Week)", "Post Kidney Transplant (12 Weeks)", "Post Kidney Transplant (2 Weeks)",
                           "Post Kidney Transplant (4 Weeks)", "Post Kidney Transplant (8 Weeks)", "Psoriasis Post Treatment (2 Weeks)", "AML/cGVHD Post Treatment (1 Month)",
                           "AML/cGVHD Post Treatment (2 Months)", "Non-Hodgkin's Lymphoma/cGVHD Post Treatment (2 Months)", "Non-Hodgkin's Lymphoma/cGVHD Post Treatment (1 Month)")
  remove.smoker.names = c("Bipolar Disorder/Smoker", "Chronic Illness/Smoker", "Schizophrenia/Smoker", "MDD/Smoker", "MDD/Former Smoker")
  critical.illness.names = c("Critical Illness w/ Infection", "Critical Illness w/ Possible Infection", "Critical Illness w/ Probable Infection", "Critical Illness w/o Infection")
  invasive.fungal.infection.names = c("Invasive Aspergillosis")
  healthy.smoker.names = c("Healthy (Moist Snuff User)", "Healthy (Smoker/Motor Complications)", "Healthy (Smoker/Acute Environmental Exposure)", "Healthy (Smoker/Environmental Exposure)",
                           "Healthy (Smoker/Minor Illness)", "High Benzene Exposure (Smoker)", "Low Benzene Exposure (Smoker)", "Very High Benzene Exposure (Smoker)")
  healthy.former.smoker.names = c("Healthy (Former Smoker/Acute Environmental Exposure)", "Healthy (Former Smoker/Environmental Exposure)", "Very Low Benzene Exposure (Smoker)")
  healthyish.names = c("Healthy (Motor Complications)", "Healthy (Post Colonoscopy)", "Healthy (Snorer)", "Chronic Granulomatous Disease Carrier", "DYT1 Dystonia Carrier",
                       "Genetic Parkinson's Disease Risk", "HTLV-1 Carrier", "Pre-Diabetic (IFG)", "Healthy (Acute Environmental Exposure)", "Healthy (Breast-Feeding)",
                       "Healthy (Minor Illness)", "Healthy (Suspicious Mammogram)", "Healthy (Weight Loss)", "Preclinical RA", "Convalescent (Sepsis/Melioidosis)")
  healthy.names = c("Healthy (Pooled)", "Very Low Benzene Exposure")
  remove.healthy.names = c("Healthy (Bipolar Disorder/Mild Pain Disorder)", "Healthy (MDD/Mild Pain Disorder)", "Healthy (Mild Pain Disorder)", "Healthy (Hypertension)", "Healthy (PTSD)",
                           "Healthy (Environmental Exposure)")
  environmental.exposure.names = c("Environmental Exposure (Nickel)", "High Radiation Exposure", "Medium Radiation Exposure", "Low Radiation Exposure", "High Benzene Exposure",
                                   "Low Benzene Exposure", "Very High Benzene Exposure")
  vte.names = c("High Risk VTE", "Moderate Risk VTE", "Low Risk VTE", "VTE (Single)", "VTE (Recurrent)")
  ibd.names = c("IBD-U")
  sleep.apnea.names = c("Moderate Sleep Apnea", "Severe Sleep Apnea", "Severe Sleep Apnea (CPAP Treated)")
  myeloproliferative.neoplasm.names = c("Myeloproliferative Disorder")
  esrd.names = c("Pre Kidney Transplant", "ESRD/Uremia")
  myelofibrosis.names = c("PMF")
  oa.names = c("OA (Knee)")
  sah.names = c("SAH (Aneurysmal)")
  bipolar.disorder.names = c("Bipolar Disorder (Manic)", "Bipolar Disorder (Euthymic)")
  benign.breast.tumor.names = c("Breast Hamartoma", "Breast Cyst", "Fibroadenoma")
  t1d.names = c("Fulminant T1D")
  non.lung.cancer = c("Post Kidney Transplant Cancer (Non-Lung)")
  
  GSElist$pheno$group9 = replace_names(GSElist$pheno$group9,benign.pulmonary.tumor.names, "Benign Pulmonary Tumor")
  GSElist$pheno$group9[GSElist$pheno$group9 %in% post.vaccination.names] = gsub(" \\(Post Vaccination\\)", "",GSElist$pheno$group9[GSElist$pheno$group9 %in% post.vaccination.names])
  GSElist$pheno$group9[GSElist$pheno$group9 %in% unknown.COPD.names] = gsub("/Unknown COPD", "",GSElist$pheno$group9[GSElist$pheno$group9 %in% unknown.COPD.names])
  GSElist$pheno$group9[GSElist$pheno$group9 %in% treatment.time.names] = gsub(" \\([0-9]+ [Month,Week,Year].*\\)", "",GSElist$pheno$group9[GSElist$pheno$group9 %in% treatment.time.names])
  GSElist$pheno$group9[GSElist$pheno$group9 %in% remove.smoker.names] = gsub("/Smoker|/Former Smoker", "",GSElist$pheno$group9[GSElist$pheno$group9 %in% remove.smoker.names])
  GSElist$pheno$group9[GSElist$pheno$group9 %in% critical.illness.names] = "Critical Illness"
  GSElist$pheno$group9 = replace_names(GSElist$pheno$group9,invasive.fungal.infection.names, "Invasive Fungal Infection")
  GSElist$pheno$group9[GSElist$pheno$group9 %in% healthy.smoker.names] = "Healthy (Smoker)"
  GSElist$pheno$group9[GSElist$pheno$group9 %in% healthy.former.smoker.names] = "Healthy (Former Smoker)"
  GSElist$pheno$group9[GSElist$pheno$group9 %in% healthyish.names] = "Healthy-ish"
  GSElist$pheno$group9[GSElist$pheno$group9 %in% healthy.names] = "Healthy"
  GSElist$pheno$group9[GSElist$pheno$group9 %in% remove.healthy.names] = gsub("Healthy \\(|\\)$", "",GSElist$pheno$group9[GSElist$pheno$group9 %in% remove.healthy.names])
  GSElist$pheno$group9 = replace_names(GSElist$pheno$group9,environmental.exposure.names, "Environmental Exposure")
  GSElist$pheno$group9 = replace_names(GSElist$pheno$group9,vte.names, "VTE")
  GSElist$pheno$group9 = replace_names(GSElist$pheno$group9,ibd.names, "IBD")
  GSElist$pheno$group9 = replace_names(GSElist$pheno$group9,sleep.apnea.names, "Sleep Apnea")
  GSElist$pheno$group9 = replace_names(GSElist$pheno$group9,myeloproliferative.neoplasm.names, "Myeloproliferative Neoplasm")
  GSElist$pheno$group9 = replace_names(GSElist$pheno$group9,esrd.names, "ESRD")
  GSElist$pheno$group9 = replace_names(GSElist$pheno$group9,myelofibrosis.names, "Myelofibrosis")
  GSElist$pheno$group9 = replace_names(GSElist$pheno$group9,oa.names, "OA")
  GSElist$pheno$group9 = replace_names(GSElist$pheno$group9,sah.names, "SAH")
  GSElist$pheno$group9 = replace_names(GSElist$pheno$group9,bipolar.disorder.names, "Bipolar Disorder")
  GSElist$pheno$group9 = replace_names(GSElist$pheno$group9,benign.breast.tumor.names, "Benign Breast Tumor")
  GSElist$pheno$group9 = replace_names(GSElist$pheno$group9, t1d.names, "T1D")
  GSElist$pheno$group9 = replace_names(GSElist$pheno$group9, non.lung.cancer, "Cancer (Non-Lung)")
  GSElist$pheno$group9 = gsub("^/|/$|//", "",GSElist$pheno$group9)
  GSElist$pheno$group9[GSElist$pheno$group9 == "Benign Breast Tumor/Benign Breast Tumor"] = "Benign Breast Tumor"
  
  
  #second stage is same as the first, but more thorough. also will collapse most severity info
  GSElist$pheno$group8 = GSElist$pheno$group9
  
  post.chemo.names = c("ALL (Post Chemo)", "AML (Post Chemo)", "Burkitt Lymphoma & Leukemia (Post Chemo)", "CLL (Post Chemo)", "CML (Post Chemo)", "Mantle Cell Lymphoma (Post Chemo)",
                       "SSc (Post Chemo)", "DLBCL (Post Chemo)", "Follicular Lymphoma (Post Chemo)", "MDS (Post Chemo)", "Precursor B-Cell Lymphoblastic Lymphoma & Leukemia (Post Chemo)",
                       "Aplastic Anemia (Post Chemo)", "Tonsil SCC (Post Chemo)", "T-Cell ALL (Post Chemo)", "Non-Hodgkin's Lymphoma (Post Chemo)", "Multiple Myeloma (Post Chemo)")
  post.treatment.names = c("ALL/cGVHD Post Treatment", "AML/cGVHD Post Treatment", "Aplastic Anemia/cGVHD Post Treatment", "CLL/cGVHD Post Treatment", "CML/cGVHD Post Treatment",
                           "HIV Post Treatment", "Myelofibrosis/cGVHD Post Treatment", "Non-Hodgkin's Lymphoma/cGVHD Post Treatment", "Psoriasis Post Treatment", "MDS/cGVHD Post Treatment")
  liver.transplant.cirrhosis.names = c("Liver Transplant/Laennec's Cirrhosis/Postnecrotic Cirrhosis Type C", "Liver Transplant/Laennec's Cirrhosis", "Liver Transplant/Cryptogenic Cirrhosis",
                                       "Liver Transplant/Primary Biliary Cirrhosis")
  liver.transplant.hepatitis.names = c("Liver Transplant/Autoimmune Hepatitis", "Liver Transplant/Hepatitis B", "Liver Transplant/Hepatitis C")
  rcc.names = c("Advanced RCC")
  pain.disorder.names = c("Mild Pain Disorder", "Severe Pain Disorder")
  stress.disorder.names = c("PTSD", "Childhood Trauma", "PTSD/Mild Pain Disorder", "PTSD/Pain Disorder", "PTSD/Severe Pain Disorder", "PTSD/Depression")
  ckd.names = c("CKD (Stage 2-3)", "CKD (Stage 1)", "CKD (Stage 2)", "CKD (Stage 3)", "CKD (Stage 4)")
  dystonia.names = c("Dopa-Responsive Dystonia", "DYT1 Dystonia", "DYT5 Dystonia")
  parkinsons.disease.names = c("Genetic Parkinson's Disease", "Atypical Parkinson's Disease", "Early Parkinson's Disease", "Parkinson's Disease Dementia")
  alzheimers.disease.names = c("Late Alzheimer's Disease")
  metastatic.melanoma.names = c("Metastatic Melanoma (Subclass M1a)", "Metastatic Melanoma (Subclass M1c)", "Metastatic Melanoma (Subclass M1c)/Pulmonary Lesion/Hepatic Lesion")
  asthma.names = c("Mild Asthma", "Mild Atopic Asthma", "Atopic Asthma", "Severe Asthma", "Severe Atopic Asthma", "Moderate Asthma")
  pregnant.names = c("MS/Pregnant (9 Months)", "Pregnant (3rd Trimester)", "Pregnant (9 Months)", "Pregnant (4 Months)", "Pregnant (6 Months)")
  vitiligo.names = c("Non-Segmental Vitiligo", "Segmental Vitiligo")
  breast.cancer.names = c("Post Surgery Breast Cancer", "Breast Cancer (IDC)", "Breast Cancer (IDC/DMS)", "Breast Cancer (IDC/ITC)", "Breast Cancer (ILC)", "Breast Cancer (ITC)",
                          "Breast Cancer (Medullary Carcinoma)", "Breast Cancer (Microinvasive DMS)", "Breast Cancer (Pure DMS)", "Breast Cancer (IDC/DCIS)",
                          "Breast Cancer (Microinvasive DCIS)", "Breast Cancer (Pure DCIS)")
  sle.names = c("pSLE", "SLE (LN-)", "SLE (LN+)")
  schizophrenia.names = c("Severe Schizophrenia", "Schizophrenia/Mild Pain Disorder", "Schizophrenia/Severe Pain Disorder")
  atll.names = c("ATLL (Acute)", "ATLL (Chronic)", "ATLL (Smoldering)")
  vasculitis.names = c("Granulomatosis with Polyangiitis")
  ms.names = c("RRMS", "PPMS", "SPMS", "CIS")
  cd.names = c("CD (Active)", "CD (Inactive)")
  uc.names = c("UC (Active)", "UC (Inactive)")
  cad.names = c("Intermediate CAD", "CAD/Diabetes", "PAD")
  chronic.myocardial.infarction.names = c("Myocardial Infarction History")
  cerebrovascular.disease.names = c("SAH", "Vascular Dementia")
  post.kidney.transplant.names = c("Post Kidney Transplant")
  hiv.names = c("HIV Treatment Failure")
  
  GSElist$pheno$group8[GSElist$pheno$group8 %in% post.chemo.names] = gsub(" \\(Post Chemo\\)", "",GSElist$pheno$group8[GSElist$pheno$group8 %in% post.chemo.names])
  GSElist$pheno$group8[GSElist$pheno$group8 %in% post.treatment.names] = gsub(" Post Treatment", "",GSElist$pheno$group8[GSElist$pheno$group8 %in% post.treatment.names])
  GSElist$pheno$group8[GSElist$pheno$group8 %in% liver.transplant.cirrhosis.names] = "Liver Transplant/Cirrhosis"
  GSElist$pheno$group8[GSElist$pheno$group8 %in% liver.transplant.hepatitis.names] = "Liver Transplant/Hepatitis"
  GSElist$pheno$group8 = replace_names(GSElist$pheno$group8,rcc.names, "RCC")
  GSElist$pheno$group8 = replace_names(GSElist$pheno$group8,pain.disorder.names, "Pain Disorder")
  GSElist$pheno$group8 = replace_names(GSElist$pheno$group8,stress.disorder.names, "Stress Disorder")
  GSElist$pheno$group8 = replace_names(GSElist$pheno$group8,ckd.names, "CKD")
  GSElist$pheno$group8 = replace_names(GSElist$pheno$group8,dystonia.names, "Dystonia")
  GSElist$pheno$group8 = replace_names(GSElist$pheno$group8,parkinsons.disease.names, "Parkinson's Disease")
  GSElist$pheno$group8 = replace_names(GSElist$pheno$group8,alzheimers.disease.names, "Alzheimer's Disease")
  GSElist$pheno$group8[GSElist$pheno$group8 %in% metastatic.melanoma.names] = gsub("Metastatic Melanoma \\(.*\\)", "Metastatic Melanoma",GSElist$pheno$group8[GSElist$pheno$group8 %in% metastatic.melanoma.names])
  GSElist$pheno$group8 = replace_names(GSElist$pheno$group8,asthma.names, "Asthma")
  GSElist$pheno$group8[GSElist$pheno$group8 %in% pregnant.names] = gsub("Pregnant \\(.*\\)", "Pregnant",GSElist$pheno$group8[GSElist$pheno$group8 %in% pregnant.names])
  GSElist$pheno$group8 = replace_names(GSElist$pheno$group8,vitiligo.names, "Vitiligo")
  GSElist$pheno$group8 = replace_names(GSElist$pheno$group8,breast.cancer.names, "Breast Cancer")
  GSElist$pheno$group8 = replace_names(GSElist$pheno$group8,sle.names, "SLE")
  GSElist$pheno$group8 = replace_names(GSElist$pheno$group8,schizophrenia.names, "Schizophrenia")
  GSElist$pheno$group8 = replace_names(GSElist$pheno$group8,atll.names, "ATLL")
  GSElist$pheno$group8 = replace_names(GSElist$pheno$group8,vasculitis.names, "Vasculitis")
  GSElist$pheno$group8 = replace_names(GSElist$pheno$group8,ms.names, "MS")
  GSElist$pheno$group8 = replace_names(GSElist$pheno$group8,cd.names, "CD")
  GSElist$pheno$group8 = replace_names(GSElist$pheno$group8,uc.names, "UC")
  GSElist$pheno$group8 = replace_names(GSElist$pheno$group8,cad.names, "CAD")
  GSElist$pheno$group8 = replace_names(GSElist$pheno$group8,chronic.myocardial.infarction.names, "Chronic Myocardial Infarction")
  GSElist$pheno$group8 = replace_names(GSElist$pheno$group8,cerebrovascular.disease.names, "Cerebrovascular Disease")
  GSElist$pheno$group8 = replace_names(GSElist$pheno$group8,post.kidney.transplant.names, "Kidney Transplant")
  GSElist$pheno$group8 = replace_names(GSElist$pheno$group8,hiv.names, "HIV")
  GSElist$pheno$group8 = gsub("^/|/$|//", "",GSElist$pheno$group8)
  
  
  #third stage is focused on collapsing similar conditions, but mostly for less important stuff
  GSElist$pheno$group7 = GSElist$pheno$group8
  
  acute.heart.failure.names = c("Acute Myocardial Infarction", "Acute Myocardial Infarction/T2D", "Acute Heart Failure/T2D", "Non-Ischemic Cardiomyopathy", "Ischemic Cardiomyopathy")
  chronic.heart.failure.names = c("Chronic Myocardial Infarction")
  cad.names = c("Atherosclerosis")
  neurodegenerative.disorder.names = c("ADLD", "CBD", "Dystonia", "MSA", "PSP", "MCI", "XDP", "ALS-like Disease")
  rheumatic.disease.names = c("Ankylosing Spondylitis", "Psoriatic Arthritis", "sJIA", "OA")
  vascular.disease.names = c("Behcet's Disease", "Vasculitis", "VTE", "Idiopathic Portal Hypertension")
  bipolar.disorder.names = c("Bipolar Disorder/Pain Disorder", "Bipolar Disorder/Depression", "Bipolar Disorder/Depressive or Anxiety Disorder")
  psychotic.disorder.names = c("Schizophrenia", "Schizoaffective Disorder", "Schizophrenia/Pain Disorder", "Schizoaffective Disorder/Pain Disorder", "Psychotic Disorder/Pain Disorder", "Psychosis",
                               "Psychosis/Pain Disorder")
  chronic.hepatitis.names = c("Chronic HCV", "Chronic Hepatitis B")
  skin.disease.names = c("Cutaneous Psoriasis", "Chronic Hives", "Hidradenitis Suppurativa", "Psoriasis", "Vitiligo")
  depression.anxiety.names = c("Depression", "MDD", "MDD/Pain Disorder", "Anxiety")
  poor.cvd.health.names = c("Hypercholesterolemia", "Hypertension", "Hypertension/Hypercholesterolemia", "Metabolic Syndrome")
  poor.cvd.health.cerebrovascular.names = c("Ischemic Stroke History/Diabetes/Hypercholesterolemia", "Ischemic Stroke History/Diabetes/Hypertension/Hypercholesterolemia",
                                            "Ischemic Stroke History/Hypercholesterolemia", "Ischemic Stroke History/Hypertension", "Ischemic Stroke History/Hypertension/Hypercholesterolemia")
  cerebrovascular.names = c("Ischemic Stroke History")
  renal.disease = c("CKD", "IgA Nephropathy", "Membranous Nephropathy")
  psychiatric.disorder.names = c("Mood Disorder", "Mood Disorder/Pain Disorder", "OCD")
  ibd.names = c("IBD", "IBS", "UC", "CD", "Diverticulitis")
  kidney.transplant.names = c("Kidney Transplant (Chronic Rejection)", "Kidney Transplant (Stable)", "Kidney Transplant (Tolerant)")
  liver.transplant.names = c("Liver Transplant (ACR)", "Liver Transplant (Non-Tolerant)", "Liver Transplant (Tolerant)", "Liver Transplant/Cirrhosis", "Liver Transplant/Fulminant Hepatic Failure",
                             "Liver Transplant/Hepatitis", "Liver Transplant/NASH", "Liver Transplant/PSC")
  metastatic.melanoma.2.names = c("Metastatic Melanoma/Pulmonary Lesion/Hepatic Lesion")
  sarcoidosis.names = c("Sarcoidosis/HCV", "Uncomplicated Sarcoidosis", "Sarcoidosis/Bronchiolitis Obliterans", "Sarcoidosis/Bronchiolitis obliterans", "Sarcoidosis/COPD")
  ssc.pah.names = c("SSc/PAH/ILD")
  diabetes.names = c("T1D", "T2D")
  chronic.illness.names = c("CFS", "XLA", "Pain Disorder")
  hsct.names = c("HSCT (Non-Tolerant)", "HSCT (Tolerant)")
  active.tb.names = c("Active TB/Diabetes", "Active TB/IRIS")
  latent.tb.names = c("Latent TB/Diabetes")
  respiratory.distress.names = c("Non-Cancer/Pneumonia Respiratory Distress")
  pdac.chronic.pancreatitis.names = c("PDAC/Chronic Pancreatitis/Diabetes")
  pdac.names = c("PDAC/Diabetes")
  
  GSElist$pheno$group7 = replace_names(GSElist$pheno$group7,acute.heart.failure.names, "Acute Heart Failure")
  GSElist$pheno$group7 = replace_names(GSElist$pheno$group7,chronic.heart.failure.names, "Chronic Heart Failure")
  GSElist$pheno$group7 = replace_names(GSElist$pheno$group7,cad.names, "CAD")
  GSElist$pheno$group7 = replace_names(GSElist$pheno$group7, neurodegenerative.disorder.names, "Neurodegenerative Disorder")
  GSElist$pheno$group7 = replace_names(GSElist$pheno$group7,rheumatic.disease.names, "Rheumatic Disease")
  GSElist$pheno$group7 = replace_names(GSElist$pheno$group7,vascular.disease.names, "Vascular Disease")
  GSElist$pheno$group7 = replace_names(GSElist$pheno$group7,bipolar.disorder.names, "Bipolar Disorder")
  GSElist$pheno$group7 = replace_names(GSElist$pheno$group7,psychotic.disorder.names, "Psychotic Disorder")
  GSElist$pheno$group7 = replace_names(GSElist$pheno$group7,chronic.hepatitis.names, "Chronic Hepatitis")
  GSElist$pheno$group7 = replace_names(GSElist$pheno$group7,skin.disease.names, "Skin Disease")
  GSElist$pheno$group7 = replace_names(GSElist$pheno$group7,depression.anxiety.names, "Depressive or Anxiety Disorder")
  GSElist$pheno$group7[GSElist$pheno$group7 %in% poor.cvd.health.names] = "Poor Cardiovascular Health"
  GSElist$pheno$group7[GSElist$pheno$group7 %in% poor.cvd.health.cerebrovascular.names] = "Poor Cardiovascular Health/Cerebrovascular Disease"
  GSElist$pheno$group7[GSElist$pheno$group7 %in% cerebrovascular.names] = "Cerebrovascular Disease"
  GSElist$pheno$group7 = replace_names(GSElist$pheno$group7,renal.disease, "Renal Disease")
  GSElist$pheno$group7 = replace_names(GSElist$pheno$group7,psychiatric.disorder.names, "Psychiatric Disorder")
  GSElist$pheno$group7 = replace_names(GSElist$pheno$group7,ibd.names, "IBD")
  GSElist$pheno$group7 = replace_names(GSElist$pheno$group7,kidney.transplant.names, "Kidney Transplant")
  GSElist$pheno$group7 = replace_names(GSElist$pheno$group7, liver.transplant.names, "Liver Transplant")
  GSElist$pheno$group7 = replace_names(GSElist$pheno$group7,metastatic.melanoma.2.names, "Metastatic Melanoma")
  GSElist$pheno$group7 = replace_names(GSElist$pheno$group7,sarcoidosis.names, "Sarcoidosis")
  GSElist$pheno$group7[GSElist$pheno$group7 %in% ssc.pah.names] = "SSc/PAH"
  GSElist$pheno$group7 = replace_names(GSElist$pheno$group7,diabetes.names, "Diabetes")
  GSElist$pheno$group7 = replace_names(GSElist$pheno$group7,chronic.illness.names, "Chronic Illness")
  GSElist$pheno$group7 = replace_names(GSElist$pheno$group7,hsct.names, "HSCT")
  GSElist$pheno$group7 = replace_names(GSElist$pheno$group7,active.tb.names, "Active TB")
  GSElist$pheno$group7 = replace_names(GSElist$pheno$group7, latent.tb.names, "Latent TB")
  GSElist$pheno$group7 = replace_names(GSElist$pheno$group7,respiratory.distress.names, "Respiratory Distress")
  GSElist$pheno$group7 = replace_names(GSElist$pheno$group7,pdac.chronic.pancreatitis.names, "PDAC/Chronic Pancreatitis")
  GSElist$pheno$group7 = replace_names(GSElist$pheno$group7,pdac.names, "PDAC")
  GSElist$pheno$group7 = gsub("^/|/$|//", "",GSElist$pheno$group7)
  GSElist$pheno$group7[GSElist$pheno$group7 == "Depressive or Anxiety Disorder/Depressive or Anxiety Disorder"] = "Depressive or Anxiety Disorder"
  GSElist$pheno$group7[GSElist$pheno$group7 == "Lung Skin Disease"] = "Lung Psoriasis"
  
  #fourth stage is collapsing infecting pathogen info to just bacterial/viral
  GSElist$pheno$group6 = GSElist$pheno$group7
  
  asymptomatic.lrti.bacterial = c()
  asymptomatic.lrti.viral = c("Asymptomatic LRTI (RSV)", "Asymptomatic LRTI (Influenza)", "Asymptomatic LRTI (Rhinovirus)", "Asymptomatic LRTI (RSV)")
  asymptomatic.lrti.mixed = c()
  asymptomatic.lrti.other = c()
  lrti.bacterial = c("LRTI (Bacterial)", "LRTI (MRSA)", "LRTI (K. Pneumoniae)",  "LRTI (S. Aureus)", "LRTI (M. Catarrhalis)", "LRTI (MSSA)",
                     "LRTI (H. Influenzae)", "LRTI (S. Pneumoniae)", "LRTI (S. Marcescens)", "LRTI (S. Pyogenes)", "LRTI (P. Aeruginosa)")
  lrti.viral = c("LRTI (Influenza/Rhinovirus)", "LRTI (Coronavirus)", "LRTI (Rhinovirus)", "LRTI (Viral)", "LRTI (Influenza)", "LRTI (Enterovirus)",
                 "LRTI (RSV)", "LRTI (HMPV)", "ARDS (Influenza)", "LRTI (Enterovirus/Rhinovirus)", "LRTI (Parainfluenza)", "LRTI (Rhinovirus/RSV)",
                 "LRTI (Rhinovirus/Enterovirus)", "LRTI (Influenza/Coronavirus)", "LRTI (Influenza/RSV)", "LRTI (Rhinovirus/Coronavirus)", "Severe LRTI (Influenza)",
                 "Moderate LRTI (Influenza)")
  lrti.mixed = c("LRTI (Bacterial/Viral)", "LRTI (Bacterial/Influenza)")
  lrti.other = c("LRTI (Unknown Febrile Illness)")
  pneumonia.bacterial = c("Lung Abscess (MSSA)", "Pneumonia (Probable Bacterial)", "Pneumonia (M. Pneumoniae)", "Pneumonia/Empyema (S. Pneumoniae)",
                          "Bacteremia/Pneumonia (S. Pneumoniae)", "Pneumonia (S. Pyogenes)", "Bacteremia/Pneumonia (MRSA)", "Pneumonia (S. Pneumoniae)" ,
                          "Pneumonia", "Pneumonia (MRSA)", "Pneumonia (Bacterial)" , "Bacteremia/Osteomyelitis/Pneumonia (MRSA)")
  pneumonia.viral = c("Pneumonia (Influenza)", "Pneumonia (Rhinovirus)", "Pneumonia (Probable Viral)", "Pneumonia (RSV)", "Pneumonia/ARDS (Influenza)" ,
                      "Pneumonia (Parainfluenza)", "Pneumonia (Rhinovirus/hMPV)", "Pneumonia (MRSA)", "Pneumonia (Mixed Viral)", "Pneumonia (Coronavirus)")
  pneumonia.mixed = c("Pneumonia (Mixed Infection)", "Pneumonia (Influenza/Bacterial)")
  pneumonia.other = c("Pneumonia/COPD", "Pneumonia/Diabetes", "Sepsis/Pneumonia")
  urti.bacterial = c("Pharyngitis (S. Pyogenes)", "Retropharyngeal Abscess (E. Corrodens)")
  urti.viral = c("URTI (Influenza)", "Bronchiolitis (RSV)")
  urti.mixed = c()
  urti.other = c("Pharyngitis")
  sepsis.bacterial = c("Sepsis (Bacterial)", "Sepsis (CoNS)", "Sepsis (Corynebacterium)", "Sepsis (E. Coli)", "Sepsis (Enterococcus)", "Sepsis (Gram-Negative Bacteria)",
                       "Sepsis (Gram-Positive Bacteria)", "Sepsis (S. Aureus)", "Sepsis (Salmonella)", "Sepsis (Streptococcus)", "Sepsis/Melioidosis")
  sepsis.viral = c()
  sepsis.mixed = c()
  sepsis.other = c("Sepsis (C. Albicans)", "Sepsis/Bacteremia", "Sepsis/UTI")
  sepsis.ards.bacterial = c("Sepsis/ARDS (A. Hydrophila)", "Sepsis/ARDS (CoNS)", "Sepsis/ARDS (E. Faecium)", "Sepsis/ARDS (Enterococcus)", "Sepsis/ARDS (K. Pneumoniae)",
                            "Sepsis/ARDS (S. Aureus)", "Sepsis/ARDS (S. Pneumoniae)", "Sepsis/ARDS/Melioidosis")
  sepsis.ards.viral = c()
  sepsis.ards.mixed = c()
  sepsis.ards.other = c()
  sepsis.and.possible.respiratory.distress.names = c("Sepsis & Possible Respiratory Distress/Diabetes")
  
  GSElist$pheno$group6[GSElist$pheno$group6 %in% asymptomatic.lrti.bacterial] = "Asymptomatic LRTI (Bacterial)"
  GSElist$pheno$group6[GSElist$pheno$group6 %in% asymptomatic.lrti.viral] = "Asymptomatic LRTI (Viral)"
  GSElist$pheno$group6[GSElist$pheno$group6 %in% asymptomatic.lrti.mixed] = "Asymptomatic LRTI (Mixed)"
  GSElist$pheno$group6[GSElist$pheno$group6 %in% asymptomatic.lrti.other] = "Asymptomatic LRTI"
  GSElist$pheno$group6[GSElist$pheno$group6 %in% lrti.bacterial] = "LRTI (Bacterial)"
  GSElist$pheno$group6[GSElist$pheno$group6 %in% lrti.viral] = "LRTI (Viral)"
  GSElist$pheno$group6[GSElist$pheno$group6 %in% lrti.mixed] = "LRTI (Mixed)"
  GSElist$pheno$group6[GSElist$pheno$group6 %in% lrti.other] = "LRTI"
  GSElist$pheno$group6[GSElist$pheno$group6 %in% pneumonia.bacterial] = "Pneumonia (Bacterial)"
  GSElist$pheno$group6[GSElist$pheno$group6 %in% pneumonia.viral] = "Pneumonia (Viral)"
  GSElist$pheno$group6[GSElist$pheno$group6 %in% pneumonia.mixed] = "Pneumonia (Mixed)"
  GSElist$pheno$group6[GSElist$pheno$group6 %in% pneumonia.other] = "Pneumonia"
  GSElist$pheno$group6[GSElist$pheno$group6 %in% urti.bacterial] = "URTI (Bacterial)"
  GSElist$pheno$group6[GSElist$pheno$group6 %in% urti.viral] = "URTI (Viral)"
  GSElist$pheno$group6[GSElist$pheno$group6 %in% urti.mixed] = "URTI (Mixed)"
  GSElist$pheno$group6[GSElist$pheno$group6 %in% urti.other] = "URTI"
  GSElist$pheno$group6[GSElist$pheno$group6 %in% sepsis.bacterial] = "Sepsis (Bacterial)"
  GSElist$pheno$group6[GSElist$pheno$group6 %in% sepsis.viral] = "Sepsis (Viral)"
  GSElist$pheno$group6[GSElist$pheno$group6 %in% sepsis.mixed] = "Sepsis (Mixed)"
  GSElist$pheno$group6[GSElist$pheno$group6 %in% sepsis.other] = "Sepsis"
  GSElist$pheno$group6[GSElist$pheno$group6 %in% sepsis.ards.bacterial] = "Sepsis/ARDS (Bacterial)"
  GSElist$pheno$group6[GSElist$pheno$group6 %in% sepsis.ards.viral] = "Sepsis/ARDS (Viral)"
  GSElist$pheno$group6[GSElist$pheno$group6 %in% sepsis.ards.mixed] = "Sepsis/ARDS (Mixed)"
  GSElist$pheno$group6[GSElist$pheno$group6 %in% sepsis.ards.other] = "Sepsis/ARDS"
  GSElist$pheno$group6[GSElist$pheno$group6 %in% sepsis.and.possible.respiratory.distress.names] = "Sepsis & Possible Respiratory Distress"
  
  
  #fifth stage is collapsing blood cancer info
  GSElist$pheno$group5 = GSElist$pheno$group6
  
  leukemia.names = c("ALL", "AML", "CML", "CLL", "AML/MDS", "Leukemia/MDS", "T-Cell ALL", "AML Remission", "LGL Leukemia")
  lymphoma.names = c("NLPHL", "Non-Hodgkin's Lymphoma", "DLBCL", "Follicular Lymphoma", "Mantle Cell Lymphoma")
  blood.cancer.names = c("MDS", "ATLL", "Multiple Myeloma", "Multiple Myeloma/ONJ", "Blood Cancer/ONJ", "Myelofibrosis", "Myeloproliferative Neoplasm",
                         "Precursor B-Cell Lymphoblastic Lymphoma & Leukemia", "ET", "PV", "Burkitt Lymphoma & Leukemia", "Aplastic Anemia")
  #blood.bm.disease = c("Aplastic Anemia")
  
  GSElist$pheno$group5 = replace_names(GSElist$pheno$group5, leukemia.names, "Leukemia")
  GSElist$pheno$group5 = replace_names(GSElist$pheno$group5, lymphoma.names, "Lymphoma")
  GSElist$pheno$group5 = replace_names(GSElist$pheno$group5,blood.cancer.names, "Blood Cancer")
  #GSElist$pheno$group5 = replace_names(GSElist$pheno$group5,blood.bm.disease, "Blood or BM Disease")
  
  
  #sixth stage is collapsing most things that aren't super relevant
  GSElist$pheno$group4 = GSElist$pheno$group5
  
  renal.disease.2.names = c("ESRD")
  respiratory.infection.bacterial.names = c("Asymptomatic LRTI (Bacterial)", "LRTI (Bacterial)", "Pneumonia (Bacterial)")
  respiratory.infection.viral.names = c("Asymptomatic LRTI (Viral)", "LRTI (Viral)", "Pneumonia (Viral)")
  respiratory.infection.names = c("LRTI", "LRTI (Mixed)", "Pneumonia (Mixed)", "URTI", "Pneumonia")
  respiratory.illness.names = c("Hypersensitivity Pneumonitis", "Sleep Apnea", "Non-Healthy (Lung)", "Respiratory Distress", "Lung Cirrhosis", "Lung Psoriasis",
                                "Lung Cirrhosis/COPD", "CF", "Granulomatous Lung Inflammation", "Granulomatous Lung Inflammation/COPD", "Sepsis/ARDS (Bacterial)",
                                "Sepsis/ARDS (Viral)", "Sepsis/ARDS (Mixed)", "Sepsis/ARDS")
  psychiatric.mental.disorder.names = c("Bipolar Disorder", "Depressive or Anxiety Disorder", "Psychiatric Disorder", "Psychotic Disorder", "Stress Disorder")
  neurodegenerative.disorder.2.names = c("Alzheimer's Disease", "Huntington's Disease", "ALS")
  blood.cancer.2.names = c("Blood or BM Disease", "Leukemia", "Lymphoma")
  gi.disease.names = c("Celiac Disease", "IBD")
  copd.names = c("COPD or Bronchiectasis", "IECOPD")
  non.respiratory.infection.names = c("Melioidosis", "Melioidosis/Diabetes", "Non-Respiratory Infection/Diabetes")
  transplant.names = c("HSCT", "Kidney Transplant", "cGVHD", "Liver Transplant")
  neuroendocrine.cancer.names = c("Metastatic Carcinoid Tumor")
  pregnant.2.names = c("MS/Pregnant")
  pancreatic.cancer.names = c("Pancreatic Cancer/Diabetes", "PDAC", "PDAC/Chronic Pancreatitis", "PDAC/Chronic Pancreatitis/Diabetes", "PDAC/Diabetes")
  sepsis.names = c("SIRS", "Sepsis (Bacterial)", "Sepsis (Viral)", "Sepsis (Mixed)")
  other.names = c("Critical Illness", "Environmental Exposure")
  benign.pulmonary.tumor.2.names = c("Benign Pulmonary Tumor/COPD")
  
  GSElist$pheno$group4 = replace_names(GSElist$pheno$group4,renal.disease.2.names, "Renal Disease")
  GSElist$pheno$group4 = replace_names(GSElist$pheno$group4,respiratory.infection.bacterial.names, "Respiratory Infection (Bacterial)")
  GSElist$pheno$group4 = replace_names(GSElist$pheno$group4,respiratory.infection.viral.names, "Respiratory Infection (Viral)")
  GSElist$pheno$group4 = replace_names(GSElist$pheno$group4,respiratory.infection.names, "Respiratory Infection")
  GSElist$pheno$group4 = replace_names(GSElist$pheno$group4,respiratory.illness.names, "Respiratory Illness")
  GSElist$pheno$group4 = replace_names(GSElist$pheno$group4,psychiatric.mental.disorder.names, "Psychiatric or Mental Disorder")
  GSElist$pheno$group4 = replace_names(GSElist$pheno$group4, neurodegenerative.disorder.2.names, "Neurodegenerative Disorder")
  GSElist$pheno$group4 = replace_names(GSElist$pheno$group4,blood.cancer.2.names, "Blood Cancer")
  GSElist$pheno$group4 = replace_names(GSElist$pheno$group4, gi.disease.names, "GI Disease")
  GSElist$pheno$group4 = replace_names(GSElist$pheno$group4,copd.names, "COPD")
  GSElist$pheno$group4 = replace_names(GSElist$pheno$group4, non.respiratory.infection.names, "Non-Respiratory Infection")
  GSElist$pheno$group4 = replace_names(GSElist$pheno$group4, transplant.names, "Transplant")
  GSElist$pheno$group4 = replace_names(GSElist$pheno$group4, neuroendocrine.cancer.names, "Neuroendocrine Cancer")
  GSElist$pheno$group4 = replace_names(GSElist$pheno$group4,pregnant.2.names, "Pregnant")
  GSElist$pheno$group4 = replace_names(GSElist$pheno$group4,pancreatic.cancer.names, "Pancreatic Cancer")
  GSElist$pheno$group4 = replace_names(GSElist$pheno$group4,sepsis.names, "Sepsis")
  GSElist$pheno$group4 = replace_names(GSElist$pheno$group4,other.names, "Other")
  GSElist$pheno$group4 = replace_names(GSElist$pheno$group4,benign.pulmonary.tumor.2.names, "Benign Pulmonary Tumor")
  GSElist$pheno$group4[GSElist$pheno$group4 == "Sepsis & Possible Respiratory Illness"] = "Sepsis & Possible Respiratory Distress"
  
  #seventh stage is collapsing relevant groups
  GSElist$pheno$group3 = GSElist$pheno$group4
  
  cardiovascular.issues.names = c("Poor Cardiovascular Health", "Acute Heart Failure", "Chronic Heart Failure", "CAD", "Diabetes", "Poor Cardiovascular Health/Cerebrovascular Disease")
  blood.cancer.3.names = c("Blood Cancer/Invasive Fungal Infection", "Blood Cancer/Transplant")
  copd.2.names = c("Severe COPD", "Mild COPD", "Moderate COPD")
  mycobacterial.infection.names = c("Active TB", "Latent TB", "pNTM", "Active TB/HIV", "Latent TB/HIV")
  cns.disease = c("MS", "Neurodegenerative Disorder", "Parkinson's Disease", "Cerebrovascular Disease")
  transplant.2.names = c("Lung Transplant", "Transplant/GI Disease")
  chronic.illness.2.names = c("Chronic Hepatitis", "Sjogren's Syndrome", "SSc", "SSc/PAH")
  other.respiratory.issue.names = c("Respiratory Illness", "PAH", "Sepsis & Possible Respiratory Distress")
  non.respiratory.infection.2.names = c("HIV", "Sepsis or SIRS", "Sepsis")
  other.2.names = c("Benign Breast Tumor", "Pregnant")
  rheumatic.disease.2.names = c("RA")
  
  GSElist$pheno$group3 = replace_names(GSElist$pheno$group3,cardiovascular.issues.names, "Cardiovascular Issues")
  GSElist$pheno$group3[GSElist$pheno$group3 %in% blood.cancer.3.names] = "Blood Cancer"
  GSElist$pheno$group3 = replace_names(GSElist$pheno$group3,copd.2.names, "COPD")
  GSElist$pheno$group3 = replace_names(GSElist$pheno$group3,mycobacterial.infection.names, "Mycobacterial Infection")
  GSElist$pheno$group3 = replace_names(GSElist$pheno$group3,cns.disease, "CNS Disease")
  GSElist$pheno$group3 = replace_names(GSElist$pheno$group3, transplant.2.names, "Transplant")
  GSElist$pheno$group3 = replace_names(GSElist$pheno$group3,chronic.illness.2.names, "Chronic Illness")
  GSElist$pheno$group3 = replace_names(GSElist$pheno$group3,other.respiratory.issue.names, "Other Respiratory Issue")
  GSElist$pheno$group3 = replace_names(GSElist$pheno$group3, non.respiratory.infection.2.names, "Non-Respiratory Infection")
  GSElist$pheno$group3 = replace_names(GSElist$pheno$group3,other.2.names, "Other")
  GSElist$pheno$group3 = replace_names(GSElist$pheno$group3,rheumatic.disease.2.names, "Rheumatic Disease")
  
  
  #eight stage is mostly for collapsing cancer
  GSElist$pheno$group2 = GSElist$pheno$group3
  
  non.lung.cancer.names = c("Blood Cancer", "Breast Cancer", "CRC", "Gastrointestinal Cancer", "HCC", "Metastatic Melanoma", "Pancreatic Cancer", "RCC", "Transplant/Neuroendocrine Cancer",
                            "Breast Cancer/Colon Cancer", "Breast Cancer/Blood Cancer", "Breast Cancer/Gastric Cancer", "Breast Cancer/Gastrointestinal Stromal Tumor",
                            "Breast Cancer/Ovarian Cancer", "Breast Cancer/Uterine Carcinoid Tumor", "Cancer (Non-Lung)", "Colon Cancer", "Colon Cancer/Gastric Cancer",
                            "Colon Cancer/Ovarian Cancer", "Colon Cancer/Prostate Cancer", "Epithelial Ovarian Cancer", "Gastric Cancer", "Ovarian Cancer", "Ovarian Cancer/Bladder Cancer",
                            "Tonsil SCC")
  lung.cancer.names = c("Lung Adenocarcinoma", "Lung SCC", "Lung LCC", "NSCLC", "SCLC", "Lung Adenocarcinoma/COPD", "Lung B-Adenocarcinoma", "Lung B-Adenocarcinoma/COPD", "Lung SCC/COPD",
                        "NSCLC/COPD", "Bronchial Carcinoma")
  
  GSElist$pheno$group2 = replace_names(GSElist$pheno$group2, non.lung.cancer.names, "Non-Lung Cancer")
  GSElist$pheno$group2 = replace_names(GSElist$pheno$group2, lung.cancer.names, "Lung Cancer")
  
  
  #ninth stage is for dividing into healthy, smoker, non-respiratory, respiratory, non-lung cancer, and lung cancer
  GSElist$pheno$group = GSElist$pheno$group2
  
  smoker.names = c("Healthy (Former Smoker)", "Healthy (Smoker)")
  respiratory.issue.names = c("Asthma", "Benign Pulmonary Tumor", "COPD", "IPF", "Mycobacterial Infection", "Other Respiratory Issue", "Respiratory Infection",
                              "Respiratory Infection (Bacterial)", "Respiratory Infection (Viral)", "Sarcoidosis")
  non.respiratory.issue.names = c("CNS Disease", "Healthy-ish", "Cardiovascular Issues", "Chronic Illness", "GI Disease", "Non-Respiratory Infection", "Other",
                                  "Psychiatric or Mental Disorder", "Renal Disease", "Rheumatic Disease", "Skin Disease", "SLE", "Transplant", "Vascular Disease")
  
  
  GSElist$pheno$group = replace_names(GSElist$pheno$group,smoker.names, "Smoker")
  GSElist$pheno$group = replace_names(GSElist$pheno$group,respiratory.issue.names, "Respiratory Issue")
  GSElist$pheno$group = replace_names(GSElist$pheno$group, non.respiratory.issue.names, "Non-Respiratory Issue")
  
  return(GSElist)
}






















