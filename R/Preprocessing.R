library("ggplot2")
library("survival")
library("survminer")
library(rms)

rm(list=ls())

MDS_data = read.csv("MDS_data.csv")
n = dim(MDS_data)[1]

# Converting date to date format
MDS_data$dodx_prev = as.Date(MDS_data$dodx_prev)
MDS_data$dodx_hsct = as.Date(MDS_data$dodx_hsct)
MDS_data$dodx_ngs = as.Date(MDS_data$dodx_ngs)
MDS_data$docyto = as.Date(MDS_data$docyto)
MDS_data$dohsct = as.Date(MDS_data$dohsct)
MDS_data$dorel = as.Date(MDS_data$dorel)
MDS_data$doaml = as.Date(MDS_data$doaml)
MDS_data$dofup = as.Date(MDS_data$dofup)

# Converting genomic_group to character
MDS_data$genomic_group <- as.character(MDS_data$genomic_group)

# Converting IPSSR to integers
MDS_data$IPSSR <- as.character(MDS_data$IPSSR)
MDS_data$IPSSR[which(MDS_data$IPSSR == "very low")] = 1
MDS_data$IPSSR[which(MDS_data$IPSSR == "low")] = 2
MDS_data$IPSSR[which(MDS_data$IPSSR == "intermediate")] = 3
MDS_data$IPSSR[which(MDS_data$IPSSR == "high")] = 4
MDS_data$IPSSR[which(MDS_data$IPSSR == "very high")] = 5
MDS_data$IPSSR = as.numeric(MDS_data$IPSSR)

# dodx: first available diagnosis
MDS_data$dodx = pmin(MDS_data$dodx_ngs, MDS_data$dodx_hsct, MDS_data$dodx_prev, na.rm = TRUE)

# date_with_age: date of first diagnosis for which we have also the age
date_with_age = pmin(MDS_data$dodx_ngs[which(!is.na(MDS_data$age_ngs))], MDS_data$dodx_hsct[which(!is.na(MDS_data$age_hsct))],  
                     MDS_data$dodx_prev[which(!is.na(MDS_data$age_prev))], na.rm = TRUE)

# age_with_date: correspondent age                                  (ASSUMPTION: there are no ages that are not associated to any diagnosis        VERIFIED!)
age_with_date = pmin(MDS_data$age_ngs, MDS_data$age_hsct, MDS_data$age_prev, na.rm = TRUE)

# age: age at dodx
MDS_data$age = age_with_date - as.integer((date_with_age - MDS_data$dodx)/365)

# time: survival time in days from first diagnosis to last follow up
MDS_data$time = MDS_data$dofup - MDS_data$dodx

# wait_time: wait time in days from first diagnosis to the transplantation
MDS_data$wait_time = MDS_data$dohsct - MDS_data$dodx


# MDS_phase: initial MDS phase ("AML" if the first diagnosis is of AML)

# if patient is DNH
#   MDS_phase = WHO_ngs, namely = MDS_{5Q- / EB1 / ...}
# if patient is HSCT
#     if WHO_hsct != AML
#       MDS_phase = WHO_hsct, namely = MDS_{5Q- / EB1 / ...}
#     else
#       if patient has _prev (we use WHO_prev since we are looking at a WHO covariate)
#         MDS_phase = WHO_prev, namely = MDS_{5Q- / EB1 / ...}
#       else
#         MDS_phase = "AML"

MDS_data$MDS_phase = NA
MDS_data$MDS_phase[which(MDS_data$cohort == "DNH")] =
  MDS_data$WHO_ngs[which(MDS_data$cohort == "DNH")]
MDS_data$MDS_phase[which(MDS_data$cohort == "HSCT" & MDS_data$WHO_hsct != "AML-MRC")] =
  MDS_data$WHO_hsct[which(MDS_data$cohort == "HSCT" & MDS_data$WHO_hsct != "AML-MRC")]
MDS_data$MDS_phase[which(MDS_data$cohort == "HSCT" & MDS_data$WHO_hsct == "AML-MRC" & !is.na(MDS_data$WHO_prev))] =
  MDS_data$WHO_prev[which(MDS_data$cohort == "HSCT" & MDS_data$WHO_hsct == "AML-MRC" & !is.na(MDS_data$WHO_prev))]
MDS_data$MDS_phase[which(MDS_data$cohort == "HSCT" & MDS_data$WHO_hsct == "AML-MRC" & is.na(MDS_data$WHO_prev))] = "AML-MRC"

# AML / do_AML: dummy variable / date associated with AML evolution 

# if patient is DNH                                            (ASSUMPTION: there are no patients in the DNH cohort with AML as first diagnosis         VERIFIED!)
#   AML_phase = AML (the original one), namely = 1 or 0
#   do_AML = doaml (the original one)
# if patient is HSCT
#     if WHO_hsct != AML
#       AML_phase = 0
#       do_AML = NA
#     else
#       AML_phase = 1                                           (ASSUMPTION: there are no AML in WHO_prev         VERIFIED!)
#       do_AML = dodx_HSCT    

MDS_data$AML_phase = MDS_data$AML
MDS_data$do_aml = MDS_data$doaml
MDS_data$AML_phase[which(MDS_data$cohort == "HSCT" & MDS_data$WHO_hsct != "AML-MRC")] = 0
MDS_data$do_aml[which(MDS_data$cohort == "HSCT" & MDS_data$WHO_hsct != "AML-MRC")] = NA
MDS_data$AML_phase[which(MDS_data$cohort == "HSCT" & MDS_data$WHO_hsct == "AML-MRC")] = 1
MDS_data$do_aml[which(MDS_data$cohort == "HSCT" & MDS_data$WHO_hsct == "AML-MRC")] =
  MDS_data$dodx_hsct[which(MDS_data$cohort == "HSCT" & MDS_data$WHO_hsct == "AML-MRC")]

# age_AML: age at AML evolution
MDS_data$age_AML = MDS_data$age + as.integer((MDS_data$do_aml - MDS_data$dodx) / 365)

# MDS_data <- cbind(MDS_data, dodx, time, IPSSR)
# MDS_data <- data.frame(MDS_data, dodx, time, IPSSR)

write.csv(MDS_data ,"MDS_data_pre.csv", row.names = T)
