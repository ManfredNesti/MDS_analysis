rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("M3C")


library(car)
library(dplyr)
library(M3C)
library(e1071)
library(tree)
library(maptree)

load("tmp_MDS.RData")
load("new_MDS.RData")
load('mcshapiro.test.RData')
attach(new_MDS)

## ANOVA (2-WAY) TIME ~ CYTO_TRT + COHORT ------ ------

# dead <- which(fup_status==1)
# data_dead <- new_MDS[dead,]
# a <- which(data_dead$time==0)
# data_dead <- data_dead[-a,]
# detach(new_MDS)
# attach(data_dead)
# 
# ## make the groups
# cyto_dead = cyto_trt[which(fup_status==1)]
# cyto_dead = ifelse(cyto_dead=='yes', 1, 0)
# cohort_dead = cohort[which(fup_status==1)]
# cohort_dead = ifelse(cohort_dead=='HSCT', 1, 0)
# 
# group_00 = survtime_dead[which(cyto_dead==0 & cohort_dead==0)]
# group_10 = survtime_dead[which(cyto_dead==1 & cohort_dead==0)]
# group_01 = survtime_dead[which(cyto_dead==0 & cohort_dead==1)]
# group_11 = survtime_dead[which(cyto_dead==1 & cohort_dead==1)]
# 
# length(group_00)
# length(group_10)
# length(group_01)
# length(group_11)
# 
# ## set the time
# survtime_dead = as.numeric(time[which(fup_status==1)])
# 
# ## plot the data
# x11()
# par(mfrow=c(1,3)) #we see that the data is highly skewed
# plot(survtime_dead)
# hist(survtime_dead)
# boxplot(survtime_dead)
# 
# ## test normality of the groups
# shapiro.test(group_00)
# shapiro.test(group_10)
# shapiro.test(group_01)
# shapiro.test(group_11)
# 
# x11()
# par(mfrow=c(2,2))
# hist(group_00)
# hist(group_10)
# hist(group_01)
# hist(group_11)
# 
# ## Box-Cox transformations
# lambda = powerTransform ( time ) 
# time.transformed = bcPower ( as.numeric(time), lambda$lambda )
# boxplot(time.transformed)
# 
# ## repeat checking with transformed time
# survtime_dead = time.transformed
# 
# ## try removing outliers
# outliers <- boxplot(time)$out
# out_time_index <- which(time %in% outliers)
# data_dead <- data_dead[-out_time_index,]
# dim(data_dead)
# detach(data_dead)
# attach(data_dead)
# 
# boxplot(time.transformed)
# outliers.T <- boxplot(time.transformed)$out
# out_time_index.T <- which(time.transformed %in% outliers.T)
# data_dead.T <- data_dead[-out_time_index.T,]
# dim(data_dead.T)
# detach(data_dead)
# attach(data_dead)
# 
# 
# ## fit anova model
# fit <- aov(time.transformed ~ cyto_trt + cohort)
# 
# detach(data_dead)


## PREPROCESSING + PCA ------- -------
attach(new_MDS)
head(new_MDS)
new_MDS$agecat69 <- cut(new_MDS$age, breaks=c(0, 69, Inf), labels=c(1, 0))

nan_report <- function(df1, df2){
  cat("Initial Lenght: ", nrow(df1))
  cat("\nLenght after removing NaN Values: ", nrow(df2))
}

drops <- c("age_AML","MAC_reg", "donor_type", "HLA_mismatch", 
           "relapse", "wait_time", "age")
data_red <- new_MDS[ , !(names(new_MDS) %in% drops)]
head(data_red)

bmb_rm <- which(bmblasts <=30) #the records of bmblasts that are out of range (corrputed data)
data_red <- data_red[bmb_rm,]

nan_report(data_red, na.omit(data_red)) # 1956
## it is ok, we'll proceed with removing the NaNs

pca_data <- na.omit(data_red[ , !(names(data_red) %in% c("IPSSR"))])      # remove the patiets for which we have NA EXCLUDING IPSSR
pca_data$IPSSR = new_MDS$IPSSR[as.numeric(row.names(pca_data))] # recover the IPSSR for the new set
sum(is.na(pca_data$IPSSR))               # check that everything is fine
dim(pca_data)


pca_data$time <- as.numeric(pca_data$time)
detach(new_MDS)
attach(pca_data)
## We proceed with making one-hot encoded

gender_ = ifelse(gender=='F', 1, 0)
cohort_ = ifelse(cohort=='HSCT', 1, 0)
ggroup1 = ifelse(genomic_group=='1', 1, 0)
ggroup2 = ifelse(genomic_group=='2', 1, 0)
ggroup3 = ifelse(genomic_group=='3', 1, 0)
ggroup4 = ifelse(genomic_group=='4', 1, 0)
ggroup5 = ifelse(genomic_group=='5', 1, 0)
ggroup6 = ifelse(genomic_group=='6', 1, 0)
ggroup7 = ifelse(genomic_group=='7', 1, 0)
cyto = ifelse(cyto_trt == 'yes', 1, 0)



## Normalizing the numerical and categorical variables (FAMD algorithm)
# cat_vars <- c("IPSSR", "genomic_group", "KT_complex","cyto_trt", 
#               "fup_status", "AML_phase", "agecat69", "gender", "cohort")

# num_vars <- c(hgb, leu, neut, plt, bmblasts, n_mutgenes, n_mutations, time)
num_vars_scaled <- scale(cbind(hgb, leu, neut, plt, bmblasts, 
                               n_mutgenes, n_mutations, time))

# 
# scaled.center.num <- attr(num_vars_scaled, "scaled:center")
# scaled.scale.num <- attr(num_vars_scaled, "scaled:scale")
# 
# tmp_2<- (cbind(hgb, leu, neut, plt, bmblasts, 
#                n_mutgenes, n_mutations, time)[1:3,] - rep(scaled.center.num)/scaled.scale.num
# colMeans(tmp_2)
# colMeans(num_vars_scaled)

scaled.center.num <- attr(num_vars_scaled, "scaled:center")
scaled.scale.num <- attr(num_vars_scaled, "scaled:scale")
scaled.center.num
scaled.scale.num
# 
# as.matrix(rep(rep(0, dim(tmp_MDS)[2]), 
#               dim(tmp_MDS)[1]),0, dim(tmp_MDS)[1], dim(tmp_MDS)[2])

# scale_cat <- function(col){
#   scale.var.cat <- 1/(sqrt((sum(col)/length(col))))
#   return (c(col*scale.var.cat, scale.var.cat))
# }
scale_cat <- function(cols){
  scaled.var.cat <- sqrt((colSums(cols)/dim(cols)[1]))   #NOTE! here we don't have NA values, so it's correct to use dim(cols)[1] 
  return (list(t(t(cols)/scaled.var.cat), scaled.var.cat))
  
}

cat_vars_NOT_scaled <- cbind(gender_, cohort_, ggroup1, 
                             ggroup2, ggroup3, ggroup4,
                             ggroup5, ggroup6, ggroup7, 
                             cyto, KT_complex,
                             fup_status, AML_phase, agecat69)

scaled.scale.cat <- scale_cat(cat_vars_NOT_scaled)[2]
scaled.scale.cat <- as.numeric(scaled.scale.cat[[1]])
scaled.data.cat <- scale_cat(cat_vars_NOT_scaled)[1]
scaled.scale.cat
scaled.data.cat


# scale_cat(gender_)
# gender_scaled <- scale_cat(gender_)
# cohort_scaled <- scale_cat(cohort_)
# ggroup1_scaled <- scale_cat(ggroup1)
# ggroup2_scaled <- scale_cat(ggroup2)
# ggroup3_scaled <- scale_cat(ggroup3)
# ggroup4_scaled <- scale_cat(ggroup4)
# ggroup5_scaled <- scale_cat(ggroup5)
# ggroup6_scaled <- scale_cat(ggroup6)
# ggroup7_scaled <- scale_cat(ggroup7)
# cyto_scaled <- scale_cat(cyto)
# KT_complex_scaled <- scale_cat(KT_complex)
# fup_status_scaled <- scale_cat(fup_status)
# AML_phase_scaled <- scale_cat(AML_phase)
# agecat69_scaled <- scale_cat(as.numeric(agecat69))
# cat_vars_scaled <- cbind(gender_scaled, cohort_scaled, ggroup1_scaled, 
#                          ggroup2_scaled, ggroup3_scaled, ggroup4_scaled,
#                          ggroup5_scaled, ggroup6_scaled, ggroup7_scaled, 
#                          cyto_scaled, KT_complex_scaled,
#                          fup_status_scaled, AML_phase_scaled, agecat69_scaled)

## create a new dataset from scaled variables
pca_scaled_data <- data.frame(num_vars_scaled, scaled.data.cat)
head(pca_scaled_data)

# pca <- princomp(pca_data, scores = TRUE)
detach(pca_data)
attach(pca_scaled_data)
cov(pca_scaled_data)

pca <- princomp(pca_scaled_data, scores = T)
load.data <- pca$loadings

summary(pca) # 12 components needed for Explained Variance ~ 0.82

#plot the loadings
x11()
nk <- 12 # number of compoments to use
par(mar = c(1,1,1,1), mfrow = c(ncol(load.data[,c(1:nk)]),1))
for(i in 1:nk) barplot(load.data[,i], ylim = c(-1, 1))

#plot the scores
scores <- pca$scores
x11()
plot(scores[,1:2], pch=19, cex=0.5,col=c('blue','green')[factor(pca_scaled_data$cohort_)])
legend(x='topleft', legend=levels(factor(pca_scaled_data$cohort_)), fill=c('blue','green') )


# Dimensionality reduction
k = 22                        # number of columns to explain at least 80% of total variance
pca_data_red <- pca$scores[,1:k]  # reduced dataset
IPSSR_red = pca_data$IPSSR  # train, test and 63 missing values -> total length 2019
length(IPSSR_red)


## BUILDING THE CLASSIFIER WITH SVM ------- 

# Construction of datasets

train_test_idxs <- which(!is.na(IPSSR_red))    # indexes of the patients with known IPSSR
# train_test_set  <- pca$scores[train_test_idxs,] # Train set and test set together ->total length 1956
train_test_set  <- pca_data_red[train_test_idxs,] # Train set and test set together ->total length 1956

train_size = 1850    # size of the train set 
set.seed(1)
train_idxs = sample(dim(train_test_set)[1], train_size) # random indexes of the train set

train_set <- train_test_set[train_idxs,]   # train set -> 1800 elements
test_set <- train_test_set[-train_idxs,]   # test set  -> 156 elements
pred_set <- scores[which(is.na(IPSSR_red)),1:k]  # test on which we will predict IPSSR -> 63 elements

dim(train_set)
dim(test_set)
dim(pred_set)


# Construction of labels

train_test_IPSSR <- IPSSR_red[train_test_idxs]  # known values of IPSSR (so train and test together)
train_IPSSR <- as.factor(train_test_IPSSR[train_idxs])    # labels of the train set (1800 elements)
test_IPSSR <- as.factor(train_test_IPSSR[-train_idxs])   # labels of the test set  (156 elements)

length(train_IPSSR)
length(test_IPSSR)

table(train_IPSSR)      # Always check to have enough samples for each level
table(test_IPSSR)


# Perform classifier through svm


# Here we perform svm with parameters tuned using tune.svm command (gamma= 0.5, cost =10)

# gammas <- c(0.0001, 0.001, 0.01)
# costs <- c(0.1, 1, 5, 10, 20)
# tab_list <- list()
# for (g in gammas){
#   for(c in costs){
#     model=svm(train_IPSSR~., data=data.frame(train_set), kernel='linear',
#               gamma= g, cost = c) 
#     
#     misclf_tab = table(true=test_IPSSR, pred=predict (model,
#                                                       newdata =test_set)) 
#     tab1<-misclf_tab
#     print(c(g,c))
#     print(tab1)
#     
#   }
# }
model=svm(train_IPSSR~., data=data.frame(train_set), kernel='linear',
            gamma= 0.5, cost = 10)  

misclf_tab = table(true=test_IPSSR, pred=predict (model,
                                               newdata =test_set)) 
tab1<-misclf_tab
tab1
     
pred=predict (model,
              newdata = pred_set)
pred



# The predicted IPSSR are all equal to 4 or 5. We can justify that noticing that all the 63
# patients for which we predicted the values are from the HSCT cohort

idxs <- which(is.na(IPSSR_red))
NA_infos <- pca_data[idxs,]
table(NA_infos$cohort)





# ## TRYING TO ENRICH THE PREDICTION SET (DIDN'T WORK) ------
# 
# attach(new_MDS)
# 
# ## just chekcing the correspondence between IPSSR and genomic groups
# # gg01 <- which(genomic_group==0 | genomic_group==1 | genomic_group==6)
# # table(IPSSR[gg01])
# # gg24 <- which(genomic_group==2 | genomic_group==4)
# # table(IPSSR[gg24])
# 
# ## we proceed with the results of PCA - the reduced data of 12 Comp
# 
# drops <- c("age_AML","MAC_reg", "donor_type", "HLA_mismatch",
#            "relapse", "wait_time", "age", "n_mutgenes")
# data_red <- new_MDS[ , !(names(new_MDS) %in% drops)]
# head(data_red)
# 
# 
# bmb_rm <- which(bmblasts <=30) #the records of bmblasts that are out of range (corrputed data)
# data_red <- data_red[bmb_rm,]
# tmp_isna <- data_red[is.na(data_red$IPSSR), ] # patients without IPSSR (need to do prediction)
# # but we have missing values for other columns!
# data_red_noIPSSR <-  tmp_isna[ , !(names(tmp_isna) %in% c("IPSSR"))] # construct the dataset without the IPSSR
# tmp_omit <- na.omit(data_red_noIPSSR) # we check how many rows are NaN-free
# nan_report(data_red_noIPSSR, tmp_omit) # we get only 63 such rows which means we can perform prediction only for this few records
# 
# ## checking if we can get bigger data with removing other columns
# sumisna <- function(col){
#   return(sum(is.na(col)))
# }
# sapply(data_red_noIPSSR, sumisna) 
# ## We tried to remove age to see if we get more data, but we're still stuck with 63 rows for prediction
# data_red_noIPSSR_noage <- data_red_noIPSSR[ , !(names(data_red_noIPSSR) %in% c("agecat69"))]
# sapply(data_red_noIPSSR_noage, sumisna)
# tmp_omit <- na.omit(data_red_noIPSSR_noage)
# nan_report(data_red_noIPSSR_noage, tmp_omit)
# ## so we proceed with data_red_noIPSSR
# 
# 
# ## REGRESSION TREES (NOT WORKING :( ) ------- 
# 
# ## HGB
# 
# attach(new_MDS)
# data_tree <- data.frame(gender, genomic_group, bmblasts, age, AML_phase, hgb)
# bmb_rm <- which(bmblasts <=30) #the records of bmblasts that are out of range (corrputed data)
# data_tree <- data_tree[bmb_rm,]
# row.names(data_tree) <- NULL    #restarting the indexes
# 
# data_tree_noHGB <- na.omit(data_tree[, !(names(data_tree) %in% c("hgb"))]) # train test and pred sets
# dim(data_tree_noHGB)[1] #2616
# hgb_idxs  <- as.numeric(row.names(data_tree_noHGB)) 
# labels_hgb <- data_tree$hgb[hgb_idxs] # train test and pred index 2616
# 
# data_full_HGB <- data.frame(data_tree_noHGB, labels_hgb) # train test and pred data for HGB 2616
# detach(new_MDS)
# attach(data_full_HGB)
# train_test_idxs_hgb <- which(!is.na(labels_hgb))
# length(train_test_idxs_hgb)[1] #2281
# 
# train_test_set_hgb <- data_full_HGB[train_test_idxs_hgb,] #2281
# pred_set_hgb <- data_full_HGB[-train_test_idxs_hgb,] #335
# set.seed(1)
# train_size_hgb <-2100
# train_idxs_hgb <- sample(length(train_test_idxs_hgb), train_size_hgb)
# train_set_hgb <- train_test_set_hgb[train_idxs_hgb,]
# test_set_hgb <- train_test_set_hgb[-train_idxs_hgb,]
# 
# dim(train_set_hgb) #2100
# dim(test_set_hgb) #181
# 
# train_labels_hgb <- train_test_set_hgb$labels_hgb[train_idxs_hgb]
# test_labels_hgb <- train_test_set_hgb$labels_hgb[-train_idxs_hgb]
# 
# length(train_labels_hgb) #2100
# length(test_labels_hgb) #181
# 
# detach(data_full_HGB)
# attach(train_set_hgb) # dont forget !!!
# model_hgb <- tree(train_labels_hgb ~., train_set_hgb[,!(names(train_set_hgb) %in% c("labels_hgb"))])
# model_hgb <- tree(train_labels_hgb ~., train_set_hgb[,!(names(train_set_hgb) %in% c("labels_hgb", "bmblasts", "age"))])
# 
# summary(model_hgb)
# draw.tree(model_hgb)
# 
# detach(train_set_hgb)
# ## LEU
# 
# attach(new_MDS)
# data_tree <- data.frame(genomic_group, bmblasts, age, AML_phase, leu)
# bmb_rm <- which(bmblasts <=30) #the records of bmblasts that are out of range (corrputed data)
# data_tree <- data_tree[bmb_rm,]
# row.names(data_tree) <- NULL    #restarting the indexes
# 
# data_tree_noleu <- na.omit(data_tree[, !(names(data_tree) %in% c("leu"))]) # train test and pred sets
# dim(data_tree_noleu)[1] #2616
# leu_idxs  <- as.numeric(row.names(data_tree_noleu)) 
# labels_leu <- data_tree$leu[leu_idxs] # train test and pred index 2616
# 
# data_full_leu <- data.frame(data_tree_noleu, labels_leu) # train test and pred data for leu 2616
# detach(new_MDS)
# attach(data_full_leu)
# train_test_idxs_leu <- which(!is.na(labels_leu))
# length(train_test_idxs_leu)[1] #2281
# 
# train_test_set_leu <- data_full_leu[train_test_idxs_leu,] #2281
# pred_set_leu <- data_full_leu[-train_test_idxs_leu,] #335
# set.seed(1)
# train_size_leu <-2000
# train_idxs_leu <- sample(length(train_test_idxs_leu), train_size_leu)
# train_set_leu <- train_test_set_leu[train_idxs_leu,]
# test_set_leu <- train_test_set_leu[-train_idxs_leu,]
# 
# dim(train_set_leu) #2100
# dim(test_set_leu) #181
# 
# train_labels_leu <- train_test_set_leu$labels_leu[train_idxs_leu]
# test_labels_leu <- train_test_set_leu$labels_leu[-train_idxs_leu]
# 
# length(train_labels_leu) #2100
# length(test_labels_leu) #181
# 
# detach(data_full_leu)
# attach(train_set_leu) # dont forget !!!
# model_leu <- tree(train_labels_leu ~., train_set_leu[,!(names(train_set_leu) %in% c("labels_leu"))])
# model_leu <- tree(train_labels_leu ~., train_set_leu[,!(names(train_set_leu) %in% c("labels_leu", "genomic_group", "gender", "AML_phase"))])
# 
# summary(model_leu)
# draw.tree(model_leu)
# 
# ## PREDICT THE VALUES FOR MISSING BLOOD CHARACTERISTICS WITH LINEAR REGRESSION (DIDN'T WORK) ---------
# 
# ## HGB
# attach(new_MDS)
# bmb_rm <- which(bmblasts <=30) #the records of bmblasts that are out of range (corrputed data)
# data_red <- new_MDS[bmb_rm,]
# detach(new_MDS)
# attach(data_red)
# 
# 
# 
# 
# data_linreg_hgb <- data.frame(gender, genomic_group, bmblasts, age, AML_phase, hgb)
# row.names(data_linreg_hgb) <- NULL    #restarting the indexes
# head(data_linreg_hgb)
# 
# data_linreg_noHGB <- na.omit(data_linreg_hgb[, !(names(data_linreg_hgb) %in% c("hgb"))]) # train test and pred sets
# dim(data_linreg_noHGB)[1] #2616
# hgb_idxs  <- as.numeric(row.names(data_linreg_noHGB)) 
# labels_hgb <- data_linreg_hgb$hgb[hgb_idxs] # train test and pred index 2616
# labels_hgb <- labels_hgb[which(!is.na(labels_hgb))]
# 
# detach(data_red)
# attach(data_linreg_noHGB)
# 
# 
# 
# gender_ = ifelse(gender=='F', 1, 0)
# ggroup1 = ifelse(genomic_group=='1', 1, 0)
# ggroup2 = ifelse(genomic_group=='2', 1, 0)
# ggroup3 = ifelse(genomic_group=='3', 1, 0)
# ggroup4 = ifelse(genomic_group=='4', 1, 0)
# ggroup5 = ifelse(genomic_group=='5', 1, 0)
# ggroup6 = ifelse(genomic_group=='6', 1, 0)
# ggroup7 = ifelse(genomic_group=='7', 1, 0)
# cat_vars <- cbind(gender_, ggroup1, ggroup2, ggroup3, ggroup4, ggroup5, 
#                   ggroup6, ggroup7, AML_phase)
# 
# 
# ## Normalizing the numerical and categorical variables (FAMD algorithm)
# # cat_vars <- c("IPSSR", "genomic_group", "KT_complex","cyto_trt", 
# #               "fup_status", "AML_phase", "agecat69", "gender", "cohort")
# 
# # num_vars <- c(hgb, leu, neut, plt, bmblasts, n_mutgenes, n_mutations, time)
# num_vars_scaled <- scale(cbind(bmblasts, age))
# 
# scale_cat <- function(col){
#   return(col/sqrt((sum(col)/length(col))))
# }
# # scale_cat(gender_)
# gender_scaled <- scale_cat(gender_)
# ggroup1_scaled <- scale_cat(ggroup1)
# ggroup2_scaled <- scale_cat(ggroup2)
# ggroup3_scaled <- scale_cat(ggroup3)
# ggroup4_scaled <- scale_cat(ggroup4)
# ggroup5_scaled <- scale_cat(ggroup5)
# ggroup6_scaled <- scale_cat(ggroup6)
# ggroup7_scaled <- scale_cat(ggroup7)
# AML_phase_scaled <- scale_cat(AML_phase)
# cat_vars_scaled <- cbind(gender_scaled, ggroup1_scaled, 
#                          ggroup2_scaled, ggroup3_scaled, ggroup4_scaled,
#                          ggroup5_scaled, ggroup6_scaled, ggroup7_scaled, 
#                         AML_phase_scaled)
# 
# ## create a new dataset from scaled variables
# scaled_data_hgb <- data.frame(num_vars_scaled, cat_vars_scaled)
# head(scaled_data_hgb)
# 
# 
# 
# 
# 
# 
# train_test_idxs_hgb <- which(!is.na(labels_hgb))
# length(train_test_idxs_hgb)[1] #2281
# 
# train_test_set_hgb <- scaled_data_hgb[train_test_idxs_hgb,] #2281
# pred_set_hgb <- scaled_data_hgb[-train_test_idxs_hgb,] #335
# set.seed(1)
# train_size_hgb <-2100
# train_idxs_hgb <- sample(length(train_test_idxs_hgb), train_size_hgb)
# train_set_hgb <- train_test_set_hgb[train_idxs_hgb,]
# test_set_hgb <- train_test_set_hgb[-train_idxs_hgb,]
# 
# dim(train_set_hgb) #2100
# dim(test_set_hgb) #181
# 
# train_labels_hgb <- labels_hgb[train_idxs_hgb]
# test_labels_hgb <- labels_hgb[-train_idxs_hgb]
# 
# length(train_labels_hgb) #2100
# length(test_labels_hgb) #181
# 
# detach(data_linreg_hgb)
# attach(train_set_hgb) # dont forget !!!
# 
# 
# pca_hgb <- princomp(train_set_hgb, scores = T)
# load.data <- pca_hgb$loadings
# 
# summary(pca_hgb) # 8 components needed for Explained Variance ~ 0.87 (the 8th is the elbow)
# 
# #plot the loadings
# x11()
# nk <- 8 # number of compoments to use
# par(mar = c(1,1,1,1), mfrow = c(ncol(load.data[,c(1:nk)]),1))
# for(i in 1:nk) barplot(load.data[,i], ylim = c(-1, 1))
# 
# k = 8                        # number of columns to explain at least 80% of total variance
# pca_data_red_hgb <- data.frame(pca_hgb$scores[,1:k])  # reduced dataset
# 
# model_hgb <- lm(train_labels_hgb~., pca_data_red_hgb )
# summary(model_hgb)
# 
# A=rbind(c(0,0,0,1,0,0,0,0,0),c(0,0,0,0,1,0,0,0,0), c(0,0,0,0,0,1,0,0,0), c(0,0,0,0,0,0,0,1,0))
# b=c(0,0,0,0)
# linearHypothesis(model_hgb, A, b)   
# 
# # we cut the components  3,4,5 and 7 basing on the previous test (pval = 0.8761)
# var_to_cut <- c(3,4,5,7)
# 
# model_hgb <- lm(train_labels_hgb~., pca_data_red_hgb[,-var_to_cut] )
# summary(model_hgb)
# 
# # we cut also the component 2 (pval = 0.10201)
# var_to_cut <- c(2,3,4,5,7)
# 
# model_hgb <- lm(sqrt(train_labels_hgb)~., pca_data_red_hgb[,-var_to_cut] )
# summary(model_hgb)   # now seems to be fine
# 
# x11()
# par(mfrow=c(2,2))
# plot(model_hgb)   # nice plots (but don't look to the p-value please)
# 
# shapiro.test(model_hgb$residuals)
# 
# library(xgboost)
# 
# xgbst <- xgboost(as.matrix(pca_data_red_hgb[,-var_to_cut]), 
#                  train_labels_hgb, objective="reg:squarederror", nrounds = 2, eval='r2')
# 
# summary(xgbst)
# detach(train_set_hgb)
# 
# ## REPLACE HGB, LEU, PLT, NEUT WITH MEAN VALUES (NOTHING ELSE WORKED) -> tmp_MDS.RData---------------
# attach(new_MDS)
# 
# # hgb 
# dnh <- which(cohort=="DNH")
# hsct <- which(cohort =="HSCT")
# 
# hgb_HSCT <- hgb[hsct]
# hgb_DNH <- hgb[dnh]
# 
# avail_hgb_HSCT <- hgb_HSCT[which(!is.na(hgb_HSCT))]
# avail_hgb_DNH <- hgb_DNH[which(!is.na(hgb_DNH))]
# 
# new_MDS$hgb[which(is.na(hgb) & cohort=="DNH")] <- mean(avail_hgb_DNH)
# new_MDS$hgb[which(is.na(hgb) & cohort=="HSCT")] <- mean(avail_hgb_HSCT)
# 
# 
# 
# # plt 
# plt_HSCT <- plt[hsct]
# plt_DNH <- plt[dnh]
# 
# avail_plt_HSCT <- plt_HSCT[which(!is.na(plt_HSCT))]
# avail_plt_DNH <- plt_DNH[which(!is.na(plt_DNH))]
# 
# new_MDS$plt[which(is.na(plt) & cohort=="DNH")] <- mean(avail_plt_DNH)
# new_MDS$plt[which(is.na(plt) & cohort=="HSCT")] <- mean(avail_plt_HSCT)
# 
# 
# 
# 
# 
# # leu
# leu_HSCT <- leu[hsct]
# leu_DNH <- leu[dnh]
# 
# avail_leu_HSCT <- leu_HSCT[which(!is.na(leu_HSCT))]
# avail_leu_DNH <- leu_DNH[which(!is.na(leu_DNH))]
# 
# new_MDS$leu[which(is.na(leu) & cohort=="DNH")] <- mean(avail_leu_DNH)
# new_MDS$leu[which(is.na(leu) & cohort=="HSCT")] <- mean(avail_leu_HSCT)
# 
# 
# 
# 
# # neut
# neut_HSCT <- neut[hsct]
# neut_DNH <- neut[dnh]
# 
# avail_neut_HSCT <- neut_HSCT[which(!is.na(neut_HSCT))]
# avail_neut_DNH <- neut_DNH[which(!is.na(neut_DNH))]
# 
# new_MDS$neut[which(is.na(neut) & cohort=="DNH")] <- mean(avail_neut_DNH)
# new_MDS$neut[which(is.na(neut) & cohort=="HSCT")] <- mean(avail_neut_HSCT)
# 
# 
# tmp_MDS <- new_MDS[which(bmblasts<=30 & !is.na(bmblasts)),]
# save(tmp_MDS, file = "tmp_MDS.RData")
# sum(is.na(tmp_MDS$hgb))
# sum(is.na(tmp_MDS$plt))
# sum(is.na(tmp_MDS$leu))
# sum(is.na(tmp_MDS$neut))
# 
# 
# 
# 
# 
# 
## REPEAT CLASSIFICATION FOR THE MISSING VALUES OF IPSSR -----------


# Preprocessing the new dataset (scaling according to the training)

attach(tmp_MDS)

tmp_MDS$agecat69 <- cut(tmp_MDS$age, breaks=c(0, 69, Inf), labels=c(1, 0))

nan_report <- function(df1, df2){
  cat("Initial Lenght: ", nrow(df1))
  cat("\nLenght after removing NaN Values: ", nrow(df2))
}

drops <- c("age_AML","MAC_reg", "donor_type", "HLA_mismatch", 
           "relapse", "wait_time", "age")
data_red_tmp <- tmp_MDS[ , !(names(tmp_MDS) %in% drops)]

gender_ = as.numeric(ifelse(tmp_MDS$gender=='F', 1, 0))
cohort_ = as.numeric(ifelse(tmp_MDS$cohort=='HSCT', 1, 0))
ggroup1 = as.numeric(ifelse(tmp_MDS$genomic_group=='1', 1, 0))
ggroup2 = as.numeric(ifelse(tmp_MDS$genomic_group=='2', 1, 0))
ggroup3 = as.numeric(ifelse(tmp_MDS$genomic_group=='3', 1, 0))
ggroup4 = as.numeric(ifelse(tmp_MDS$genomic_group=='4', 1, 0))
ggroup5 = as.numeric(ifelse(tmp_MDS$genomic_group=='5', 1, 0))
ggroup6 = as.numeric(ifelse(tmp_MDS$genomic_group=='6', 1, 0))
ggroup7 = as.numeric(ifelse(tmp_MDS$genomic_group=='7', 1, 0))
cyto = as.numeric(ifelse(tmp_MDS$cyto_trt == 'yes', 1, 0))

# note: the following command works properly because no values are missing in the entire dataset.
#       indeed hgb, leu, neut and plt were filled, bmblasts were removed, while for the others 3 variables
#       we had no missing values from the beginning
tmp_num_vars <- cbind(tmp_MDS$hgb, tmp_MDS$leu, tmp_MDS$neut, tmp_MDS$plt, tmp_MDS$bmblasts, 
                  tmp_MDS$n_mutgenes, tmp_MDS$n_mutations, tmp_MDS$time)


# Here it's not the case for the last four variables, so we have to remove the missing values
tmp_cat_vars <- cbind(gender_, cohort_, ggroup1, 
                  ggroup2, ggroup3, ggroup4,
                  ggroup5, ggroup6, ggroup7, 
                  cyto, tmp_MDS$KT_complex,
                  tmp_MDS$fup_status, tmp_MDS$AML_phase, as.numeric(tmp_MDS$agecat69))

tmp_cat_vars[which(tmp_cat_vars[,14]==2),14]=0  # wtf?



idxs_not_NA <- which(!is.na( as.numeric(tmp_MDS$KT_complex) + as.numeric(tmp_MDS$fup_status) + as.numeric(tmp_MDS$AML_phase) + as.numeric(tmp_MDS$agecat69)))

tmp_cat_vars_red <- tmp_cat_vars[idxs_not_NA,]  # final dataset of categorical variables
IPSSR_final <- tmp_MDS$IPSSR[idxs_not_NA]
tmp_num_vars_red <- tmp_num_vars[idxs_not_NA,]  # final dataset of numerical variables


scaled.center.num
scaled.scale.num
scaled.scale.cat

tmp_num_vars_scaled <- tmp_num_vars_red
tmp_cat_vars_scaled <- tmp_cat_vars_red


for ( i in c(1:dim(tmp_num_vars_scaled)[2])) {
  tmp_num_vars_scaled[,i] <- (tmp_num_vars_scaled[,i] - scaled.center.num[i])/scaled.scale.num[i]
}


for ( i in c(1:dim(tmp_cat_vars_scaled)[2])) {
    
  tmp_cat_vars_scaled[,i] <- tmp_cat_vars_scaled[,i]/scaled.scale.cat[i]
  
}

# transform the data with the pca loadings obtained before (load.data)
data_pred <- cbind(tmp_num_vars_scaled, tmp_cat_vars_scaled)
dim(data_pred) #2333x22
k <- 22
pca_transform_matrix <- load.data[,1:k]
dim(pca_transform_matrix) #22x12

data_pred_transformed <- data_pred %*% pca_transform_matrix
dim(data_pred_transformed) #2333x12


# make predictions with the SVM clf we've trained
pred_IPSSR <- predict(model, newdata = data_pred_transformed[which(is.na(IPSSR_final)),])

IPSSR_final[which(is.na(IPSSR_final))] = pred_IPSSR
final_MDS <- tmp_MDS[idxs_not_NA,]
final_MDS$IPSSR <- IPSSR_final


# For Manfred and Luca

cat_vars_names <- (names(final_MDS) %in% c("genomic_group", "KT_complex", "cyto_trt", 
                    "donor_type", "HLA_mismatch", "MAC_reg",
                    "relapse", "AML_phase", 'agecat69')) + 0.0

sum(is.na(final_MDS$HLA_mismatch))

for (i in c(1:length(cat_vars_names))) {
  if(cat_vars_names[i] == 1 & sum(is.na(final_MDS[,i]))>0) {
    final_MDS[,i] <- sapply(final_MDS[,i], as.character)
    final_MDS[which(is.na(final_MDS[,i])),i] = "none"
  }
}

sum(is.na(final_MDS$HLA_mismatch))
sum(final_MDS$HLA_mismatch=="none")

