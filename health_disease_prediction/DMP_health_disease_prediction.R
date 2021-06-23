# ==================================================
# By: A.Kurilshikov, R.Gacesa, UMCG (2021)
# DMP project, training and testing of
# models for prediction of health and diseases
#
# NOTE: 
# - code is implemented on Mock data, 
# real DMP data can be obtained from EGA and Lifelines, 
# please see the DMP manuscript for details
# ==================================================

## libraries and functions
library(glmnet)
library(pROC)
library(vegan)
library(doSNOW)
library(corrplot)
library(foreach)

## function do_clr_extrnalWeighting
#
# function performs CLR transformation 
#  using external weighting
# ===========================================
do_clr_externalWeighting = function(interest_matrix, core_matrix){
  if(any(interest_matrix==0)) interest_matrix = interest_matrix + min(interest_matrix[interest_matrix>0])/2
  if(any(core_matrix==0)) core_matrix = core_matrix + min(core_matrix[core_matrix>0])/2
  
  # estimate weighting parameter
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x), na.rm=na.rm) / length(x))
  }
  Gmean_core = apply(core_matrix, 1, gm_mean)
  
  #do transformation
  data_prepared = cbind(Gmean_core,interest_matrix)
  data_transformed = t(apply(data_prepared,1,function(x){
    log(x / x[1])[-1]
  }))
  colnames(data_transformed) = colnames(data_transformed)
  rownames(data_transformed) = rownames(data_transformed)
  data_transformed
}

# ===========================
# ===========================
# MAIN
# ===========================
# ===========================

print(' >> starting microbiome-disease prediction')
# setwd: should be the root of github repo (DMP folder)
# example: 
#setwd('C:/Users/Ranko/Documents/ub_shared/dag/gits/git_23_06/DMP') 
setwd('.')

print(' >> loading data ...')
# load all necessary data
covar = read.table("Mock_data/covariates.txt")
Ncovar = ncol(covar)

taxa = read.table("Mock_data/taxa.txt")
shannon.div = diversity(taxa[,grep("[.]s__",colnames(taxa))],index = "shannon")
pathways = read.table("Mock_data/pathways.txt")

# select diseases
pheno_disease = read.table("Mock_data/diseases.txt")
Ndisease = ncol(pheno_disease) - 1 # the last disease is 'no disease' phenotype

print(' >> normalizing microbiome data')
# transform data using CLR
taxa_transformed = do_clr_externalWeighting(taxa,taxa[,grep("[.]s__",colnames(taxa))])
taxa_transformed = taxa_transformed[,colSums(taxa>0)>nrow(taxa) * 0.05]
pathways_transformed = do_clr_externalWeighting(pathways,pathways)
pathways_transformed = pathways_transformed[,colSums(pathways>0)>nrow(pathways) * 0.05]

## split data into training and test sets in cross-validation setup
print(' >> prepping training and test sets')
set.seed(12348)
randomized_samples = sample(1:nrow(taxa_transformed))
batches = split(randomized_samples,cut(seq_along(randomized_samples),5,labels = F))

data.pred = data.frame(covar,shannon = shannon.div,taxa_transformed,pathways_transformed)

# training/test sets, predictors
train_X = lapply(batches, function(x) data.pred[-x,])
test_X = lapply(batches, function(x) data.pred[x,])

#training/test sets, with microbiome data nullified, used in calculation of AUCs             
train_X.clin = lapply(train_X, function(x) {x[,4:ncol(x)] = 0;x})
test_X.clin = lapply(test_X, function(x) {x[,4:ncol(x)] = 0;x})

#training/test sets, with antropometric data nullified, used in calculation of AUCs                
train_X.microb = lapply(train_X, function(x) {x[,1:3] = 0;x})
test_X.microb = lapply(test_X, function(x) {x[,1:3] = 0;x})

#training/test sets, outcomes
train_Y = lapply(batches, function(x) pheno_disease[-x,])
test_Y = lapply(batches, function(x) pheno_disease[x,])
                

## Build prediction models
# =================================================================
all_prediction_models = list()
print(' >> training models')
for(i in 1:5) {
  print(paste0('   > X-val set ',i,' ...'))
  for(j in 1:ncol(pheno_disease)){
    print(paste0('    > Disease ',j,' ...'))
    current_X = train_X[[i]]
    current_Y = train_Y[[i]]
    current_y = current_Y[,j]
    #all_prediction_models[["pred.batch",i,".pheno.", j]] = cv.glmnet(as.matrix(complete.X),complete.Y,alpha = 0.5,nfolds = 10,family = "binomial")
    all_prediction_models[[paste0("pred.batch",i,".pheno.", j)]] = cv.glmnet(x = as.matrix(current_X),
                                                                     y = current_y,
                                                                     alpha = 0.5,nfolds = 10,family = "binomial")
  }
}

## calculate training set and test set predictions
# =========================================                
# predictions on training set, complete model (microbiome + clinical parameters)
train.predictions.completeModel = list()
for (i in 1:5){
  train.predictions.completeModel[[i]] = foreach (j = 1:37,.combine = cbind)%do%{
    # predict(all_prediction_models[[paste0("pred.batch",i,".pheno.",j,".RData")]],
    #          newx = as.matrix(train_X[[i]]),s="lambda.min")[,1]
    predict(all_prediction_models[[paste0("pred.batch",i,".pheno.",j)]],
            newx = as.matrix(train_X[[i]]),s="lambda.min")[,1]
    
    
  }
}
# predictions on training set, clinical model (age, sex, BMI only) 
train.predictions.clin = list()
for (i in 1:5){
  train.predictions.clin[[i]] = foreach (j = 1:37,.combine = cbind)%do%{
    # predict(all_prediction_models[[paste0("pred.batch",i,".pheno.",j,".RData")]],
    #         newx = as.matrix(train_X.clin[[i]]),s="lambda.min")[,1]
    predict(all_prediction_models[[paste0("pred.batch",i,".pheno.",j)]],
            newx = as.matrix(train_X.clin[[i]]),s="lambda.min")[,1]
  }
}
# predictions on training set, microbiome only model 
train.predictions.microb = list()
for (i in 1:5){
  train.predictions.microb[[i]] = foreach (j = 1:37,.combine = cbind)%do%{
    # predict(all_prediction_models[[paste0("pred.batch",i,".pheno.",j,".RData")]],
    #         newx = as.matrix(train_X.microb[[i]]),s="lambda.min")[,1]
    predict(all_prediction_models[[paste0("pred.batch",i,".pheno.",j)]],
            newx = as.matrix(train_X.microb[[i]]),s="lambda.min")[,1]
    
  }
}   

# predictions on test set
# =================================================
# > complete model
test.predictions.completeModel = list()
for (i in 1:5){
  test.predictions.completeModel[[i]] = foreach (j = 1:37,.combine = cbind)%do%{
    # predict(all_prediction_models[[paste0("pred.batch",i,".pheno.",j,".RData")]],
    #         newx = as.matrix(test_X[[i]]),s="lambda.min")[,1]
    predict(all_prediction_models[[paste0("pred.batch",i,".pheno.",j)]],
            newx = as.matrix(test_X[[i]]),s="lambda.min")[,1]
  }
}
# clinical characteristics model
test.predictions.clin = list()
for (i in 1:5){
  test.predictions.clin[[i]] = foreach (j = 1:37,.combine = cbind)%do%{
    # predict(all_prediction_models[[paste0("pred.batch",i,".pheno.",j,".RData")]],
    #         newx = as.matrix(test_X.clin[[i]]),s="lambda.min")[,1]
    predict(all_prediction_models[[paste0("pred.batch",i,".pheno.",j)]],
            newx = as.matrix(test_X.clin[[i]]),s="lambda.min")[,1]
  }
}
# microbiome model
test.predictions.microb = list()
for (i in 1:5){
  test.predictions.microb[[i]] = foreach (j = 1:37,.combine = cbind)%do%{
    # predict(all_prediction_models[[paste0("pred.batch",i,".pheno.",j,".RData")]],
    #         newx = as.matrix(test_X.microb[[i]]),s="lambda.min")[,1]
    predict(all_prediction_models[[paste0("pred.batch",i,".pheno.",j)]],
            newx = as.matrix(test_X.microb[[i]]),s="lambda.min")[,1]
  }
}   

## calculate AUCs 
# =================================================
auc.train.completeModel = foreach(i = 1:5,.combine = rbind)%:%
  foreach(j = 1:37,.combine = cbind) %do%{
    auc(roc(train_Y[[i]][,j],train.predictions.completeModel[[i]][,j],direction = "<",levels = c(0,1)))
  }
auc.test.completeModel = foreach(i = 1:5,.combine = rbind)%:%
  foreach(j = 1:37,.combine = cbind) %do%{
    auc(roc(test_Y[[i]][,j],test.predictions.completeModel[[i]][,j],direction = "<",levels = c(0,1)))
  }

#clinical prediction
auc.train.clin = foreach(i = 1:5,.combine = rbind)%:%
  foreach(j = 1:37,.combine = cbind) %do%{
    auc(roc(train_Y[[i]][,j],train.predictions.clin[[i]][,j],direction = "<",levels = c(0,1)))
  }
auc.test.clin = foreach(i = 1:5,.combine = rbind)%:%
  foreach(j = 1:37,.combine = cbind) %do%{
    auc(roc(test_Y[[i]][,j],test.predictions.clin[[i]][,j],direction = "<",levels = c(0,1)))
  }

#microbial prediction
auc.train.microb = foreach(i = 1:5,.combine = rbind)%:%
  foreach(j = 1:37,.combine = cbind) %do%{
    auc(roc(train_Y[[i]][,j],train.predictions.microb[[i]][,j],direction = "<",levels = c(0,1)))
  }
auc.test.microb = foreach(i = 1:5,.combine = rbind)%:%
  foreach(j = 1:37,.combine = cbind) %do%{
    auc(roc(test_Y[[i]][,j],test.predictions.microb[[i]][,j],direction = "<",levels = c(0,1)))
  }

## Make signature correlation plot
# =================================================================

#merge all predictions together
all_diseases.test= Reduce(rbind,test_Y)
all_diseases.train= Reduce(rbind,train_Y)

all_predictions.train.microb = Reduce(rbind,train.predictions.microb)
all_predictions.train.clin = Reduce(rbind,train.predictions.clin)
all_predictions.train.completeModel = Reduce(rbind,train.predictions.completeModel)
all_predictions.test.microb = Reduce(rbind,test.predictions.microb)
all_predictions.test.clin = Reduce(rbind,test.predictions.clin)
all_predictions.test.completeModel = Reduce(rbind,test.predictions.completeModel)

# plot 

CorMat.forPlot = cor(pheno_disease,use = "pairwise.complete.obs")
CorMat.forPlot[lower.tri(CorMat.forPlot)] = cor(all_predictions.test.microb,use = "pairwise.complete.obs")[lower.tri(CorMat.forPlot)]

col1<-colorRampPalette(c( "#9B2226", "#AE2021","#BB3E03","#CA6702","#FFFFFF",
                          "#90e0ef", "#0077b6","#023e8a","#001219"))

rownames(CorMat.forPlot) = c("1.Blood.Anemia",
                             "2.Blood.Thrombosis",
                             "3.Cancer.Any",
                             "4.Cardiovascular.Arrythmia.MedDiagnosed",
                             "5.Cardiovascular.Colesterol.high",
                             "6.Cardiovascular.Heart.Attack",
                             "7.Cardiovascular.Heart.Failure.Disorder",
                             "8.Cardiovascular.Heartrate.complains",
                             "9.Cardiovascular.Hypertension",
                             "10.Endocrine.DiabetesT2",
                             "11.Gastrointestinal.Stomach.Ulcer",
                             "12.Hepatologic.Gallstones",
                             "13.Mental.Any",
                             "14.Mental.Burn.Out",
                             "15.Mental.Depression",
                             "16.Mental.Other.anxiety",
                             "17.Mental.Panic.disorder",
                             "18.Neurological.Dizziness.Falling",
                             "19.Neurological.Mental.Fibromyalgia",
                             "20.Neurological.Migraine",
                             "21.Other.Autoimmune.Rheumatoid.Artritis",
                             "22.Other.Chronic.cystitis",
                             "23.Other.Chronic.Inflammation.Throatnose",
                             "24.Other.Chronic.Muscle.Weakness",
                             "25.Other.Fractures",
                             "26.Other.Incontinence",
                             "27.Other.Kidney.Stones",
                             "28.Other.Osteoarthritis",
                             "29.Other.Osteoporosis",
                             "30.Other.RSI",
                             "31.Pulmonary.Autoimmune.Asthma",
                             "32.Pulmonary.COPD",
                             "33.Skin.Autoimmune.Atopic.dermatitis",
                             "34.Skin.Autoimmune.Psoriasis",
                             "35.Skin.Autoimmune.Severe.acne",
                             "36.Gastrointestinal.Rome3_IBS.Any",
                             "37.None.NoDiseases"
)
colnames(CorMat.forPlot) = 1:37
pdf(file = "health_disease_prediction/Mockdata_Fig4b.pdf",width = 9,height = 9)
corrplot(CorMat.forPlot,method = "square",col = col1(20),tl.cex = 0.85,cl.cex = 0.85,tl.col = "black")
dev.off()

