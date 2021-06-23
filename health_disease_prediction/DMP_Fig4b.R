# ==================================================
# By: A.Kurilshikov, R.Gacesa, UMCG (2021)
#
# Script for plotting Fig 4/b from summary statistics
# provided in this github repo file: 
#  DMP/health_disease_prediction/diseases_signature_correlation.txt

# ==================================================

## libraries and functions
library(corrplot)

# set WD:
setwd('.')
# plot 
# =============
# load data
plotCorMat <- read.table('health_disease_prediction/diseases_signature_correlation.txt',header = T,row.names = 1,stringsAsFactors = F)
# set colors
col1<-colorRampPalette(c( "#9B2226", "#AE2021","#BB3E03","#CA6702","#FFFFFF",
                          "#90e0ef", "#0077b6","#023e8a","#001219"))
# rename
rownames(plotCorMat) = c("1.Blood.Anemia",
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
colnames(plotCorMat) = 1:37
# plot
pdf(file = "health_disease_prediction/DMP_Fig4b.pdf",width = 9,height = 9)
corrplot(as.matrix(plotCorMat),method = "square",tl.cex = 0.85,cl.cex = 0.85,col = col1(20),tl.col = "black")
dev.off()

