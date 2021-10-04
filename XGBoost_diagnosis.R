
library(caret)
library(corrplot)			# plot correlations
library(doParallel)		# parallel processing
library(dplyr)        # Used by caret		
library(pROC)				  # plot the ROC curve
library(xgboost)      # Extreme Gradient Boosting
library(verification)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("paste", "base")
conflict_prefer("as.matrix", "base")
conflict_prefer("rename", "dplyr")

my_comprisons=list(c("AMS4k","nAMS4k"))

#### MRM protein XGBoost ####

print(paste(my_comprisons[[1]],collapse = "_vs_"))

# input MRM rawdata (not log2rawdata) and q_value
load(paste("../data/",my_comprisons[[1]][1],"_",my_comprisons[[1]][2],"_MRM_raw_result.RData",sep=""))

print(1)

data <- MRM_raw_re  %>% filter(q_value < 0.05) %>%
  select(-c(log2FC,q_value)) %>% tibble::column_to_rownames( var ='symbol')

print(nrow(data))
td<- t(data)

td <-data.frame(td, stringsAsFactors = F, check.names = FALSE)
 
rawdata <- td %>% tibble::add_column(sampleID = rownames(td)) %>%
  mutate(group =sapply(strsplit(sampleID, "\\_"), `[`, 2 )) %>%
  subset(select= -c(sampleID)) %>% filter(!is.na(group))
rawdata$group <- factor(rawdata$group, levels= c(my_comprisons[[1]][1],my_comprisons[[1]][2]))

sub_raw <- rawdata[,c("group", rownames(data))]
colnames(sub_raw)[1] <- "target"


#xgboost model
set.seed(10)
trainIndex <- createDataPartition(sub_raw$target,p=.6,list=FALSE)
trainData_MRM <- sub_raw[trainIndex,]
testData_MRM  <- sub_raw[-trainIndex,]

trainX <-trainData_MRM[,!names(trainData_MRM) %in% c("target")]   # Pull out the dependent variable
testX <- testData_MRM[,!names(testData_MRM) %in% c("target")]

set.seed(301) 
ControlParamteres <- trainControl(method = "repeatedcv",
                                  number = 10,
                                  summaryFunction=twoClassSummary,
                                  savePredictions = TRUE,
                                  classProbs = TRUE)

parametersGrid <-  expand.grid(eta = c(0.01,0.02,0.05,0.06,0.1,0.15), #0.01 or
                               colsample_bytree=c(0.6,0.7,0.8),
                               max_depth=c(2,3,4),
                               nrounds=100,
                               gamma=c(0,1),
                               min_child_weight=c(1,2,3),
                               subsample=c(0.6,0.7,0.8,0.9)
)

set.seed(608)
modelxgboost_MRM <- train(target~., 
                          data = trainData_MRM,
                          method = "xgbTree",
                          metric = "ROC",
                          trControl = ControlParamteres,
                          tuneGrid=parametersGrid)



#modelxgboost_MRM
tmp_MRM <- varImp(modelxgboost_MRM)
imp_MRM <-data.frame(tmp_MRM$importance , stringsAsFactors = F) %>% 
  filter(Overall > 0) %>% arrange( var= "Overall" ) %>% 
  tibble::rownames_to_column(var = "symbol")
write.csv(imp_MRM, file = paste(paste(my_comprisons[[1]],collapse = "_vs_")," MRM variables important feature.csv",sep = ""))

pdf(file = paste(paste(my_comprisons[[1]],collapse = "_vs_"),"MRM variables important feature.pdf", sep = ""),
    width = 6, height = 4)

p_MRM <- ggplot(data= imp_MRM, aes(x=reorder(symbol,Overall),y=Overall ))+
  geom_bar(position="dodge",stat="identity",width = 0.5 ,fill ="#F3574A") + 
 xlab("Variables")+ 
  ylab(" Importance Score") + theme_bw()

p_MRM +theme(axis.title= element_text(size=14,  color="black", face= "bold"))+
  theme(plot.title= element_text(size=14,  color="black", face= "bold",  hjust=0.5))+
  theme(axis.text.x= element_text(size=14,color="black",  vjust=0.7, hjust=0.7,angle = 45))+
  labs(title="Proteins variables importance" )
dev.off()

predictions_MRM<-predict(modelxgboost_MRM,testData_MRM)
confusionMatrix(table(predictions_MRM,testData_MRM$target))  

predictions_MRM<-predict(modelxgboost_MRM,testData_MRM,type="prob")

xgb.ROC_MRM <- roc(predictor=predictions_MRM[,1],
                   response=testData_MRM$target,
                   levels=levels(as.factor(testData_MRM$target)))
xgb.ROC_MRM$auc 
levels(testData_MRM$target)<- c('1','0')

pvalue_MRM <- roc.area(as.numeric(as.vector(testData_MRM$target)),predictions_MRM[,1])$p.value

save.image(file = paste(paste(my_comprisons[[1]],collapse = "_vs_"),"_MRM.RData",sep = ""))

#### Clinical indexes XGBoosts ####
load(paste("../data/",my_comprisons[[1]][1],"_",my_comprisons[[1]][2],"_clinic_raw_result.RData",sep=""))
rawmetadata=read.csv("../data/mrm_sample_list_v4.csv",sep= ',', check.names = F)
clinc_info <- rawmetadata %>% select(c(sample_id, group)) 

merge_clinc_data <- clin_raw_re %>% filter(p_value < 0.05) %>%
  tibble::column_to_rownames(var ='clinic') %>% select(-c(log2FC,p_value, q_value))
merge_clinc_data_temp <- t(merge_clinc_data)
merge_clinc_data_temp <- data.frame(merge_clinc_data_temp, stringsAsFactors = F, check.names = F)
clinic_data_merge <- merge_clinc_data_temp %>% tibble::rownames_to_column(var ='sample_id') %>%
  left_join(clinc_info, by = 'sample_id') %>% tibble::column_to_rownames(var = 'sample_id') %>%
  select(c(group, colnames(merge_clinc_data_temp)))

clinic_data_merge$group <- factor(clinic_data_merge$group, levels=c(my_comprisons[[1]][1],my_comprisons[[1]][2] ))

colnames(clinic_data_merge)[1]="target"
#xgboost model
set.seed(301)
trainIndex <- createDataPartition(clinic_data_merge$target,p=.6,list=FALSE)
trainData_clin <- clinic_data_merge[trainIndex,]
testData_clin  <- clinic_data_merge[-trainIndex,]

trainX <-trainData_clin[,!names(trainData_clin) %in% c("target")]   # Pull out the dependent variable
testX <- testData_clin[,!names(testData_clin) %in% c("target")]

set.seed(301)
ControlParamteres <- trainControl(method = "repeatedcv",
                                  number = 10,
                                  summaryFunction=twoClassSummary,
                                  savePredictions = TRUE,
                                  classProbs = TRUE)
parametersGrid <-  expand.grid(eta = c(0.05,0.1,0.2,0.3,0.4,0.5), #0.01 or
                               colsample_bytree=c(0.6,0.7,0.8),
                               max_depth=c(2,3,4),
                               nrounds=100,
                               gamma=c(0,1),
                               min_child_weight=c(1,2,3),
                               subsample=c(0.6,0.7,0.8,0.9)
)

set.seed(123)


modelxgboost_clinic <- train(target~., 
                             data = trainData_clin,
                             method = "xgbTree",
                             metric = "ROC",
                             trControl = ControlParamteres,
                             tuneGrid=parametersGrid)

modelxgboost_clinic

tmp <- varImp(modelxgboost_clinic)
imp <-data.frame(tmp$importance , stringsAsFactors = F) %>% 
  filter(Overall > 0) %>% arrange( var= "Overall" ) %>% 
  tibble::rownames_to_column(var = "symbol")
write.csv(imp, file = paste(paste(my_comprisons[[1]],collapse = "_vs_")," clinic variables important feature.csv",sep = "") )

pdf(file = paste(paste(my_comprisons[[1]],collapse = "_vs_")," clinic variables important feature.pdf",sep = ""),
    width = 5, height = 6)
p_clinic <- ggplot(data= imp, aes(x=reorder(symbol,Overall),y=Overall ))+
  geom_bar(position="dodge",stat="identity",width = 0.5 , fill ="#01569B") + 
  coord_flip() + xlab("Variables")+ 
  ylab(" Importance Score") + theme_bw()

p_clinic +theme(axis.title= element_text(size=14,  color="black", face= "bold"))+
  theme(plot.title= element_text(size=14,  color="black", face= "bold",  hjust=0.5))+
  theme(axis.text= element_text(size=14,color="black",  vjust=0.5, hjust=0.5))+
  labs(title="Clinical indexes variables importance" )
dev.off()
predictions_clinic<-predict(modelxgboost_clinic,testData_clin)
confusionMatrix(table(predictions_clinic,testData_clin$target))  

predictions_clinic<-predict(modelxgboost_clinic,testData_clin,type="prob")

xgb.ROC_clinic<- roc(predictor=predictions_clinic[,1],
                     response=testData_clin$target,
                     levels=levels(as.factor(testData_clin$target)))
xgb.ROC_clinic$auc 
levels(testData_clin$target)<- c('1','0')

pvalue_clinc <- roc.area(as.numeric(as.vector(testData_clin$target)),predictions_clinic[,1])$p.value
save.image(file = paste(paste(my_comprisons[[1]],collapse = "_vs_"),"_clin.RData",sep = ""))

#### Combined data XGBoost ####

colnames(rawmetadata)
sub_raw2 <- sub_raw %>% tibble::rownames_to_column(var ='sampleID')
merge_1=merge(sub_raw2,rawmetadata,by="sampleID",all.x=T,sort = F)
#remove group
clinic_data_merge2 <- clinic_data_merge %>% select(-c(target)) %>% 
  tibble::rownames_to_column(var='sample_id')

merge_clinic_MRM <- clinic_data_merge2 %>% left_join(merge_1, by = 'sample_id') %>%
  tibble::column_to_rownames(var ='sample_id')%>%
  select(c(colnames(sub_raw2), colnames(clinic_data_merge))) %>% select(-c( sampleID))

#xgboost model
set.seed(301)
trainIndex <- createDataPartition(merge_clinic_MRM$target,p=.6,list=FALSE)
trainData_comb <- merge_clinic_MRM[trainIndex,]
testData_comb  <- merge_clinic_MRM[-trainIndex,]

trainX <-trainData_comb[,!names(trainData_comb) %in% c("target")]   # Pull out the dependent variable
testX <- testData_comb[,!names(testData_comb) %in% c("target")]

set.seed(301)
ControlParamteres <- trainControl(method = "repeatedcv",
                                  number = 10,
                                  summaryFunction=twoClassSummary,
                                  savePredictions = TRUE,
                                  classProbs = TRUE)
parametersGrid <-  expand.grid(eta = c(0.05,0.1,0.2,0.3,0.4,0.5), #0.01 or
                               colsample_bytree=c(0.6,0.7,0.8),
                               max_depth=c(2,3,4),
                               nrounds=100,
                               gamma=c(0,1),
                               min_child_weight=c(1,2,3),
                               subsample=c(0.6,0.7,0.8,0.9)
)

set.seed(123)

modelxgboost_clinic_MRM <- train(target~., 
                                 data = trainData_comb,
                                 method = "xgbTree",
                                 metric = "ROC",
                                 trControl = ControlParamteres,
                                 tuneGrid=parametersGrid)

#modelxgboost_clinic_MRM

tmp <- varImp(modelxgboost_clinic_MRM)
imp_comnine <-data.frame(tmp$importance , stringsAsFactors = F) %>% 
  filter(Overall > 0) %>% arrange( var= "Overall" ) %>% 
  tibble::rownames_to_column(var = "symbol")
write.csv(imp_comnine,file = paste(paste(my_comprisons[[1]],collapse = "_vs_")," clinic&MRM variables important feature.pdf",sep = "") )

pdf(file = paste(paste(my_comprisons[[1]],collapse = "_vs_")," clinic&MRM variables important feature.pdf",sep = ""),
    width = 5, height = 6)

p_clinic_MRM <- ggplot(data= imp_comnine, aes(x=reorder(symbol,Overall),y=Overall ))+
  geom_bar(position="dodge",stat="identity",width = 0.5 ,fill ="#FFA500") + 
  coord_flip() + xlab("Variables")+ 
  ylab(" Importance Score") + theme_bw()

p_clinic_MRM +theme(axis.title= element_text(size=14,  color="black", face= "bold"))+
  theme(plot.title= element_text(size=14,  color="black", face= "bold",  hjust=0.5))+
  theme(axis.text= element_text(size=14,color="black",  vjust=0.5, hjust=0.5))+
  labs(title="Combined variables importance")
dev.off()

predictions_clinic_MRM<-predict(modelxgboost_clinic_MRM,testData_comb)
confusionMatrix(table(predictions_clinic_MRM,testData_comb$target))  

predictions_clinic_MRM<-predict(modelxgboost_clinic_MRM,testData_comb,type="prob")

xgb.ROC_clinic_MRM <- roc(predictor=predictions_clinic_MRM[,1],
                          response=testData_comb$target,
                          levels=levels(as.factor(testData_comb$target)))
xgb.ROC_clinic_MRM$auc 

 levels(testData_comb$target)<- c('1','0')
 pvalue_combine <- roc.area(as.numeric(as.vector(testData_comb$target)),predictions_clinic_MRM[,1])$p.value
save.image(file = paste(paste(my_comprisons[[1]],collapse = "_vs_"),"_combine.RData",sep = ""))

xgb.ROC=list(xgb.ROC_clinic,xgb.ROC_MRM,xgb.ROC_clinic_MRM)
save(xgb.ROC,file=paste(paste(my_comprisons[[1]],collapse = "_vs_"),"_ROC.Rdata",sep = ""))


pvalue_MRM <- sprintf("%.4f", pvalue_MRM)
pvalue_clinc <- sprintf("%.4f", pvalue_clinc)
pvalue_combine <- sprintf("%.4f", pvalue_combine)
xgb.ROC_clinic_MRM$auc <- sprintf("%.2f", xgb.ROC_clinic_MRM$auc)
xgb.ROC_clinic$auc <- sprintf("%.2f", xgb.ROC_clinic$auc)
xgb.ROC_MRM$auc <- sprintf("%.2f", xgb.ROC_MRM$auc)

print(pvalue_MRM)
print(pvalue_clinc)
print(pvalue_combine)

pdf(file = paste(paste(my_comprisons[[1]],collapse = "_vs_"),"_ROC.pdf",sep = ""),
    width = 5.5, height = 5.5)

plot(xgb.ROC_MRM, main= "Diagnosis ROC curve", lwd=3 , col= "#F3574A")
p <- plot.roc(xgb.ROC_clinic, lwd = 3, add=TRUE, col="#01569B")

p2 <- plot.roc(xgb.ROC_clinic_MRM, lwd = 3, add=TRUE, col="#FFA500")

legend(0.8,0.3,  
       bty = "n",  
       title="",   
       legend=c(paste("AUC:",xgb.ROC_MRM$auc,"(",pvalue_MRM,") ,Proteins",sep = ""),
                paste("AUC:",xgb.ROC_clinic$auc,"(",pvalue_clinc,") ,Clinical indexes",sep = ""),
                paste("AUC:",xgb.ROC_clinic_MRM$auc,"(",pvalue_combine,") ,Combined",sep = "")),
       col=c("#F3574A","#01569B","#FFA500"), 
       lwd=3,
       y.intersp = 1)  

dev.off()     
save.image(file = "AMS4k_nAMS4k_3group.RData")
