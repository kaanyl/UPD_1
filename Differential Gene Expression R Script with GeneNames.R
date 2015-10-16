#Working Code for MArray CEL File Processig Dif Exp Analysis
######################################################################

targettt <- readTargets("targettt.txt", sep="")

library(affy)
library(hgu133plus2.db)
abatch <- ReadAffy()
eset <- rma(abatch)

f <- paste(targettt$Name,targettt$Target,sep="")
f <- factor(f)
f
design <- model.matrix(~0+f)
colnames(design) <- levels(f)
design
library(limma)
fit <- lmFit(eset, design)
names(fit)
cont.matrix <- makeContrasts(C="aC-bS", levels=design)
cont.matrix
fit2  <- contrasts.fit(fit, cont.matrix)
fit2  <- eBayes(fit2)


colnames(fit2)
topTable(fit2,coef=1)
# topTable(fit2,coef=2,adjust="fdr")
# topTable(fit2,coef=2)
# topTable(fit2,coef=2,adjust="fdr")

Selected  <- (fit2$p.value) <0.001
esetSel <- eset[Selected, ]
heatmap(exprs(esetSel))
unigeneTopTable <- topTable(fit2,coef=1,n=54675)

#Pathway and GO Analysis Gene Name Adding to Analysis 
#########################################################################

ls("package:hgu133plus2.db")
probes=row.names(unigeneTopTable)

Symbols = unlist(mget(probes, hgu133plus2GENENAME, ifnotfound=NA))
unigeneTopTable=cbind(probes,Symbols, unigeneTopTable)

ENTREZ_ID=unlist(mget(probes, hgu133plus2ENTREZID, ifnotfound=NA))
unigeneTopTable=cbind(probes, ENTREZ_ID, unigeneTopTable)

Genenames=unlist(mget(probes, hgu133plus2SYMBOL, ifnotfound=NA))
unigeneTopTable=cbind(probes, Genenames, unigeneTopTable)


#########################################################################


library(mlbench)
library(caret)
set.seed(998)
inTraining <- createDataPartition(BAA$ID, p = .75, list = FALSE)
training <- BAA[ inTraining,]
testing  <- BAA[-inTraining,]
##############
fitControl <- trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 10,
  ## repeated ten times
  repeats = 10)
##############
set.seed(825)
RfFit1 <- train(ID ~ ., data = training,
                method = "rf",
                trControl = fitControl,
                verbose = FALSE)
RfFit1
plot(RfFit1)

PRED<-predict(RfFit1,newdata=testing)


#SVM
##################################
set.seed(825)
svmFit <- train(ID ~ ., data = training,
                method = "svmRadial",
                trControl = fitControl,
                preProc = c("center", "scale"),
                tuneLength = 8)
svmFit


PRED<-predict(svmFit,newdata=testing)
PRED
