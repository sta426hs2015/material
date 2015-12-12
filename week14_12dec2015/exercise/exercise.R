### R code from vignette source 'exercise'

library(ALL)
data(ALL)
show(ALL)


 bCellSamples = grep("^B", ALL$BT)
 BcrAndNegSamples = which(ALL$mol.biol %in% c("BCR/ABL", "NEG"))
 samplesToUse = intersect(bCellSamples, BcrAndNegSamples)
 dataMatrix = exprs(ALL[ ,samplesToUse])
 classLabels = factor(ALL$mol.biol[samplesToUse])


 dim(dataMatrix)


##   highVarGenes = ...
##   dataMatrix = ...


highVarGenes = order(apply(dataMatrix, 1, var), decreasing=TRUE)[1:1000]
dataMatrix = dataMatrix[ highVarGenes, ]


library(e1071)


  model = svm(t(dataMatrix), classLabels, kernel = "linear")


predicted = predict(model, t(dataMatrix))
table(true = classLabels, pred = predicted)


model.cv = svm(t(dataMatrix), classLabels, kernel = "linear", cross = length(classLabels))
summary(model.cv)


model.cv$tot.accuracy


## library(tspair)


##   tspResult = tspcalc...


  library(tspair)
  tspResult = tspcalc(dataMatrix, classLabels)


 predictedLabels = predict(tspResult, dataMatrix)
 table(predictedLabels, classLabels)


##   for( i in 1:ncol(dataMatrix)){
##    ...
##   }


