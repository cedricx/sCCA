
## load original clinical data
med <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreezeDec2016/clinical/n1601_goassess_112_itemwise_vars_20161214.csv")
med.item <- med[,3:114]

## visualize the missing data
med.vis<-aggr(med.item,sortVars=TRUE,col=c('navyblue','red'),sortCombs = TRUE, numbers = TRUE, cex.numbers = 0.1)
med.item.rm <- med.item[,-which(colnames(med.item) == "scr008")]
## test MCAR
mcar.test <- TestMCARNormality(data = med.item.rm, del.lesscases = 1)
### MCAR is rejected

## impute the data using MICE (single thread)

numImp <- 15
med.item.imp.rm <- mice(med.item.rm, m = numImp, maxit = 20, method= "pmm",visitSequence = 'monotone',printFlag = FALSE)
med.item.imp.all <- mice(med.item, m = numImp, maxit = 20, method= "pmm",printFlag = FALSE)


med.cpl.rm<-lapply(1:numImp,function(x) complete(med.item.imp.rm,action = x))
med.cpl.rm.out<-do.call(abind, c(med.cpl.rm, along = 3))
med.cpl.rm.out.mean<-apply(med.cpl.rm.out, c(1,2), mean)
med.cpl.rm.out.round <- round(med.cpl.rm.out.mean)

save(med.cpl.rm.out.round,med.item.imp.all,med.item.imp.rm,med.item,med.item.rm,med.vis,file = './result/med_impute.RData')

