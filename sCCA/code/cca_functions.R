
cv_score <- function(uhat_tr,vhat_tr,X_ts,Y_ts){
  score <- cor(X_ts %*% uhat_tr , Y_ts %*% vhat_tr )
} 

ccaDW <- function(X,Y,pen_x,pen_y,rank){
  mode<-PMA::CCA(x=X, z=Y, typex=c("standard"),typez=c("standard"), penaltyx=pen_x, penaltyz=pen_y, K=rank, niter=100, trace=FALSE)
}

ccaDWcv <- function(X,Y,pen_x,pen_y,rank,foldid) {  
  # run cca on the training set
  cca_train <- mclapply(foldid, function(id)  ccaDW( X[id,] , Y[id,] , pen_x , pen_y , rank))
  
  # apply the loadings onto testing data and calculate the CV score on each fold
  cv_results <- mclapply(seq_along(foldid), function(i) cv_score(cca_train[[i]]$u, cca_train[[i]]$v, X[ -foldid[[i]] , ] , Y[ -foldid[[i]] , ]))
  
  # take the mean CV
  cv_mean <- mean(unlist(cv_results))
  cv_sd <- sd(unlist(cv_results))
  cv_all <- unlist(cv_results)
  
  # take ave of the correlation estimated
  fold_cor <- mclapply(seq_along(foldid), function(i) cca_train[[i]]$cors)
  cor_mean <- mean(unlist(fold_cor))
  
  
  # output
  out<-list(CV_SCORE=cv_mean,CV_SD = cv_sd, CV_ALL = cv_all, COR_MEAN=cor_mean, PEN_X = pen_x, PEN_Y = pen_y,RANK=rank)
}

ccaDWgs <-function(X,Y,pen_xseq,pen_yseq, rank, foldid){
  #loop through x penalty sequence
  gs_x<-function(X,Y,pen_xseq,pen_yseq, rank,foldid){
    gs_x <- mclapply(pen_xseq,function(x) ccaDWcv(X,Y,x,pen_yseq,rank,foldid))
  }
  #loop through y penalty sequence
  gs_y <- mclapply(pen_yseq,function(y) gs_x(X,Y,pen_xseq, y,rank,foldid))
  
  #summarize the result
  gs2 <- unlist(gs_y,recursive = F)
  gs.out <- t(sapply(gs2, function(x) unlist(x)))
  
  #find the parameters with the best CV score
  cv.max <- max(gs.out[,'CV_SCORE'])
  cv.max.loc <- which.max(gs.out[,'CV_SCORE'])
  
  best.para = gs.out[cv.max.loc,]
  final.cca <- ccaDW(X,Y,as.numeric(best.para['PEN_X']),as.numeric(best.para['PEN_Y']),rank)
  out <- list(GRIDSEARCH=gs.out, BEST_PARA=best.para , FINAL_CCA = final.cca)
}

# permutation test in the test data
ccaDWpm <- function(Xtest,Ytest,uhat,vhat,nperm){
  Xtest.perm <- rlply(nperm,Xtest[sample(nrow(Xtest)),])
  
  Xtest.perm.cor <- mclapply(Xtest.perm, function(xperm) cv_score(uhat, vhat, xperm, Ytest) )
  
  Xtest.real.cor <- cv_score(uhat,vhat,Xtest,Ytest)
  
  p.val <- length(which(abs(unlist(Xtest.perm.cor)) >= abs(as.numeric(unlist(Xtest.real.cor)))))/nperm
  out <- list(PERMCOR = Xtest.perm.cor, REALCOR = Xtest.real.cor  , PVAL = p.val )
}

# permutation test within the training data
ccaDWpermtr <- function(Xtrain,Ytrain,nperm,x_pen,y_pen,rank){
  real.cor <- ccaDW(Xtrain,Ytrain,x_pen,y_pen,rank)
  
  Xtrain.perm <- rlply(nperm,Xtrain[sample(nrow(Xtrain)),])
  Xtrain.perm <- mclapply(Xtrain.perm, function(xperm) ccaDW(xperm,Ytrain,x_pen,y_pen,rank))
  
  perm.cor <- sapply(seq_along(Xtrain.perm), function(i) Xtrain.perm[[i]]$cors)
  p.val <- length(which(perm.cor >= real.cor$cors))/nperm
  out <- list(REALCOR = real.cor$cors, PVAL = p.val, PEN_X = x_pen, PEN_Y = y_pen)
}

ccaDWpermgs <-function(X,Y,pen_xseq,pen_yseq,rank, nperm){
  #loop through x penalty sequence
  gs_x<-function(X,Y,pen_xseq,pen_yseq, rank,nperm){
    gs_x <- mclapply(pen_xseq,function(x) ccaDWpermtr(X,Y,nperm,x,pen_yseq,rank))
  }
  #loop through y penalty sequence
  gs_y <- mclapply(pen_yseq,function(y) gs_x(X,Y,pen_xseq, y,rank,nperm))
  
  #summarize the result
  gs2 <- unlist(gs_y,recursive = F)
  gs.out <- t(sapply(gs2, function(x) unlist(x)))
  
  #find the parameters with the smallest perm p-val
  p.min <- min(gs.out[,'PVAL'])
  
  p.min.loc <- which.min(gs.out[,'PVAL'])
  
  best.para = gs.out[p.min.loc,]
  final.cca <- ccaDW(X,Y,as.numeric(best.para['PEN_X']),as.numeric(best.para['PEN_Y']),rank)
  out <- list(GRIDSEARCH=gs.out, BEST_PARA=best.para , FINAL_CCA = final.cca)
}

ccaDWfoldgs <-function(Xlist,Ylist,pen_xseq,pen_yseq){
  #loop through x penalty sequence
  gs_x<-function(X,Y,pen_xseq,pen_yseq){
    gs_x <- mclapply(pen_xseq,function(pen_x) ccaDWfold(X,Y,pen_x,pen_yseq))
  }
  #loop through y penalty sequence
  gs_y <- mclapply(pen_yseq,function(pen_y) gs_x(Xlist,Ylist,pen_xseq,pen_y))
  
  #summarize the result
  gs2 <- unlist(gs_y,recursive = F)
  gs.out <- t(sapply(gs2, function(x) unlist(x)))
  
  #find the parameters with the largest correlations
  best.para <- gs.out[which.max(gs.out[,'COR_MEAN']),]
  best.cor <- as.numeric(best.para[1])
  best.penx <- as.numeric(best.para[2])
  best.peny <- as.numeric(best.para[3])
  out <-list(GS = gs.out,COR= best.cor, PENX = best.penx, PENY = best.peny)
}

ccaDWfold <- function(Xlist,Ylist,pen_x,pen_y){
result.cca<-mclapply(seq_along(Xlist),function(i) ccaDW(Xlist[[i]],Ylist[[i]],pen_x,pen_y,2))
result.cor <- sapply(seq_along(Xlist), function(i) result.cca[[i]]$cors[2])
result.cor.mean <- mean(result.cor)

#sapply(seq_along(Xlist),function(i) {Xlist[[i]] %*% })

out <-list(COR_MEAN = result.cor.mean, PEN_X= pen_x, PEN_Y = pen_y)
}