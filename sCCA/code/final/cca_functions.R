load_power_net <- function(sample_qa){
  n_sample <- dim(sample_qa)[1]
  sample_net<-array(NA, c(264, 264, n_sample))
  for (i in 1:n_sample){
    scanid <- sample_qa$scanid[i]
    netpath<- paste("./studies/pnc/n1601_dataFreeze2016/neuroimaging/rest/restNetwork_264PowerPNC/264PowerPNCNetworks/",scanid,"_264PowerPNC_network.txt",sep="")
    sample_net[,,i] <- as.matrix(read.table(netpath))
    print(paste(i,"."," copying ",scanid,"_","Power",sep=""))
  }
  
  #make feature table of the matrix
  print("making feature table")
  net_ft <-t(apply(sample_net,c(3),function(x) x[upper.tri(x, diag = F)]))
  
  print("making inclusion net mat")
  
  inc_net <- net_ft[,pwr.3k.train.idx]
}


load_reg_power_net <- function(sample_qa){
  #load net
  inc_net <- load_power_net(sample_qa)
  
  # Compile covariates
  print("Compile covariates")
  sample_qa$sex <- as.factor(sample_qa$sex)
  sample_qa$race2 <- as.factor(sample_qa$race2)
  
  # Regress
  print("Regress out covariates")
  inc_net[which(is.na(inc_net))] <- 0
  net.rgr <- apply(inc_net, 2, function(x) residuals.glm(glm(x ~ ageAtScan1 + 
                                                               sex + race2 + restRelMeanRMSMotion, data = sample_qa), type = "response"))
}

cv_score <- function(uhat_tr,vhat_tr,X_ts,Y_ts){
  score <- CCA(X_ts %*% uhat_tr , Y_ts %*% vhat_tr )
} 

ccaDW <- function(X,Y,pen_x,pen_y,rank){
  mode<-PMA::CCA(x=X, z=Y, typex=c("standard"),typez=c("standard"), penaltyx=pen_x, penaltyz=pen_y, K=rank, niter=100, trace=FALSE)
}

ccaDWpermorder <- function(X,Y,pen_x,pen_y,rank,cca_org){
  perm.mode<-PMA::CCA(x=X, z=Y, typex=c("standard"),typez=c("standard"), penaltyx=pen_x, penaltyz=pen_y, K=rank, niter=20, trace=FALSE)
  perm.mode.reorder<-reorderCCA(perm.mode,cca_org,rank)
}

ccaDWpermorg <- function(X,Y,pen_x,pen_y,rank){
  perm.mode<-PMA::CCA(x=X, z=Y, typex=c("standard"),typez=c("standard"), penaltyx=pen_x, penaltyz=pen_y, K=rank, niter=20, trace=FALSE)
}

ccaDWpermrank <- function(X,Y,pen_x,pen_y,rank){
  perm.mode<-PMA::CCA(x=X, z=Y, typex=c("standard"),typez=c("standard"), penaltyx=pen_x, penaltyz=pen_y, K=rank, niter=20, trace=FALSE)
  cor.order <- order(-perm.mode$cors)
  perm.mode$cors <- perm.mode$cors[cor.order]
  perm.mode$u <- perm.mode$u[,cor.order]
  perm.mode$v <- perm.mode$v[,cor.order]
  perm.mode
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
  result.cca<-mclapply(seq_along(Xlist),function(i) ccaDW(Xlist[[i]],Ylist[[i]],pen_x,pen_y,1))
  result.cor <- sapply(seq_along(Xlist), function(i) result.cca[[i]]$cors)
  result.cor.mean <- mean(result.cor)
  out <-list(COR_MEAN = result.cor.mean, PEN_X= pen_x, PEN_Y = pen_y)
}

perm.plot <-function(perm_file,cor_file,pval,dim){
  perm_file <- data.frame(cor_perm = perm.cor.df[,dim])
  p <- ggplot(perm_file,aes(cor_perm))+
    geom_histogram(binwidth = 0.005, fill = "blue", alpha = 0.5)  +
    geom_vline(xintercept = cor_file$cor[dim], colour = "red", linetype = "longdash") +
    labs(x = "Correlations") +
    annotate("text", x = median(perm_file$cor_perm,na.rm = T), y = c(10,5),label = c("Permuted Data","(1000 times)"),size =10,colour = "black" ) +
    annotate("text",x = -Inf, y = Inf, hjust = -0.1,vjust = 1,label = paste("Cor=",round(cor_file$cor[dim],2),", p<",round(pval[dim],3)), size = 10, colour = "red" ) +
    theme_classic(base_size = 35) + 
    scale_y_continuous(expand = c(0, 0))
  p
}

reorderCCA <- function(res,org,k){
  u <- org$u
  v <- org$v
  
  res.match<-apply(abs(cor(v,res$v[,1:k])),1,which.max)
  res.cor <- apply(abs(cor(v,res$v[,1:k])),1,max)
  res.match.count <- length(unique(res.match))
  
  u.od <- res$u[,res.match]
  v.od <- res$v[,res.match]
  cors.final <- res$cors[res.match]
  
  org.sign <- sign(colMeans(sign(v)))
  res.sign <- sign(colMeans(sign(v.od)))
  
  sign.prod <- org.sign *res.sign
  
  u.final <- t(t(u.od) * sign.prod)
  v.final <- t(t(v.od) * sign.prod)
  
  res.one.reorder <- list(u= u.final ,v= v.final, cors = cors.final, pos = res.match, dimcor = res.cor)
  
  if (res.match.count < dim(u)[2]  ) {
    u.na <- array(NA,dim= c(dim(u)[1],dim(u)[2]))
    v.na <- array(NA,dim= c(dim(v)[1],dim(v)[2]))
    cors.na <- rep(NA,dim(u)[2])
    match.na <- cors.na
    dimcor.na <- cors.na
    res.one.reorder <- list(u= u.na,v= v.na, cors = cors.na, pos = res.match, dimcor = dimcor.na) 
  }
  out <- res.one.reorder
}

bootstats <- function(bootdata){
  boot.stats <- data.frame(mean = rowMeans(bootdata,na.rm = T), se= rowSds(bootdata,na.rm = T) )
  boot.stats$me = boot.stats$se * qt(.95, df=dim(bootdata)[2]-1) 
  boot.stats$up <- boot.stats$mean + boot.stats$me
  boot.stats$lw <- boot.stats$mean - boot.stats$me
  boot.stats
}

bootstats2 <- function(org,bootdata, uplim,lwlim){
  boot.stats <- sapply(seq_along(1:length(org)), function(i) org[i] - quantile(bootdata[i,] - org[i],c(uplim,lwlim),na.rm =T))
  boot.stats <- as.data.frame(t(boot.stats))
  colnames(boot.stats)<-c("low","high")
  boot.stats$ci <- boot.stats$high - boot.stats$low
  boot.stats$load <- org
  boot.stats$fea <- 1:length(org)
  boot.stats
  }

bootplot <- function(org, boot){
  btst <- bootstats2(org,boot,0.975,0.025)
  
  p <- ggplot(btst,aes(fea,load))+
    geom_point(aes(colour = low * high > 0)) +
    geom_errorbar(aes(ymax = high, ymin = low),  width=0.25) +
    theme_classic()
  
  fea <- which(btst$low * btst$high >0)
  load <- org[fea]
  fea <- fea[order(-abs(load))]
  load <- load[order(-abs(load))]
  out <- list(plot = p, fea = fea, load = load)
}

## write an alternative bootplot function for the u's so it doesn't involve the absolute values

bootplot_u <- function(org, boot){
  btst <- bootstats2(org,boot,0.995,0.005)
  
  p <- ggplot(btst,aes(fea,load))+
    geom_point(aes(colour = low * high > 0)) +
    geom_errorbar(aes(ymax = high, ymin = low),  width=0.25) +
    theme_classic()
  
  fea <- which(btst$low * btst$high >0)
  load <- org[fea]
  #fea <- fea[order(-abs(load))]
  #load <- load[order(-abs(load))]
  load_nm <- org / btst$ci
  load_nm_fea <- load_nm[fea]
  out <- list(plot = p, fea = fea, load = load, ci = btst$ci[fea], high = btst$high[fea], low = btst$low[fea],load_nm = load_nm, load_nm_fea = load_nm_fea)
}

med_vis <- function(mddim,title){
  mddim$feaitem <- item$Question[mddim$fea]
  mddim$code <- item$clinicalcode[mddim$fea]
  mddim$color <- item$color[mddim$fea]
  mddim$cate <- item$cate[mddim$fea]
  #mddim$cate <- gsub("[[:digit:]]","",mddim$code)
  mddim.vis <- as.data.frame(within(mddim,rm("plot")))
  cats <- unique(mddim.vis$cate)
  cols <- unique(mddim.vis$color)
  ggplot(mddim.vis,aes(rep(1,length(mddim$load)),abs(load), label = feaitem)) +
    #geom_point(aes(size= 0,alpha = 0, color = color )) +
    geom_point(aes(size= abs(load*10),alpha =abs(load) )) +
    scale_alpha(guide = 'none') +
    scale_size(guide = "none") +
    coord_cartesian(xlim = c(1, 20)) +
    geom_label_repel(data = mddim.vis, fill = mddim.vis$color,
                     colour = "white", force = 3, nudge_x=10, show.legend = TRUE) +
    #scale_colour_manual(name='', values=c('Important line'='grey', 'Point values'='red')) +
    theme_classic(base_size = 15) +
    theme(axis.ticks = element_blank(), axis.text.x = element_blank()) +
    theme(legend.position="top") +
    theme(plot.title = element_text(hjust=0.5)) +
    scale_fill_discrete(name="") +
    labs(title = title ,x='',y= "CCA Loadings")
}

brain_vis <- function(braindim,title,sign,idx, parcels){
  ccamat <- array(0,c(264,264)) #create a blank matrix
  ccamat.up<-ccamat[upper.tri(ccamat,diag=F)] 
  ccamat.up[idx[braindim$fea]] <- braindim$load
  
  #ccamat.up[idx] <- braindim$load_nm
  
  ccamat[upper.tri(ccamat,diag=F)]  <- ccamat.up
  ccamat <- sna::symmetrize(ccamat,rule = 'upper')
  ccamat <- sign * ccamat
  ccamat <- ccamat[parcels$ROI,parcels$ROI]
  ccamat.lim <- max(abs(ccamat),na.rm=T)
  p <- levelplot(ccamat,par.settings = BuRdTheme(), at = seq(-ccamat.lim,ccamat.lim,length.out = 100), xlab="",ylab = "", main=title)
  
  our <- list(plot = p, mat = ccamat)
}

load_by_ft_plot <- function(load_mat,ave_mat,title){
  #emp_mat <- array(0, c(264,264))
  #ave_ft<-colMeans(pwr_train_net_ft,na.rm = T)
  #emp_mat[upper.tri(emp_mat,diag = F)] <- ave_ft
  #ave_ft_mat <- sna::symmetrize(emp_mat,rule = "upper") 
  #ave_ft_mat <- ave_ft_mat[parcelsTR$ROI,parcelsTR$ROI]
  
  #ave_ft_lim <- max(abs(ave_ft_mat),na.rm=T)
  #levelplot(ave_ft_mat,par.settings = BuRdTheme(), at = seq(-ave_ft_lim,ave_ft_lim,length.out = 100))
  
  loadbyft <- sign(ave_mat) * load_mat
  loadbyft_lim <- max(abs(loadbyft),na.rm=T)
  
  p <- levelplot(loadbyft,par.settings = BuRdTheme(), at = seq(-loadbyft_lim,loadbyft_lim,length.out = 100),xlab="",ylab = "",main=title )
  
  out <- list(plot = p, mat = loadbyft)
}

load_norm <- function(load_mat,ft_mat,ft_rgr_mat,title){
  #emp_mat <- array(0, c(264,264))
  #ave_ft<-colMeans(pwr_train_net_ft,na.rm = T)
  #emp_mat[upper.tri(emp_mat,diag = F)] <- ave_ft
  #ave_ft_mat <- sna::symmetrize(emp_mat,rule = "upper") 
  #ave_ft_mat <- ave_ft_mat[parcelsTR$ROI,parcelsTR$ROI]
  
  #ave_ft_lim <- max(abs(ave_ft_mat),na.rm=T)
  #levelplot(ave_ft_mat,par.settings = BuRdTheme(), at = seq(-ave_ft_lim,ave_ft_lim,length.out = 100))
  
  #u_by_x <- load_mat$mat * ft_rgr_mat$ave_mat
  
  #ft_hat <- ft_mat$ave_mat - ft_rgr_mat$ave_mat
  
  #prod.sign <- sign(u_by_x) * sign(ft_hat)
  
  #load.final <- abs(load_mat$mat) * prod.sign
  #load.final <- abs(load_mat$mat) * sign(u_by_x)
  #load.final <- u_by_x
  #load.final.lim <- max(abs(load.final),na.rm=T)
  
  #p <- levelplot(load.final,par.settings = BuRdTheme(), at = seq(-load.final.lim,load.final.lim,length.out = 100),xlab="",ylab = "",main=title )
  
  u.by.x <- load_mat$mat * ft_rgr_mat$ave_mat
  delta.percent <- u.by.x / (ft_mat$ave_mat - ft_rgr_mat$ave_mat)
  u.by.x.sign <- u.by.x * sign(ft_mat$ave_mat - ft_rgr_mat$ave_mat)
  u.by.sign <- load_mat$mat * sign(ft_rgr_mat$ave_mat)
  lim = max(abs(u.by.sign))
  p<-levelplot(u.by.sign, par.settings = BuRdTheme(), at = seq(-lim,lim,length.out = 100))
  load.final <- u.by.sign
  
  out <- list(plot = p, mat = load.final, ux = u.by.x, uxs = u.by.x.sign)
}



mod_calc_within <- function (bi_brain_mat, mod_assign){
 within<-sapply(1:length(unique(mod_assign)),function(x) mean(bi_brain_mat[mod_assign == unique(mod_assign)[x],mod_assign == unique(mod_assign)[x]]))
 
}


mod_pair<-function(bi_brain_mat, mod_assign, seed){
  
  seed_mod <- unique(mod_assign)[seed]
  pair_seq <- (1+seed):length(unique(mod_assign))
  sapply(pair_seq, function(x) mean(bi_brain_mat[mod_assign == seed_mod, mod_assign == unique(mod_assign)[x] ], na.rm = T))
}

mod_calc_between <- function (bi_brain_mat, mod_assign){
between<-sapply(1:(length(unique(mod_assign))-1),function(x) mod_pair(bi_brain_mat,mod_assign,x))
between <- unlist(between)
}

mod_pair_name<-function(mod_name, seed){
  seed_mod <- unique(mod_name)[seed]
  pair_seq <- (1+seed):length(mod_name)
  sapply(pair_seq, function(x) paste(seed_mod,'-',mod_name[x]))
}

mod_calc_between_name <- function (mod_name){
  between<-sapply(1:(length(mod_name)-1),function(x) mod_pair_name(mod_name,x))
  between <- unlist(between)
}

mod_rich_within<-function(brbft) {
  #binarize the brain by ft matrix
  brbft[is.na(brbft)] <-0
  brbft[which(brbft != 0)] <- 1
  
  
  #permute the membership of modules
  nperm = 5000
  mod_mem_perm <- rlply(nperm,sample(parcelsTR$Community))
  
  #calculate mod enrichment
  modrich <- mod_calc_within(brbft,parcelsTR$Community)
  
  #calc perm mod enrich
  mod_perm<-sapply(mod_mem_perm,function(x) mod_calc_within(brbft, x))
  
  #calc p val
  mod_pval <- sapply(1:dim(mod_perm)[1], function(x) length(which(mod_perm[x,]>=modrich[x])) / nperm )
  sig_mod_idx <- which(mod_pval < 0.05)
  sig_mod <-unique(parcelsTR$System)[sig_mod_idx]
  sig_mod_idx <- unique(parcelsTR$Community)[sig_mod_idx]
  
  out <- list(RICH = modrich, PVAL = mod_pval, MOD = sig_mod, MODid = sig_mod_idx)
}

mod_rich_between<-function(brbft) {
  #calc between mod mean loading
  modmean <- mod_calc_between(brbft,parcelsTR$Community)
  
  #binarize the brain by ft matrix
  brbft[is.na(brbft)] <-0
  brbft[which(brbft != 0)] <- 1
  
  
  #permute the membership of modules
  nperm = 1000
  mod_mem_perm <- rlply(nperm,sample(parcelsTR$Community))
  
  #calculate mod enrichment
  modrich <- mod_calc_between(brbft,parcelsTR$Community)
  
  #calc perm mod enrich
  mod_perm<-sapply(mod_mem_perm,function(x) mod_calc_between(brbft, x))
  
  #calc p val
  mod_pval <- sapply(1:dim(mod_perm)[1], function(x) length(which(mod_perm[x,]>=modrich[x])) / nperm )
  mod_pval_fdr <- p.adjust(mod_pval,method = "fdr")

  sig_mod_idx <- which(mod_pval_fdr <0.05)
  sig_mod <-between_name[sig_mod_idx]

  out <- list(MEAN = modmean, RICH = modrich, PVAL = mod_pval,PVAL_fdr = mod_pval_fdr, MOD = sig_mod, MODid = sig_mod_idx)
}


within_mod_plot <- function(netmat,mod){
  brainmat <- netmat
  brainlim = max(abs(brainmat),na.rm = T)
  mod <- unique(parcelsTR$Community)[mod]
  title <- toString(unique(parcelsTR$System[which(parcelsTR$Community==mod)]))
  levelplot(brainmat[parcelsTR$Community==mod,parcelsTR$Community==mod],par.settings = BuRdTheme(), at = seq(-brainlim,brainlim,length.out=100), xlab="", ylab="",main = title )
}

within_mod_calc <- function(netmat,mod){
  title <- modname[mod]
  brainmat <- netmat
  brainlim = max(abs(brainmat),na.rm = T)
  mod <- unique(parcelsTR$Community)[mod]
  list(mod = title, mean = mean(brainmat[parcelsTR$Community==mod,parcelsTR$Community==mod],na.rm=T))
}

listplot<-function(list_of_plot,subname) {
  for (i in seq_along(list_of_plot)) {
    plotname <- paste('./projects/xiaNetworkCca/sCCA/aim1/figure/201701/',subname,'_',i,'.pdf',sep="")
    pdf(file = plotname)
    print(list_of_plot[[i]])
    dev.off()
  }
}


between_mod_plot<-function(betweenmod,withinmod, title) {
  mod_mat <- array(0,c(12,12))
  mod_mat[lower.tri(mod_mat,diag = F)][betweenmod$MODid] <- 1
  withinmod_id<-rapply(withinmod[1,], function(x) which(modname == x))
  diag(mod_mat)[withinmod_id] <- 1
  mod_mat <- sna::symmetrize(mod_mat,rule = "lower")
  colnames(mod_mat) <- modname
  rownames(mod_mat) <- modname
  p<-levelplot(mod_mat, par.settings = BuRdTheme(), at = seq(-1,1,length.out = 10), xlab="",ylab = "",scales=list(x=list(rot=90)),colorkey = FALSE, main = title)
  print(p)
}

between_mod_load_plot<-function(betweenmod,withinmod,title,i) {
  sysdim <- read.csv("~/Google Drive/TDSlab/CEDRIC/PNC_CCA/Figure Resources/systemdim.csv")
  keycolor = c(rev(brewer.pal(9,"Blues")[3:9]),"#FFFFFF",brewer.pal(9,"OrRd")[3:9])
  md_dim_color = c("#8E24AA","#1E88E5","#FF6F00","#D32F2F","#D32F2F")
  mod_mat <- array(0,c(12,12))
  mod_mat[lower.tri(mod_mat,diag = F)][betweenmod$MODid] <- betweenmod$MEAN[betweenmod$MODid]
  withinmod_id<-rapply(withinmod[1,], function(x) which(modname == x))
  diag(mod_mat)[withinmod_id] <- simplify2array(withinmod[2,])
  mod_mat <- sna::symmetrize(mod_mat,rule = "lower")
  colnames(mod_mat) <- modname
  rownames(mod_mat) <- modname
  absmax <- max(abs(mod_mat),na.rm = T)
  mod_mat <- mod_mat/absmax
  #mod_mat <- mod_mat/0.002607361
  modlim <- max(abs(mod_mat),na.rm = T)
  #modlim <- 1
  plotrange <- c(5:8,10:11)
  #col.regions = keycolor
  #par.settings = BuRdTheme()
  lattice.options(axis.padding=list(factor=0.5))
  
  p<-levelplot(mod_mat[plotrange,plotrange], col.regions = keycolor , colorkey = TRUE,
               at = seq(-modlim,modlim,length.out = 16), 
               xlab="",ylab = "",main = list(title, cex = 2, col = md_dim_color[i], just = 0.425),
               scales=list(x=list(rot=90,tck = 0, cex = 1.5, font=2,col = (as.character(sysdim$Color[plotrange]))),
                           y=list(tck = 0, cex = 1.5, font=2,col = (as.character(sysdim$Color[plotrange])))))
  
  out <- list(MAT = mod_mat, plot = p, max = absmax)
}

bt_mod_diff <- function(dim1, dim2, title){
  dimdiff <- dim1$MAT - dim2$MAT
  difflim <- max(abs(dimdiff),na.rm = T)
  p<-levelplot(dimdiff, par.settings = BuRdTheme(), at = seq(-difflim,difflim,length.out = 100), xlab="",ylab = "",scales=list(x=list(rot=90)), main = title)
  print(p)
}

mask_mat<-function(sample_net_ft_rg,idx) {
  sample_net_ave <- colMeans(sample_net_ft_rg,na.rm = TRUE)
  ave_mat <- array(0,c(264,264))
  ave_mat[upper.tri(ave_mat,diag=F)][idx] <- sample_net_ave
  ave_mat <- sna::symmetrize(ave_mat,rule = "upper")
  ave_mat <- ave_mat[parcelsTR$ROI,parcelsTR$ROI]
  mat_lim <- max(abs(ave_mat))
  p<-levelplot(ave_mat, par.settings = BuRdTheme(),at = seq(-mat_lim,mat_lim,length.out = 100))
  out <- list(ave_mat = ave_mat, plot = p)
}

mask_mat_mod<-function(sample_net_ft_rg,mod) {
  sample_net_ave <- colMeans(sample_net_ft_rg)
  ave_mat <- array(0,c(264,264))
  ave_mat[upper.tri(ave_mat,diag=F)][pwr.3k.train.idx] <- sample_net_ave
  ave_mat <- sna::symmetrize(ave_mat,rule = "upper")
  ave_mat <- ave_mat[parcelsTR$ROI,parcelsTR$ROI]
  ave_mat <- ave_mat[parcelsTR$Community==mod,parcelsTR$Community==mod]
  mat_lim <- max(abs(ave_mat))
  p<-levelplot(ave_mat, par.settings = BuRdTheme(),at = seq(-mat_lim,mat_lim,length.out = 100))
  out <- list(ave_mat = ave_mat, plot = p)
}

mask_mat_ind<-function(ftmat) {
  mat <- array(0,c(264,264))
  mat[upper.tri(mat,diag=F)][pwr.3k.train.idx] <- ftmat
  mat <- sna::symmetrize(mat,rule = "upper")
  mat <- mat[parcelsTR$ROI,parcelsTR$ROI]
  mat_lim <- max(abs(mat))
  p<-levelplot(mat, par.settings = BuRdTheme(),at = seq(-mat_lim,mat_lim,length.out = 100))
  out <- list(mat = mat, plot = p)
}
