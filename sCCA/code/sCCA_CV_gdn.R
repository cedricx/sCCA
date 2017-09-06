## Load necessary libraries
require('PMA')
require('Matrix')
require('parallel')
require('emdbook')
require('caret')
require('R.matlab')
require('MASS')
require('permute')
require('matrixStats')
require('scales')

## Load Data
setwd('/data/joy/BBL/projects/xiaNetworkCca/sCCA/aim1/')
load("/data/joy/BBL/projects/xiaNetworkCca/sCCA/aim1/result/201701/gdn_regr_data.RData")
load("/data/joy/BBL/projects/xiaNetworkCca/sCCA/aim1/result/201701/gdn_med_rgr_data.RData")
net.data <- unname(gdn.rgr.train)
med.data <- unname(gdn_train_med_rgr)
data <- list(brain = net.data, behavior = med.data)

## create sub-CV splits
load("/data/joy/BBL/projects/xiaNetworkCca/sCCA/aim1/result/201701/gdn_train_qa.RData")
subjid <- gdn_train_qa

source('./../code/cca_functions.R')

# create 3 fold CV sets 10 times in the sub-training set

trainid <- createDataPartition(subjid$overall_psychopathology_4factor, p = 0.667, list =T,times=10)

brain_train <- mclapply(trainid, function(id) data$brain[id,])
brain_test <- mclapply(trainid, function(id) data$brain[-id,])

behavior_train <- mclapply(trainid, function(id) data$behavior[id,])
behavior_test <- mclapply(trainid, function(id) data$behavior[-id,])

# select best parameters
x_pen <- seq(0.1,1,length.out=10)
y_pen <- seq(0.1,1,length.out=10)

p3Km111.gs<-ccaDWfoldgs(brain_train,behavior_train,x_pen,y_pen)
p3Km111.cca<-mclapply(seq_along(trainid),function(i) ccaDW(brain_train[[i]],behavior_train[[i]],p3Km111.gs$PENX,p3Km111.gs$PENY,111))
p3Km111.cca<-mclapply(seq_along(trainid),function(i) ccaDW(brain_train[[i]],behavior_train[[i]],0.8,0.4,111))


p3Km111.ve <- sapply(seq_along(trainid), function(i) p3Km111.cca[[i]]$d^2/sum(p3Km111.cca[[i]]$d^2))

# plot ve
# variance explained plots
varE <- data.frame( modenum = as.factor(1:111),ve_mean = rowMeans(p3Km111.ve),ve_se = rowSds(p3Km111.ve)/sqrt(dim(p3Km111.ve)[2]),vediff = c(abs(diff(rowMeans(p3Km111.ve))),0))
varlimits <- aes(ymax = varE$ve_mean + varE$ve_se, ymin=varE$ve_mean - varE$ve_se)
p.varE<-ggplot(varE,aes(modenum,ve_mean))+ 
  geom_bar(stat = 'identity',aes(fill = ve_mean >= ve_mean[3] )) +
  geom_errorbar(varlimits,  width=0.25) +
  geom_hline(yintercept = 1/111,linetype="dashed") +
  scale_x_discrete(name ="Mode", limits=c(1:20),breaks =  seq(1,20,2)) +
  scale_y_continuous(expand = c(0,0),limits=c(0,0.075),labels = percent,name = "Variance Explained", breaks=seq(0,0.075,length=5)) +
  #guides(fill = guide_legend(title = "diff >= 4SD")) +
  theme_classic(base_size = 30) +
  theme(legend.position = "none") 
p.varE

k=3;

# correlation
p3Km111.cor <- sapply(seq_along(trainid), function(i) p3Km111.cca[[i]]$cors)
p3Km111.cor.mean <- rowMeans(p3Km111.cor[1:k,])

# cca corelation plots
cor <- data.frame(modenum = as.factor(1:111),cor_mean = rowMeans(p3Km111.cor), cor_se = rowSds(p3Km111.cor)/sqrt(dim(p3Km111.cor)[2]) )
corlimits <- aes(ymax = cor$cor_mean + cor$cor_se, ymin = cor$cor_mean - cor$cor_se)
p.cor<-ggplot(cor,aes(modenum,cor_mean,label=round(cor_mean,2)))+ 
  geom_bar(width=.5,stat = 'identity', fill = '#80CBC4') +
  geom_errorbar(corlimits,  width=0.25) +
  geom_text(size = 5, position = position_dodge(width=0.9), vjust= -1,color='grey')+
  scale_x_discrete(name ="Mode", limits=c(1:k)) +
  ylab('CCA Correlation') +
  coord_cartesian(ylim=c(0.2,0.8)) +
  theme_classic(base_size = 30) +
  theme(legend.position="none")
p.cor

# calculate the CV
p3Km111.cv <- sapply(seq_along(trainid), function(i) { behav_hat<-(brain_test[[i]] %*% p3Km111.cca[[i]]$u[,1:k]) %*% ginv(p3Km111.cca[[i]]$v[,1:k]); out<-cancor(behav_hat,behavior_test[[i]]); out$cor })
# cross-validation plots
k=3
cv <- data.frame(modenum = as.factor(1:k),cv_mean = rowMeans(p3Km111.cv), cv_se = rowSds(p3Km111.cv)/sqrt(dim(p3Km111.cv)[2]) )
cvlimits <- aes(ymax = cv$cv_mean + cv$cv_se, ymin = cv$cv_mean - cv$cv_se)
p<-ggplot(cv,aes(modenum,cv_mean,label=round(cv_mean,2)))+ 
  geom_bar(width=.5,stat = 'identity', fill = '#80CBC4') +
  geom_errorbar(cvlimits,  width=0.25) +
  geom_text(size = 5, position = position_dodge(width=0.9), vjust= -0.8,color='grey')+
  scale_x_discrete(name ="Mode", limits=c(1:3)) +
  ylab('Cross-Validation Accuracy') +
  coord_cartesian(ylim=c(0.2,0.9)) +
  theme_classic(base_size = 24) +
  theme(legend.position="none")
p

# build a null model on all data
behavior.perm <- rlply(1000,data$behavior[sample(nrow(data$behavior)),])
p3Km111.perm.cor<-sapply(behavior.perm, function(y_perm){ out<-ccaDW(data$brain,y_perm,p3Km111.gs$PENX,p3Km111.gs$PENY,k); out$cor} )
p3Km111.perm.cor<-sapply(behavior.perm, function(y_perm){ out<-ccaDW(data$brain,y_perm,0.7,0.5,k); out$cor} )

## bootstrap
bootid<-createResample(subjid$overall_psychopathology_4factor, list = T, times = 10)
brain_boot <- mclapply(bootid, function(id) data$brain[id,])
behavior_boot <- mclapply(bootid, function(id) data$behavior[id,])

behavior.perm <- rlply(10,data$behavior[sample(nrow(data$behavior)),])

p3Km111.org<-ccaDW(data$brain,data$behavior,0.8,0.5,3)

p3Km111.boot<- mclapply(seq_along(bootid),function(i) ccaDW(brain_boot[[i]],behavior_boot[[i]],0.8,0.5,10))
p3Km111.perm<- mclapply(seq_along(behavior.perm),function(i) ccaDW(data$brain,behavior.perm[[i]],0.8,0.5,10))

p3Km111.boot.ro<- mclapply(seq_along(bootid),function(i) reorderCCA(p3Km111.boot[[i]],p3Km111.org) )

p3Km111.perm.ro<- mclapply(seq_along(bootid),function(i) reorderCCA(p3Km111.perm[[i]],p3Km111.org) )

