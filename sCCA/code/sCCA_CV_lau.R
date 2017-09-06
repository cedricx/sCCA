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
load("/data/joy/BBL/projects/xiaNetworkCca/sCCA/aim1/result/201701/lau_regr_data.RData")
load("/data/joy/BBL/projects/xiaNetworkCca/sCCA/aim1/result/201701/pwr_med_rgr_data.RData")
net.data <- unname(lau.rgr.train)
med.data <- unname(pwr_train_med_rgr)
data <- list(brain = net.data, behavior = med.data)

## create sub-CV splits
load("/data/joy/BBL/projects/xiaNetworkCca/sCCA/aim1/result/201701/lau_train_qa.RData")
subjid <- lau_train_qa

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
p3Km111.cca<-mclapply(seq_along(trainid),function(i) ccaDW(brain_train[[i]],behavior_train[[i]],0.8,0.4,111))

p3Km111.ve <- sapply(seq_along(trainid), function(i) p3Km111.cca[[i]]$d^2/sum(p3Km111.cca[[i]]$d^2))

# plot ve
# variance explained plots
varE <- data.frame( modenum = as.factor(1:111),ve_mean = rowMeans(p3Km111.ve),ve_se = rowSds(p3Km111.ve)/sqrt(dim(p3Km111.ve)[2]),vediff = c(abs(diff(rowMeans(p3Km111.ve))),0))
varlimits <- aes(ymax = varE$ve_mean + varE$ve_se, ymin=varE$ve_mean - varE$ve_se)
p.var<-ggplot(varE,aes(modenum,ve_mean))+ 
  geom_bar(stat = 'identity',aes(fill = ve_mean >= ve_mean[3] )) +
  geom_errorbar(varlimits,  width=0.25) +
  geom_hline(yintercept = 1/111,linetype="dashed") +
  scale_x_discrete(name ="Mode", limits=c(1:20),breaks =  seq(1,20,2)) +
  scale_y_continuous(expand = c(0, 0),limits=c(0,0.075),labels = percent,name = "Variance Explained", breaks=seq(0,0.075,length=5)) +
  #guides(fill = guide_legend(title = "diff >= 4SD")) +
  theme_classic(base_size = 30) +
  theme(legend.position = "none") 
p.var

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
  geom_text(size = 5, position = position_dodge(width=0.9), vjust= -3,color='grey')+
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
p.cv<-ggplot(cv,aes(modenum,cv_mean,label=round(cv_mean,2)))+ 
  geom_bar(width=.5,stat = 'identity', fill = '#80CBC4') +
  geom_errorbar(cvlimits,  width=0.25) +
  geom_text(size = 5, position = position_dodge(width=0.9), vjust= -0.8,color='grey')+
  scale_x_discrete(name ="Mode", limits=c(1:3)) +
  ylab('Cross-Validation Accuracy') +
  coord_cartesian(ylim=c(0.2,0.8)) +
  theme_classic(base_size = 30) +
  theme(legend.position="none")
p.cv

# build a null model 
behavior.perm <- rlply(1000,data$behavior[sample(nrow(data$behavior)),])
p3Km111.perm.cor<-sapply(behavior.perm, function(y_perm){ out<-ccaDWpermorder(data$brain,y_perm,p3Km111.gs$PENX,p3Km111.gs$PENY,k,data$behavior); out$cor} )
p3Km111.perm.cor<-sapply(behavior.perm, function(y_perm){ out<-ccaDWpermrank(data$brain,y_perm,p3Km111.gs$PENX,p3Km111.gs$PENY,k); out$cors} )

# permutation plots
perm <- as.data.frame(t(p3Km111.perm.cor))

## dim1
p.val <- (1+length(which(perm$V1 >= p3Km111.cor.mean[1]))) /1000
p <- ggplot(perm,aes(V1))+
  geom_histogram(binwidth = 0.005, fill = "blue", alpha = 0.5)  +
  geom_vline(xintercept = p3Km111.cor.mean[1], colour = "red", linetype = "longdash") +
  labs(x = "Correlations") +
  annotate("text", x = median(perm$V1), y = c(10,5),label = c("Permuted Data","(1000 times)"),size =10,colour = "black" ) +
  annotate("text",x = -Inf, y = Inf, hjust = -0.1,vjust = 1,label = paste("Cor=",round(p3Km111.cor.mean[1],2),", p<",p.val), size = 10, colour = "red" ) +
  theme_classic(base_size = 35) 
p

## dim2
p.val <- (1+length(which(perm$V2 >= p3Km111.cor.mean[2]))) /1000
p <- ggplot(perm,aes(V2))+
  geom_histogram(binwidth = 0.005, fill = "blue", alpha = 0.5)  +
  geom_vline(xintercept = p3Km111.cor.mean[2], colour = "red", linetype = "longdash") +
  labs(x = "Correlations") +
  annotate("text", x = median(perm$V2), y = c(10,5),label = c("Permuted Data","(1000 times)"),size =10,colour = "black" ) +
  annotate("text",x = -Inf, y = Inf, hjust = -0.1,vjust = 1, label = paste("Cor=",round(p3Km111.cor.mean[2],2),", p<",p.val), size = 10, colour = "red" ) +
  theme_classic(base_size = 35) 
p

## dim3
p.val <- (1+length(which(perm$V3 >= p3Km111.cor.mean[3]))) /1000
p <- ggplot(perm,aes(V3))+
  geom_histogram(binwidth = 0.005, fill = "blue", alpha = 0.5)  +
  geom_vline(xintercept = p3Km111.cor.mean[3], colour = "red", linetype = "longdash") +
  labs(x = "Correlations") +
  annotate("text", x = median(perm$V3), y = c(10,5),label = c("Permuted Data","(1000 times)"),size =10,colour = "black" ) +
  annotate("text",x = -Inf, y = Inf, hjust = -0.1,vjust = 1, label = paste("Cor=",round(p3Km111.cor.mean[3],2),", p<",p.val), size = 10, colour = "red" ) +
  theme_classic(base_size = 35) 
p


## bootstrap
bootid<-createResample(subjid$overall_psychopathology_4factor, list = T, times = 10)
brain_boot <- mclapply(bootid, function(id) data$brain[id,])
behavior_boot <- mclapply(bootid, function(id) data$behavior[id,])


p3Km111.org<-ccaDW(data$brain,data$behavior,0.8,0.4,3)

p3Km111.boot<- mclapply(seq_along(bootid),function(i) ccaDW(brain_boot[[i]],behavior_boot[[i]],0.8,0.4,3))

p3Km111.boot.ro<- mclapply(seq_along(bootid),function(i) reorderCCA(p3Km111.boot[[i]],p3Km111.org) )


