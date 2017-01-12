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
load("/data/joy/BBL/projects/xiaNetworkCca/sCCA/aim1/result/201701/power_regr_data.RData")
load("/data/joy/BBL/projects/xiaNetworkCca/sCCA/aim1/result/201701/med_cv_split.RData")
net.data <- unname(power.rgr.train)
med.data <- as.matrix(med.train[,2:112])
data <- list(brain = net.data, behavior = med.data)

## create sub-CV splits
subjid <- read.csv('./result/201701/trainsample.csv')

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
p3Km111.ve <- sapply(seq_along(trainid), function(i) p3Km111.cca[[i]]$d^2/sum(p3Km111.cca[[i]]$d^2))

# plot ve
# variance explained plots
varE <- data.frame( modenum = as.factor(1:111),ve_mean = rowMeans(p3Km111.ve),ve_se = rowSds(p3Km111.ve)/sqrt(dim(p3Km111.ve)[2]),vediff = c(abs(diff(rowMeans(p3Km111.ve))),0))
varlimits <- aes(ymax = varE$ve_mean + varE$ve_se, ymin=varE$ve_mean - varE$ve_se)
p<-ggplot(varE,aes(modenum,ve_mean))+ 
  geom_bar(stat = 'identity',aes(fill = vediff >= mean(vediff)+3*sd(vediff) )) +
  geom_errorbar(varlimits,  width=0.25) +
  geom_hline(yintercept = 1/111,linetype="dashed") +
  scale_x_discrete(name ="Mode", limits=c(1:20),breaks =  seq(1,20,2)) +
  scale_y_continuous(limits=c(0,0.02),labels = percent,name = "Variance Explained", breaks=seq(0,0.07,length=5)) +
  guides(fill = guide_legend(title = "diff >= 3SD")) +
  theme_classic(base_size = 24) +
  theme(legend.position = c(.5, 0.95),legend.direction = 'horizontal') 
p

boxplot(t(p3Km111.ve))
k=4;


