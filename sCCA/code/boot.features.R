

data.train.boot.brain <- lapply(u.boot.plot, function(x) data.train$brain[,x$fea])
data.train.boot.behav <- lapply(v.boot.plot, function(x) data.train$behavior[,x$fea])
boot.sCCA.train <- lapply(1:3, function(i) cc(data.train.boot.brain[[i]], data.train.boot.behav[[i]]))
sCCA.train <- ccaDW(data.train$brain, data.train$behavior, 0.4, 0.4, 5)

data.test.boot.brain <- lapply(u.boot.plot, function(x) data.test$brain[,x$fea])
data.test.boot.behav <- lapply(v.boot.plot, function(x) data.test$behavior[,x$fea])                          
boot.sCCA.test <- lapply(1:3, function(i) cc(data.test.boot.brain[[i]], data.test.boot.behav[[i]]))
sCCA.test <- ccaDW(data.test$brain, data.test$behavior, 0.4, 0.4, 5)

cv_u <- (lapply(1:3, function(i) cc(boot.sCCA.train[[i]]$xcoef,boot.sCCA.test[[i]]$xcoef)))
cv_u <- (lapply(1:3, function(i) cor(boot.sCCA.train[[i]]$xcoef[,1],boot.sCCA.test[[i]]$xcoef[,1])))

cv_v <- (sapply(1:3, function(i) cor.test(boot.sCCA.train[[i]]$ycoef[,1],boot.sCCA.test[[i]]$ycoef[,1])))

#permutation
cv_u.perm <- lapply(1:100, 
                    function(x) {
#boot.sCCA.test.perm <- lapply(1:3, function(i) cc(data.test.boot.brain[[i]][sample(nrow(data.test.boot.brain[[i]])),], data.test.boot.behav[[i]][sample(nrow(data.test.boot.behav[[i]])),]  ))
                      #out <- sapply(1:3, function(i) cor(boot.sCCA.train[[i]]$xcoef[,1], boot.sCCA.test[[i]]$xcoef[,1] [sample(length(boot.sCCA.test[[i]]$xcoef[,1]))]   )   )
                      out <- sapply(1:3, function(i) {cc.perm<-cc(boot.sCCA.train[[i]]$xcoef, boot.sCCA.test[[i]]$xcoef[sample(length(boot.sCCA.test[[i]]$xcoef[,1])),]); cc.perm$cor }   )
                      
}
)

## more selection
data.train.boot.brain <- lapply(u.boot.plot, function(x) data.train$brain[,x$fea[order(-abs(x$high))[1:10]]])
data.train.boot.behav <- lapply(v.boot.plot, function(x) data.train$behavior[,x$fea])
boot.sCCA.train <- lapply(1:3, function(i) cc(data.train.boot.brain[[i]], data.train.boot.behav[[i]]))
sCCA.train <- ccaDW(data.train$brain, data.train$behavior, 0.4, 0.4, 5)

data.test.boot.brain <- lapply(u.boot.plot, function(x) data.test$brain[,x$fea[order(-abs(x$high))[1:10]]])
data.test.boot.behav <- lapply(v.boot.plot, function(x) data.test$behavior[,x$fea])                          
boot.sCCA.test <- lapply(1:3, function(i) cc(data.test.boot.brain[[i]], data.test.boot.behav[[i]]))
sCCA.test <- ccaDW(data.test$brain, data.test$behavior, 0.4, 0.4, 5)

cv_u <- (sapply(1:3, function(i) cor.test(boot.sCCA.train[[i]]$xcoef[,1],boot.sCCA.test[[i]]$xcoef[,1])))
cv_v <- (sapply(1:3, function(i) cor.test(boot.sCCA.train[[i]]$ycoef[,1],boot.sCCA.test[[i]]$ycoef[,1])))




candmod = c(1,5,6,2)
#set up the BT samples.
bootnum <- 1000
bootid<-createResample(pwr_test_qa$overall_psychopathology_4factor, list = T, times = bootnum)
brain_boot_test <- lapply(bootid, function(id) data.test$brain[id,])
behavior_boot_test <- lapply(bootid, function(id) data.test$behavior[id,])


predictY <-function(train_model,candmod, testX, testY) {
  Y_hat <- (testX %*% train_model$u[,candmod]) %*% ginv(train_model$v[,candmod])
  cv_score<- cancor(Y_hat,testY)
}

predictY_boot <- t(sapply( seq_along(brain_boot_test), function(i) {out<-predictY(sCCA.train,candmod,brain_boot_test[[i]],behavior_boot_test[[i]]); out$cor} ) )
boxplot(predictY_boot)
  

#permutation
behav_test_perm <- rlply(1000,data.test$behavior[sample(nrow(data.test$behavior)),])
behav_test_hat <- (data.test$brain %*% sCCA.train$u[,candmod]) %*% ginv(sCCA.train$v[,candmod])

cv.perm <- t(sapply(behav_test_perm,function(x) {out <- cancor(x, behav_test_hat); out$cor}))

p.val <- sapply(1:4, function(i) length(which(cv.perm[,i]  >= mean(predictY_boot[,i]))))

#cv plot
cv <- data.frame(modenum = as.factor(1:4), cv_acc = apply(predictY_boot,2,mean))
cv.long <- melt(predictY_boot)
colnames(cv.long) <- c("bootid","mode","cv_acc")
cv.long$mode <- as.factor(cv.long$mode)

cv.perm.long <- melt(cv.perm)
colnames(cv.perm.long) <- c("bootid","mode","cv_perm")
cv.perm.long$mode <- as.factor(cv.perm.long$mode)

p.cv <- ggplot()+ 
  geom_boxplot(aes(cv.long$mode,cv.long$cv_acc),fill =c("#8E24AA40","#1E88E540","#FF6F0040","#D32F2F40"), color =c("#8E24AA","#1E88E5","#FF6F00","#D32F2F") ) + 
  scale_x_discrete(name ="Mode", limits=c(1:4)) +
  scale_color_manual(values = c("#8E24AA","#1E88E5","#FF6F00","#D32F2F")) +
  ylab('Yscore & Yscore_hat Correlation') +
  coord_cartesian(ylim=c(0.45,0.9)) +
  theme_classic(base_size = 25) +
  
  geom_dotplot(aes(cv.perm.long$mode, cv.perm.long$cv_perm), color = 'grey', stackdir = 'centerwhole',binaxis = 'y', dotsize = 0.1, stackratio = 0.3)
  #geom_boxplot(aes(cv.perm.long$mode, cv.perm.long$cv_perm), color = 'grey')
p.cv
plotname <- paste('~/Google Drive/CEDRIC/PNC_CCA/Figure Resources/','Fig8_','cv_real','.pdf',sep="")
pdf(file = plotname, width = 6, height = 5,useDingbats=F)
print(p.cv)
dev.off()

#cv plot verson 2
predictY_boot_null <- predictY_boot - cv.perm
cv.long <- melt(predictY_boot_null)
colnames(cv.long) <- c("bootid","mode","cv_acc")
cv.long$mode <- as.factor(cv.long$mode)

p.cv.null <- ggplot()+ 
  geom_boxplot(aes(cv.long$mode,cv.long$cv_acc),fill =c("#8E24AA40","#1E88E540","#FF6F0040","#D32F2F40"), color =c("#8E24AA","#1E88E5","#FF6F00","#D32F2F") ) + 
  scale_x_discrete(name ="Mode", limits=c(1:4)) +
  scale_color_manual(values = c("#8E24AA","#1E88E5","#FF6F00","#D32F2F")) +
  ylab('Yscore & Yscore_hat Correlation') +
  coord_cartesian(ylim=c(0,0.4)) +
  theme_classic(base_size = 25)
p.cv.null
plotname <- paste('~/Google Drive/CEDRIC/PNC_CCA/Figure Resources/','Fig8_','cv_real_null','.pdf',sep="")
pdf(file = plotname, width = 6, height = 5,useDingbats=F)
print(p.cv.null)
dev.off()






#perm plots
pval.perm.adj <- p.val + 1/1000
for (i in 1:4){
  colpal <- c("Purples","Blues","Oranges","Reds")
  colkey <- brewer.pal(8,colpal[i])[c(2:3,7:8)]
  
  cv_cor = mean(predictY_boot[,i])
  perm_file <- data.frame(cor_perm = cv.perm[,i])
  p.perm <- ggplot(perm_file,aes(cor_perm,fill = cor_perm > cv_cor))+
    geom_histogram(binwidth = 0.01, alpha = 0.5) +
    scale_fill_manual(values = colkey[2:3]) +
    geom_vline(xintercept = cv_cor, colour = colkey[4], linetype = "longdash") +
    labs(x = "Y_test and Y_predict Correlations") +
    #annotate("text", x = median(perm_file$cor_perm,na.rm = T), y = c(15,5),label = c("Permuted Data","(1000 times)"),size =6,colour = "black" ) +
    annotate("text",x = -Inf, y = Inf, hjust = -0.1,vjust = 1,label = paste("p<",round(pval.perm.adj[i],3)), size = 8, colour = "black",fontface ="italic" ) +
    theme_classic(base_size = 30) + 
    scale_y_continuous(expand = c(0, 0)) +
    theme(legend.position="none")
  
  plotname <- paste('~/Google Drive/CEDRIC/PNC_CCA/Figure Resources/','Fig8_','cv_perm_',i,'.pdf',sep="")
  pdf(file = plotname, width = 6, height = 5,useDingbats=F)
  print(p.perm)
  dev.off()
}
