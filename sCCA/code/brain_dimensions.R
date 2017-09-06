u <- array(0,dim = 3410)
u[u.boot.plot[[1]]$fea] <- u.boot.plot[[1]]$load

dim1 <- t(apply(net.data,1, function(x) x * -u))
org.mat <- mask_mat(net.data)
dim1.ubyrid.mat <- mask_mat(dim1)

pwr_ft <- pwr_train_net_ft[,pwr.3k.train.idx]
pwr_ft.mat <- mask_mat(pwr_ft)

pwr_ft_hat <- pwr_train_net_ft[,pwr.3k.train.idx] - net.data # get the y_hat from original regression
pwr_ft_hat.mat <- mask_mat(pwr_ft_hat) #  ave mat and its plot
pwr_ft_hat.sign <- sign(pwr_ft_hat.mat$ave_mat) # get the y_hat signs

dim1.sign <- sign(dim1.ubyrid.mat$ave_mat) # get the u*x signs

prod.sign <- pwr_ft_hat.sign * dim1.sign # get the discordenant signs between y_hat and the residual projections


u2 <- matrix(data = rep(-u,2), nrow =2)

u2.mat <- mask_mat(u2)

u.final <- prod.sign * u2.mat$ave_mat     # apply the discordinant signs to u


dim1.sign.final <- dim1.sign * prod.sign # 

dim1.final <- dim1.sign.final * dim1.ubyrid.mat$ave_mat

#################################

u_by_x <- brain.plots[[1]]$mat * org.mat$ave_mat
levelplot(u_by_x,par.settings = BuRdTheme, at = seq(-max(u_by_x),max(u_by_x),length.out = 100))
levelplot(sign(u_by_x),par.settings = BuRdTheme())

prod.sign <- sign(u_by_x) * sign(pwr_ft.mat$ave_mat)
u1.final <- abs(brain.plots[[1]]$mat) * prod.sign

u_by_x_signed <- abs(u_by_x) * prod.sign


###################################
load_mat <- brain.plots[[3]]
ft_rgr_mat <- mask_mat(net.data)
ft_mat <-  mask_mat(pwr_train_net_ft[,pwr.3k.train.idx])

u.by.x <- load_mat$mat * ft_rgr_mat$ave_mat
delta.percent <- u.by.x / (ft_mat$ave_mat - ft_rgr_mat$ave_mat)
levelplot(delta.percent, par.settings = BuRdTheme(), at = seq(-10e-17,10e-17,length.out = 100))

