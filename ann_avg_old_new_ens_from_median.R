##Monthly and annual averages based on old and new data
rm(list=ls())

ann_avg_lk_old<- read.table(paste("D:/files_18feb/ann_avg_2decdata.txt")) #old data
med_old<- apply(ann_avg_lk_old,1, median) ##median

## annual averages above HadCRUT4 baseline
ANN_HDCLM_OLD<- matrix(0, nrow=167, ncol=10000)

		for(i in 1: 10000){
            	ANN_HDCLM_OLD[,i]<- as.numeric(ann_avg_lk_old[,i])- as.numeric(med_old) #new data
			}

ann_hdclm_old<- round(ANN_HDCLM_OLD,2)
	med_hdclm_old<- apply(ann_hdclm_old,1,median)
	LL_hdclm_old<- apply(ann_hdclm_old,1,function(x)quantile(x, prob=0.025)) 
	UL_hdclm_old<- apply(ann_hdclm_old,1,function(x)quantile(x, prob=0.975)) 

##Median ts and 95% confidence interval
ann_avg_lk<- read.table(paste("C:/Users/maryam/Dropbox/h_ths_ens/ann_avg_quntiles_med.txt"))
	med_hdclm<-(ann_avg_lk[,1])[1:167]
	LL_hdclm<- (ann_avg_lk[,2])[1:167]    
	UL_hdclm<- (ann_avg_lk[,3])[1:167]

X<- c(1:167)

X11()
#e<- paste("C:/Users/maryam/Dropbox/thesis/latex_thesis/figures/com_ens_ann_avg_new_two.png",sep="")
#png(e,width = 8, height = 10, units = 'in', res = 300)
## (a)
par(oma=c(0.7,0.7,0,0),mar=c(2,4,1,1), mfrow=c(2,1))
	plot(X,med_hdclm_old, type='n',ylab=substitute(paste(~Delta,"T wrt median")),
	xlab="", ylim=c(-0.4, 0.6),xaxt="n",yaxt="n",main="")                          
	axis(side = 1,at = seq(1,length(X), by=10), labels = seq(1850,2016, by=10),tck= -.02)
	axis(side = 2,at = seq(-0.4,0.5,by=0.1))
	abline(h=0, lty=2, col='black')
mtext("a",3, line=0.3, adj=-0.13, cex=1.2)

##old data 10 samples obs+cvg
polygon(c(X,rev(X)),c(LL_hdclm_old,rev(UL_hdclm_old)),col = adjustcolor("grey55", alpha.f=0.5), border = FALSE)
lines(X,med_hdclm_old, col="black")
legend("topleft", bty='n',col= c("black",adjustcolor("grey55", alpha.f=0.5),adjustcolor("purple", alpha.f=0.5)),
c("Median","Uncertainty (observational, coverage with 10 samples)", "Uncertainty (observational, coverage with 100 samples)"),lty=c(1,1,1), lwd=c(1,4,4))

##old data 100 samples obs+cvg
polygon(c(X,rev(X)),c(LL_hdclm_old,rev(UL_hdclm_old)),col = adjustcolor("purple", alpha.f=0.5), border = FALSE)
lines(X,med_hdclm_old, col="black")

	plot(X,med_hdclm, type='n',ylab=substitute(paste(~Delta,"T wrt median")),
	xlab="", ylim=c(-0.4, 0.6),xaxt="n",yaxt="n",main="")                          
	axis(side = 1,at = seq(1,length(X), by=10), labels = seq(1850,2016, by=10),tck= -.02)
	axis(side = 2,at = seq(-0.4,0.5,by=0.1))
	abline(h=0, lty=2, col='gray')
mtext("b",3, line=0.3, adj=-0.13, cex=1.2)
## (b)
##old data 10 samples obs+cvg
polygon(c(X,rev(X)),c(LL_hdclm_old,rev(UL_hdclm_old)),col = adjustcolor("grey55", alpha.f=0.5), border = FALSE)
lines(X,med_hdclm_old, col="black")

##new data
polygon(c(X,rev(X)),c(LL_hdclm[-c(168,169)],rev(UL_hdclm[-c(168,169)])),col = adjustcolor("purple", alpha.f=0.5), border = FALSE)
lines(X,med_hdclm[-c(168,169)], col="black")
legend("topleft", bty='n',col= c("black",adjustcolor("grey55", alpha.f=0.5),adjustcolor("purple", alpha.f=0.5)),
c("Median","Uncertainty (observational, coverage with 10 samples)", "Uncertainty (observational, coverage with 10 samples, parametric)"),lty=c(1,1,1), lwd=c(1,4,4))

dev.off()





ann_avg_lk<- read.table(paste("C:/Users/maryam/Dropbox/h_ths_ens/ann_avg_quntiles_med.txt"))


ann_avg_lk<- read.table("D:/lk_par_ens_new_data/ann_avg_new.txt")         #new data
##years
yr<- 1850:2016

##median
med_old<- apply(ann_avg_lk_old,1, median) 
med_new<- apply(ann_avg_lk,1, median)

## annual averages above HadCRUT4 baseline
ANN_HDCLM<- matrix(0, nrow=167, ncol=10000)
ANN_HDCLM_OLD<- matrix(0, nrow=167, ncol=10000)

		for(i in 1: 10000){
            	ANN_HDCLM[,i]<- as.numeric(ann_avg_lk[,i])- as.numeric(med_new)         #new data
			ANN_HDCLM_OLD[,i]<- as.numeric(ann_avg_lk_old[,i])- as.numeric(med_old) #new data
			}

ann_hdclm<- round(ANN_HDCLM,2)
ann_hdclm_old<- round(ANN_HDCLM_OLD,2)

##Median ts and 95% confidence interval
	med_hdclm<- apply(ann_hdclm,1,median)
	LL_hdclm<- apply(ann_hdclm,1,function(x)quantile(x, prob=0.025))
	UL_hdclm<- apply(ann_hdclm,1,function(x)quantile(x, prob=0.975))

	med_hdclm_old<- apply(ann_hdclm_old,1,median)
	LL_hdclm_old<- apply(ann_hdclm_old,1,function(x)quantile(x, prob=0.025)) 
	UL_hdclm_old<- apply(ann_hdclm_old,1,function(x)quantile(x, prob=0.975)) 

##Subset take first 10 LK ensembles
##observational coverage with 10 samples

AABB<- matrix(1:10000, ncol=100)
ind_matrix<- AABB[1:10,]
ind<- c(ind_matrix)

ann_avg_lk_s<- as.matrix(ann_avg_lk)[, ind]         #167x1000
ann_avg_lk_old_s<- as.matrix(ann_avg_lk_old)[, ind] #167x1000

##median
med_old_s<- apply(ann_avg_lk_old_s,1, median) 
med_new_s<- apply(ann_avg_lk_s,1, median)

## annual averages above HadCRUT4 baseline
ANN_HDCLM_s<- matrix(0, nrow=167, ncol=1000)
ANN_HDCLM_OLD_s<- matrix(0, nrow=167, ncol=1000)

		for(i in 1: 1000){
            	ANN_HDCLM_s[,i]<- as.numeric(ann_avg_lk_s[,i])- as.numeric(med_new_s)         #new data
			ANN_HDCLM_OLD_s[,i]<- as.numeric(ann_avg_lk_old_s[,i])- as.numeric(med_old_s) #new data
			}
ann_hdclm_s<- round(ANN_HDCLM_s,2)
ann_hdclm_old_s<- round(ANN_HDCLM_OLD_s,2)

##Median ts and 95% confidence interval
	med_hdclm_s<- apply(ann_hdclm_s,1,median)
	LL_hdclm_s<- apply(ann_hdclm_s,1,function(x)quantile(x, prob=0.025))
	UL_hdclm_s<- apply(ann_hdclm_s,1,function(x)quantile(x, prob=0.975))

	med_hdclm_old_s<- apply(ann_hdclm_old_s,1,median)
	LL_hdclm_old_s<- apply(ann_hdclm_old_s,1,function(x)quantile(x, prob=0.025)) 
	UL_hdclm_old_s<- apply(ann_hdclm_old_s,1,function(x)quantile(x, prob=0.975)) 
X<- c(1:nrow(ann_avg_lk))

X11()
e<- paste("C:/Users/maryam/Dropbox/thesis/latex_thesis/figures/com_ens_ann_avg_new_two.png",sep="")
png(e,width = 8, height = 10, units = 'in', res = 300)
## (a)
par(oma=c(0.7,0.7,0,0),mar=c(2,4,1,1), mfrow=c(2,1))
	plot(X,med_hdclm_s, type='n',ylab=substitute(paste(~Delta,"T wrt median")),
	xlab="", ylim=c(-0.4, 0.6),xaxt="n",yaxt="n",main="")                          
	axis(side = 1,at = seq(1,length(X), by=10), labels = seq(1850,2016, by=10),tck= -.02)
	axis(side = 2,at = seq(-0.4,0.5,by=0.1))
	abline(h=0, lty=2, col='gray')
mtext("a",3, line=0.3, adj=-0.13, cex=1.2)

##old data 10 samples obs+cvg
polygon(c(X,rev(X)),c(LL_hdclm_old_s,rev(UL_hdclm_old_s)),col = adjustcolor("grey55", alpha.f=0.5), border = FALSE)
lines(X,med_hdclm_old_s, col="black")
legend("topleft", bty='n',col= c("black",adjustcolor("grey55", alpha.f=0.5),adjustcolor("purple", alpha.f=0.5)),
c("Median","Uncertainty (observational, coverage with 10 samples)", "Uncertainty (observational, coverage with 100 samples)"),lty=c(1,1,1), lwd=c(1,4,4))

##old data 100 samples obs+cvg
polygon(c(X,rev(X)),c(LL_hdclm_old,rev(UL_hdclm_old)),col = adjustcolor("purple", alpha.f=0.5), border = FALSE)
lines(X,med_hdclm_old, col="black")

	plot(X,med_hdclm_s, type='n',ylab=substitute(paste(~Delta,"T wrt median")),
	xlab="", ylim=c(-0.4, 0.6),xaxt="n",yaxt="n",main="")                          
	axis(side = 1,at = seq(1,length(X), by=10), labels = seq(1850,2016, by=10),tck= -.02)
	axis(side = 2,at = seq(-0.4,0.5,by=0.1))
	abline(h=0, lty=2, col='gray')
mtext("b",3, line=0.3, adj=-0.13, cex=1.2)
## (b)
##old data 10 samples obs+cvg
polygon(c(X,rev(X)),c(LL_hdclm_old_s,rev(UL_hdclm_old_s)),col = adjustcolor("grey55", alpha.f=0.5), border = FALSE)
lines(X,med_hdclm_old_s, col="black")

##new data
polygon(c(X,rev(X)),c(LL_hdclm,rev(UL_hdclm)),col = adjustcolor("purple", alpha.f=0.5), border = FALSE)
lines(X,med_hdclm, col="black")
legend("topleft", bty='n',col= c("black",adjustcolor("grey55", alpha.f=0.5),adjustcolor("purple", alpha.f=0.5)),
c("Median","Uncertainty (observational, coverage with 10 samples)", "Uncertainty (observational, coverage with 10 samples, parametric)"),lty=c(1,1,1), lwd=c(1,4,4))

dev.off()


###as a sanity check , check range
DD<- t(apply(ann_hdclm_old,1,range))
EE<- t(apply(ann_hdclm_old_s,1,range))

dd<-(DD[,2] - DD[,1])
ee<- (EE[,2] - EE[,1])
ind<- abs(dd)>=abs(ee)
hh<- which(ind==FALSE)
FFRR<- cbind(DD,EE)[hh,]

##only these are slightly different maybe due to the rounding off
> FFRR[,2]-FFRR[,1]
[1] 0.53 0.43 0.24 0.23 0.21 0.17
> FFRR[,4]-FFRR[,3]
[1] 0.54 0.43 0.24 0.23 0.21 0.18