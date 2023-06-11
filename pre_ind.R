##Figure-1 for new data
rm(list=ls())
library(sROC)

#################################################################################
#####~~~~Figure-1a, annual averages above HadCRUT4 baseline i.e. 1961 to 1990####
#################################################################################
ann_2016<- read.table("D:/lk_par_ens_new_data/del/2016_ann_avg_new.txt")

##annual averages 167x10000
ann_avg_lk<- read.table("D:/lk_par_ens_new_data/ann_avg_new.txt") 

##years
yr<- 1850:2016

##index of baseline years
ind_1961_1990<- which(yr>=1961 & yr<=1990)

##annual averages of baseline years
ann_avg_1961_1990<- ann_avg_lk[ind_1961_1990,] #30 by 10000
 
##HadCRUT baseline
CLM_HD<- apply(ann_avg_1961_1990,2, mean)  #vector of length 10000
clm_hd<- as.numeric(CLM_HD)                #get rid of column names 

## annual averages above HadCRUT4 baseline
ANN_HDCLM<- matrix(0, nrow=167, ncol=10000)
		for(i in 1: 167){
            	ANN_HDCLM[i,]<- as.numeric(ann_avg_lk[i,])- clm_hd
			}
##Median ts and 95% confidence interval
	med_hdclm<- apply(ANN_HDCLM,1,median)
	LL_hdclm<- apply(ANN_HDCLM,1,function(x)quantile(x, prob=0.025))
	UL_hdclm<- apply(ANN_HDCLM,1,function(x)quantile(x, prob=0.975))

######################################################################################################################################
######Figure-1b~~~~~Draw samples from the CDF of the difference between a "pre-industrial" 1400-1800 baseline and 1850-1900~~~~~######
######################################################################################################################################
##~~CDF of the difference between a "pre-industrial" 1400-1800 baseline and 1850-1900. (details below)
T<- read.table("C:/Users/maryam/Dropbox/PhD UCL/kriging_and_FEOFs/post_upgrade/paper1.5oC_r_files/data_andrew/cdf_All_Forcings.txt")
			
set.seed(102)
##~~i)Obtain samples having the CDF (T)above (details below)
	
	#A probability density function of the difference between a "pre-industrial" 1400-1800 baseline and 1850-1900.
	pdf_a<- diff(T[,2])    #/diff(T[,1])
	freq<- 23*pdf_a        #multiply with 51 so that after rounding we have 10,000
				     #sum(round(freq))
		                 #range(T[,1])

	##Frequency distribution
	breaks<- seq(min(T[,1]),max(T[,1]), by=0.01)
	Cl<- cbind(seq(breaks[1],breaks[length(breaks)-1],by=(T[,1][2]-T[,1][1])), seq(breaks[2],breaks[length(breaks)],by=(T[,1][2]-T[,1][1])))
	XX<- apply(Cl,1,mean)
	YY<- cbind(XX,round(freq))

##samples of the difference of climatology 1400-1800 to 1850-1900 
clm_a<- rep(XX,round(freq))   	#median(clm_a)=0.075 as is the median of CDF above 
						#hist(clm_a)
median(clm_a)
length(clm_a)
			X11()
			par(mfrow=c(1,2))
			plot(T[,1],T[,2])
			hist(clm_a, main=paste("median=",round(median(clm_a),4)))
			abline(v=median(clm_a),col=2)

##~~ii)Propagate uncertainties in the empirical CDF

##consider samples from the empirical CDF (of the difference)
x<- clm_a

##obtain smooth kernel estimate
x.CDF <- kCDF(x, ngrid=length(T[,1]))
ll<- CI.CDF(x.CDF)$Fhat.lower
uu<- CI.CDF(x.CDF)$Fhat.upper
xhat<- x.CDF$x
fhat<- x.CDF$Fhat

	##data based on xhat and fhat
	pdf_1<- diff(fhat)
	freq_1<- 23*pdf_1    #so that the final length is 100   

	breaks_1<- seq(min(xhat),max(xhat), by=c(xhat[2]-xhat[1]))
	Cl_1<- cbind(seq(breaks_1[1],breaks_1[length(breaks_1)-1],by=c(xhat[2]-xhat[1])), seq(breaks_1[2],breaks_1[length(breaks_1)],by=c(xhat[2]-xhat[1])))
	XX_1<- apply(Cl_1,1,mean)
	YY_1<- cbind(XX_1,round(freq_1))
	clm_1<- rep(XX_1,round(freq_1))
	median(clm_1)
	length(clm_1)

##Generate bootstrap smooth kernel estimate (because we need data which would be available this way) 
B<- 200
n<- length(clm_1)
x<- clm_1

xb<- matrix(0, nrow=n,ncol=B)

mpf<- matrix(0,nrow=B,ncol=length(XX_1))
      mn<-function(x,m){
      a<- table(x)
      b<-rep(0,m)
      b[as.integer(names(a))]<-a
      return(b)
      }

for(b in 1:B){
    xb[,b]<- sample(1:length(XX_1),n,replace=TRUE,prob=c(pdf_1))
    mpf[b,]<- mn(xb[,b],length(XX_1))/n
    }

MPF<- apply(mpf,1,cumsum)
med<- apply(MPF,1,median)
LL<- apply(MPF,1,function(x)quantile(x,probs=0.025))
UU<- apply(MPF,1,function(x)quantile(x,probs=0.975))

yy3a<- XX_1[c(xb)]   			#median(yy3a): 0.076, quantile(yy3a, 0.025):-0.13, quantile(yy3a, 0.975):0.29

	X11()
	par(mfrow=c(2,2))
	plot(T[,1],T[,2])
	hist(clm_a, main=paste("median=",round(median(clm_a),4)))
	abline(v=median(clm_a),col=2)
	plot(XX_1,LL,ylim=c(0,1),type="n", xlab="",ylab="",main="")
	lines(XX_1,LL)
	lines(XX_1,UU)
	lines(XX_1,med, col=2)
	legend("topleft", bty="n", c("95% CI","median"),col=c(1,2),lty=c(1,1) )
	plot(x.CDF)

####~~~~global annual averages above pre-industrial

##########################################################################
######Figure-1c~~~~~global annual averages above pre-industrial~~~~~######
##########################################################################

##Baseline 1850-1900
ind_1850_1900<- which(yr>=1850 & yr<=1900)
CLM_O<- ann_avg_lk[ind_1850_1900,]
clm_o<- as.numeric(apply(CLM_O,2, mean))

ANN_ABV_1850_1900<- matrix(0, nrow=167, ncol=10000)

		for(i in 1: 167){
            	ANN_ABV_1850_1900[i,]<- as.numeric(ann_avg_lk[i,])- clm_o
			}
smp<- 150

ANN_PREIND<- matrix(0, nrow=167, ncol=(smp*10000) )

	for(i in 1: 167){
			M1<- ANN_ABV_1850_1900[i,]
			M2<- rep(M1, each=smp)
			M3<- sample(yy3a,smp)
			ANN_PREIND[i,]<- rep(M3,10000)+M2
			}

ANN_cvg<- matrix(0, nrow=167, ncol=(smp*10000) )

	for(i in 1: 167){
			M1<- ANN_ABV_1850_1900[i,]
			M2<- rep(M1, each=smp)
			M3<- sample(yy3a,smp)
			ANN_PREIND[i,]<- rep(M3,10000)+M2
			}

#> summary(ANN_PREIND[1,])
#      Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-0.48409 -0.06658  0.02264  0.02859  0.11785  0.65023 

#> summary(ANN_HDCLM[1,])
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-0.637874 -0.491312 -0.439524 -0.420520 -0.383053  0.006207 

#quantile(ANN_PREIND[167,],0.025)  #1.005806 
#quantile(ANN_PREIND[167,],0.975)  #1.530255 
#sum(ANN_PREIND[167,]>=1.5)/length(ANN_PREIND[167,]) #0.04031733

	##Median ts and 95% CI
	med_ann_preind<- apply(ANN_PREIND,1,median)
	LL_ann_preind<- apply(ANN_PREIND,1,function(x)quantile(x, prob=0.025))
	UL_ann_preind<- apply(ANN_PREIND,1,function(x)quantile(x, prob=0.975))

X<- c(1:nrow(ann_avg_lk))

##For thesis
X11()
e<- paste("C:/Users/maryam/Dropbox/thesis/latex_thesis/figures/pre_ind_KCDF.png",sep="")
png(e,width = 8.5, height = 10, units = 'in', res = 300)
par(mfrow=c(2,1),mar=c(4,4,1.5,1.5), oma=c(0.5,0.5,0,0))  
plot(T[,1],T[,2], ylab="", xlab="", main="(a)", cex.main=1.5, cex.axis=1.5)
plot(x.CDF, ylab="", xlab="", main="(b)", cex.main=1.5, cex.axis=1.5)
dev.off()

X11()
e<- paste("C:/Users/maryam/Dropbox/thesis/latex_thesis/figures/bs_pre_ind.png", sep="")
png(e,width = 5, height = 5, units = 'in', res = 300)
par(oma=c(0,0,0,0),mar=c(2.5,4,0.5,0.5))

hist(yy3a, breaks=XX_1, prob=TRUE,xaxt="n",yaxt="n",xlim=range(-0.3,0.5),col="blue",
main="",ylab="probability",  xlab="") 
axis(1, at=  round(seq(-0.3,0.5, by=0.1),1))
axis(2, at= seq(0,3.5,by=0.5))
dev.off()


X11()
e<- paste("C:/Users/maryam/Dropbox/thesis/latex_thesis/figures/2016_pre_ind.png", sep="")
png(e,width = 5, height = 5, units = 'in', res = 300)
par(oma=c(0,0,0,0),mar=c(4,4,0.5,0.5))

hist(as.numeric(ANN_PREIND[167,]), prob=TRUE,breaks=seq(0.8,1.78,by=0.02),xaxt="n",yaxt="n",ylim=c(0,3.5),
col="blue", xlab=substitute(paste(~Delta,"T wrt preindustrial")), ylab="probability", 
main= "") 
axis(1, at= seq(0.8,1.76, by=0.1))
axis(2, at= seq(0,3.5,by=0.5))
dev.off()


X11()
e<- paste("C:/Users/maryam/Dropbox/thesis/latex_thesis/figures/ann_pre_ind.png", sep="")
png(e,width = 7, height = 8, units = 'in', res = 300)
par(mfrow=c(2,1),mar=c(4,4,1.5,0.5), oma=c(0.5,0.5,0,0))  

##Figure-1a
	plot(X,med_hdclm, type='n',ylab=substitute(paste(~Delta,"T wrt 1961-90")),
	xlab="", ylim=c(-1, 2),xaxt="n",yaxt="n",main="(a)")                          
	axis(side = 1,at = seq(1,length(X), by=10), labels = seq(1850,2016, by=10),tck= -.02)
	axis(side = 2,at = seq(-1,2,by=0.2))
	abline(h=0, lty=2, col='gray')

polygon(c(X,rev(X)),c(LL_hdclm,rev(UL_hdclm)),col = "grey75", border = FALSE)
lines(X,med_hdclm, col="black")
legend(-7,2.1, bty='n',col= c("black","grey75"),c("Median","Uncertainty"),lty=c(1,1,1), lwd=c(1,4,4))
#mtext("a",3, line=0.3, adj=-0.13, cex=1.2)

##Figure-1c
plot(X,med_ann_preind, type='n',ylab=substitute(paste(~Delta,"T wrt preindustrial")),
	 xlab="Time series of global annual mean temperatures", 
       ylim=c(-1, 2),xaxt="n",yaxt="n",main="(b)")                          
	 axis(side = 1,at = seq(1,length(X), by=10), labels = seq(1850,2016, by=10),tck= -.02)
	 axis(side = 2,at = seq(-1,2,by=0.2))
	 abline(h=0, lty=2, col='gray')
	 abline(h=1.5, col='red')
	 abline(h=2, col='red')
polygon(c(X,rev(X)),c(LL_ann_preind,rev(UL_ann_preind)),col = "grey75", border = FALSE)
lines(X,med_ann_preind, col="black")
#mtext("b",3, line=1.5, adj=-.13, cex=1.2)
dev.off()




X11()
e<- paste("C:/Users/maryam/Dropbox/PhD UCL/kriging_and_FEOFs/post_upgrade/paper1.5oC_r_files/figures/fig1.png",sep="")
png(e,width = 12, height = 8, units = 'in', res = 300)
	grid<- matrix(c(1,1,1,3,3, 
			    1,1,1,3,3,
                	    1,1,1,3,3,
		    	    1,1,1,3,3,
			    2,2,2,4,4,
		    	    2,2,2,4,4,
		    	    2,2,2,4,4,
		    	    2,2,2,4,4),nrow=8, ncol=5, byrow=TRUE)
layout(grid)
par(mar=c(4,4,2,0), oma=c(0.5,0.5,0,0))  

##Figure-1a
	plot(X,med_hdclm, type='n',ylab=substitute(paste(~Delta,"T wrt 1961-90")),
	xlab="", ylim=c(-1, 2),xaxt="n",yaxt="n", main="(a)")                          
	axis(side = 1,at = seq(1,length(X), by=10), labels = seq(1850,2016, by=10),tck= -.02)
	axis(side = 2,at = seq(-1,2,by=0.2))
	abline(h=0, lty=2, col='gray')

polygon(c(X,rev(X)),c(LL_hdclm,rev(UL_hdclm)),col = "grey75", border = FALSE)
lines(X,med_hdclm, col="black")
legend(-7,2.1, bty='n',col= c("black","grey75"),c("Median","Uncertainty"),lty=c(1,1,1), lwd=c(1,4,4))
#mtext("a",3, line=0.3, adj=-0.13, cex=1.2)

##Figure-1c
plot(X,med_ann_preind, type='n',ylab=substitute(paste(~Delta,"T wrt preindustrial")),
	 xlab="Time series of global annual mean temperatures", 
       ylim=c(-1, 2),xaxt="n",yaxt="n",main="(c)")                          
	 axis(side = 1,at = seq(1,length(X), by=10), labels = seq(1850,2016, by=10),tck= -.02)
	 axis(side = 2,at = seq(-1,2,by=0.2))
	 abline(h=0, lty=2, col='gray')
	 abline(h=1.5, col='red')
	 abline(h=2, col='red')
polygon(c(X,rev(X)),c(LL_ann_preind,rev(UL_ann_preind)),col = "grey75", border = FALSE)
lines(X,med_ann_preind, col="black")
#mtext("c",3, line=1.5, adj=-.13, cex=1.2)

##Figure-1c
hist(yy3a, breaks=XX_1, prob=TRUE,xaxt="n",yaxt="n",xlim=range(-0.3,0.5),col="blue",
main="(b)",ylab="probability",  xlab="") #substitute(paste("(b)Preindustrial correction"))
axis(1, at=  round(seq(-0.3,0.5, by=0.1),1))
axis(2, at= seq(0,3.5,by=0.5))
	#abline(v=median(yy3a), col=2)
#mtext("b",3, line=0.3, adj=-0.17, cex=1.2)

#hist(as.numeric(ANN_PREIND[167,]), prob=TRUE,breaks=seq(0.7,1.85,by=0.05),xaxt="n",
hist(as.numeric(ANN_PREIND[167,]), prob=TRUE,breaks=seq(0.8,1.76,by=0.02),xaxt="n",yaxt="n",ylim=c(0,3.5),
col="blue", xlab=substitute(paste(~Delta,"T wrt preindustrial")), ylab="probability", 
main= "(d)") #substitute(paste(2016~Delta,"T wrt preindustrial"))
axis(1, at= seq(0.8,1.76, by=0.1))
axis(2, at= seq(0,3.5,by=0.5))
#mtext("d",3, line=1.5, adj=-0.17, cex=1.2)
dev.off()


	

##For news item
X11()
e<- paste("C:/Users/maryam/Dropbox/PhD UCL/kriging_and_FEOFs/post_upgrade/paper1.5oC_r_files/figures/fig1_news_item.png",sep="")
png(e,width = 7, height = 8, units = 'in', res = 300)
	
par(mfrow=c(2,1),mar=c(4,4,1.5,1.5), oma=c(0.5,0.5,0,0))  

##Figure-1a
	plot(X,med_hdclm, type='n',ylab=substitute(paste(~Delta,"T wrt 1961-90")),
	xlab="", ylim=c(-1, 2),xaxt="n",yaxt="n",main="")                          
	axis(side = 1,at = seq(1,length(X), by=10), labels = seq(1850,2016, by=10),tck= -.02)
	axis(side = 2,at = seq(-1,2,by=0.2))
	abline(h=0, lty=2, col='gray')

polygon(c(X,rev(X)),c(LL_hdclm,rev(UL_hdclm)),col = "grey75", border = FALSE)
lines(X,med_hdclm, col="black")
legend(-7,2.1, bty='n',col= c("black","grey75"),c("Median","Uncertainty"),lty=c(1,1,1), lwd=c(1,4,4))
mtext("a",3, line=0.3, adj=-0.13, cex=1.2)

##Figure-1c
plot(X,med_ann_preind, type='n',ylab=substitute(paste(~Delta,"T wrt preindustrial")),
	 xlab="Time series of global annual mean temperatures", 
       ylim=c(-1, 2),xaxt="n",yaxt="n",main="")                          
	 axis(side = 1,at = seq(1,length(X), by=10), labels = seq(1850,2016, by=10),tck= -.02)
	 axis(side = 2,at = seq(-1,2,by=0.2))
	 abline(h=0, lty=2, col='gray')
	 abline(h=1.5, col='red')
	 abline(h=2, col='red')
polygon(c(X,rev(X)),c(LL_ann_preind,rev(UL_ann_preind)),col = "grey75", border = FALSE)
lines(X,med_ann_preind, col="black")
mtext("b",3, line=1.5, adj=-.13, cex=1.2)
dev.off()


##For Recover brochure
X11()
e<- paste("C:/Users/maryam/Dropbox/PhD UCL/kriging_and_FEOFs/post_upgrade/paper1.5oC_r_files/figures/fig2_.png",sep="")
png(e,width = 7, height = 8, units = 'in', res = 300)
	
par(mfrow=c(2,1),mar=c(4,4,1.5,1.5), oma=c(0.5,0.5,0,0))  

##Figure-1a
	plot(X,med_hdclm, type='n',ylab=substitute(paste(~Delta,"T wrt 1961-90")),
	xlab="", ylim=c(-1, 2),xaxt="n",yaxt="n",main="")                          
	axis(side = 1,at = seq(1,length(X), by=10), labels = seq(1850,2016, by=10),tck= -.02)
	axis(side = 2,at = seq(-1,2,by=0.2))
	abline(h=0, lty=2, col='gray')

polygon(c(X,rev(X)),c(LL_hdclm,rev(UL_hdclm)),col = "grey75", border = FALSE)
lines(X,med_hdclm, col="black")
legend(-7,2.1, bty='n',col= c("black","grey75"),c("Median","Uncertainty"),lty=c(1,1,1), lwd=c(1,4,4))
mtext("a",3, line=0.3, adj=-0.13, cex=1.2)
dev.off()



#########################################################################################################



ANN_HDCLM<- ann_avg_lk - matrix(rep(clm_hd, each=nrow(ann_avg_lk)),ncol=10000)
ann_hdclm<- round(ANN_HDCLM,2)

	##Median ts and 95% confidence interval
	med_hdclm<- apply(ann_hdclm,1,median)
	LL_hdclm<- apply(ann_hdclm,1,function(x)quantile(x, prob=0.025))
	UL_hdclm<- apply(ann_hdclm,1,function(x)quantile(x, prob=0.975))






ann_avg_lk<- read.table(paste("D:/files_18feb/ann_avg_2decdata.txt")) 		#167x10000
yr<- 1850:2016     										#years

ind_1850_1900<- which(yr>=1850 & yr<=1900)
CLM_O<- ann_avg_lk[ind_1850_1900,]
clm_o<- as.numeric(apply(CLM_O,2, mean))

ind_1961_1990<- which(yr>=1961 & yr<=1990)
CLM_HD<- ann_avg_lk[ind_1961_1990,] 
clm_hd<- apply(CLM_HD,2, mean)




ann_avg_2016<- read.table("D:/new_data/2016_ann_avg_new.txt")
mon_avg_lk<- read.table("D:/new_data/upto668months_new.txt")   #TIME.t[668]= "1905-08-16 12:00:00 GMT"
MON_AVG_LK<- mon_avg_lk[1:612,]					   #TIME.t[612]= "1900-12-16 12:00:00 GMT"

a_a<- matrix(NA, nrow= (612)/12)
S_y<- 1850
for(i in 1:ncol(MON_AVG_LK)){
		 M_A<- MON_AVG_LK[,i]
		 A_A<- ann_avg(M_A,S_y)
 		 a_a<- cbind(a_a, A_A)	
}

CLM_O<- a_a[,-1]
clm_o<- as.numeric(apply(CLM_O,2, mean))















	
library(sROC)

#####################################################################################################################################
######Sectio-1~~~~~Draw samples from the CDF of the difference between a "pre-industrial" 1400-1800 baseline and 1850-1900~~~~~######
#####################################################################################################################################

##~~CDF of the difference between a "pre-industrial" 1400-1800 baseline and 1850-1900. (details below)
T<- read.table("C:/Users/maryam/Dropbox/PhD UCL/kriging_and_FEOFs/post_upgrade/paper1.5oC_r_files/data_andrew/cdf_All_Forcings.txt")
			
set.seed(102)
##~~Obtain samples having the CDF (T)above (details below)
		#A probability density function of the difference between a "pre-industrial" 1400-1800 baseline and 1850-1900.
	pdf_a<- diff(T[,2])    #/diff(T[,1])
	freq<- 23*pdf_a        #multiply with 51 so that after rounding we have 10,000
		#sum(round(freq))
		#range(T[,1])

	##Frequency distribution
	breaks<- seq(min(T[,1]),max(T[,1]), by=0.01)
	Cl<- cbind(seq(breaks[1],breaks[length(breaks)-1],by=(T[,1][2]-T[,1][1])), seq(breaks[2],breaks[length(breaks)],by=(T[,1][2]-T[,1][1])))
	XX<- apply(Cl,1,mean)
	YY<- cbind(XX,round(freq))

##samples of the difference of climatology 1400-1800 to 1850-1900 
clm_a<- rep(XX,round(freq))   	#median(clm_a)=0.075 as is the median of CDF above 
						#hist(clm_a)
median(clm_a)
length(clm_a)

		X11()
		par(mfrow=c(1,2))
		plot(T[,1],T[,2])
		hist(clm_a, main=paste("median=",round(median(clm_a),4)))
		abline(v=median(clm_a),col=2)

x<- clm_a
x.CDF <- kCDF(x, ngrid=length(T[,1]))
x.CDF
ll<- CI.CDF(x.CDF)$Fhat.lower
uu<- CI.CDF(x.CDF)$Fhat.upper
xhat<- x.CDF$x
fhat<- x.CDF$Fhat

	##data based on xhat and fhat
	pdf_1<- diff(fhat)
	freq_1<- 23*pdf_1       

	breaks_1<- seq(min(xhat),max(xhat), by=c(xhat[2]-xhat[1]))
	Cl_1<- cbind(seq(breaks_1[1],breaks_1[length(breaks_1)-1],by=c(xhat[2]-xhat[1])), seq(breaks_1[2],breaks_1[length(breaks_1)],by=c(xhat[2]-xhat[1])))
	XX_1<- apply(Cl_1,1,mean)
	YY_1<- cbind(XX_1,round(freq_1))
	clm_1<- rep(XX_1,round(freq_1))
	median(clm_1)
	length(clm_1)

B<- 10000
n<- length(clm_1)
x<- clm_1

xb<- matrix(0, nrow=n,ncol=B)

mpf<- matrix(0,nrow=B,ncol=length(XX_1))
      mn<-function(x,m){
      a<- table(x)
      b<-rep(0,m)
      b[as.integer(names(a))]<-a
      return(b)
      }

for(b in 1:B){
    xb[,b]<- sample(1:length(XX_1),n,replace=TRUE,prob=c(pdf_1))
    mpf[b,]<- mn(xb[,b],length(XX_1))/n
    }

MPF<- apply(mpf,1,cumsum)
med<- apply(MPF,1,median)
LL<- apply(MPF,1,function(x)quantile(x,probs=0.025))
UU<- apply(MPF,1,function(x)quantile(x,probs=0.975))

#############################################################################################################################
#####~~~~~Section-2: Compute global annual averages above pre-industrial and above HadCRUT4 baseline i.e, 1961-1990~~~~~#####
#############################################################################################################################
ann_avg_lk<- read.table(paste("D:/files_18feb/ann_avg_2decdata.txt")) 		#167x10000
yr<- 1850:2016     										#years

ind_1850_1900<- which(yr>=1850 & yr<=1900)
CLM_O<- ann_avg_lk[ind_1850_1900,]
clm_o<- as.numeric(apply(CLM_O,2, mean))

ind_1961_1990<- which(yr>=1961 & yr<=1990)
CLM_HD<- ann_avg_lk[ind_1961_1990,] 
clm_hd<- apply(CLM_HD,2, mean)

####~~~~global annual averages above pre-industrial
yy3a<- XX_1[c(xb)]   			#median(yy3a): 0.076, quantile(yy3a, 0.025):-0.13, quantile(yy3a, 0.975):0.29


ANN_PREIND<- matrix(0, nrow=167, ncol=(n*10000))
		for(i in 1: 167){
            	YY2a<- as.numeric(ann_avg_lk[i,]-clm_o)
			yy3c<- rep(YY2a, each=n)
			ANN_PREIND[i,]<- yy3a+yy3c
		      }
ann_preind<- round(ANN_PREIND,2)
	
	##Median ts and 95% CI
	med_ann_preind<- apply(ann_preind,1,median)
	LL_ann_preind<- apply(ann_preind,1,function(x)quantile(x, prob=0.025))
	UL_ann_preind<- apply(ann_preind,1,function(x)quantile(x, prob=0.975))

####~~~~global annual averages above hadCRUT4 baseline i.e. 1961-1990
ANN_HDCLM<- ann_avg_lk- matrix(rep(clm_hd, each=nrow(ann_avg_lk)),ncol=10000)
ann_hdclm<- round(ANN_HDCLM,2)

	##Median ts and 95% confidence interval
	med_hdclm<- apply(ann_hdclm,1,median)
	LL_hdclm<- apply(ann_hdclm,1,function(x)quantile(x, prob=0.025))
	UL_hdclm<- apply(ann_hdclm,1,function(x)quantile(x, prob=0.975))

X<- c(1:nrow(ann_avg_lk))

X11()
e<- paste("C:/Users/maryam/Dropbox/PhD UCL/kriging_and_FEOFs/post_upgrade/paper1.5oC_r_files/figures/fig1.png",sep="")
png(e,width = 8.5, height = 6, units = 'in', res = 300)
	grid<- matrix(c(1,1,1,3,3, 
			    1,1,1,3,3,
                	    1,1,1,3,3,
		    	    1,1,1,3,3,
			    2,2,2,4,4,
		    	    2,2,2,4,4,
		    	    2,2,2,4,4,
		    	    2,2,2,4,4),nrow=8, ncol=5, byrow=TRUE)
layout(grid)
par(mar=c(4,4,2,0), oma=c(0.5,0.5,0,0))  

##panel(a)
	plot(X,med_hdclm, type='n',ylab=substitute(paste(~Delta,"T wrt 1961-90")),
	xlab="", ylim=c(-1, 2),xaxt="n",yaxt="n",main="")                          
	axis(side = 1,at = seq(1,length(X), by=10), labels = seq(1850,2016, by=10),tck= -.02)
	axis(side = 2,at = seq(-1,2,by=0.2))
	abline(h=0, lty=2, col='gray')

polygon(c(X,rev(X)),c(LL_hdclm,rev(UL_hdclm)),col = "grey75", border = FALSE)
lines(X,med_hdclm, col="black")
legend(-7,2.1, bty='n',col= c("black","grey75"),c("Median","Uncertainty"),lty=c(1,1,1), lwd=c(1,4,4))
mtext("a",3, line=0.3, adj=-0.13, cex=1.2)

##panel(b)
	 plot(X,med_ann_preind, type='n',ylab=substitute(paste(~Delta,"T wrt preindustrial")),
	 xlab="Time series of global annual mean temperatures", 
       ylim=c(-1, 2),xaxt="n",yaxt="n",main="")                          
	 axis(side = 1,at = seq(1,length(X), by=10), labels = seq(1850,2016, by=10),tck= -.02)
	 axis(side = 2,at = seq(-1,2,by=0.2))
	 abline(h=0, lty=2, col='gray')
	 abline(h=1.5, col='red')
	 abline(h=2, col='red')
polygon(c(X,rev(X)),c(LL_ann_preind,rev(UL_ann_preind)),col = "grey75", border = FALSE)
lines(X,med_ann_preind, col="black")
mtext("c",3, line=1.5, adj=-.13, cex=1.2)

##panel (c)
	#hist(yy3a,breaks=seq(-0.3,0.5, by=0.05), prob=TRUE,
	#main=substitute(paste("Preindustrial correction")),ylab="probability",  xlab="")

hist(yy3a, breaks=XX_1, prob=TRUE,xaxt="n",yaxt="n",xlim=range(-0.3,0.5),
main=substitute(paste("Preindustrial correction")),ylab="probability",  xlab="")
axis(1, at=  round(seq(-0.3,0.5, by=0.1),1))
axis(2, at= seq(0,3.5,by=0.5))
	#abline(v=median(yy3a), col=2)
mtext("b",3, line=0.3, adj=-0.17, cex=1.2)

#hist(as.numeric(ann_preind[167,]), prob=TRUE,breaks=seq(0.7,1.85,by=0.05),xaxt="n",
hist(as.numeric(ann_preind[167,]), prob=TRUE,breaks=seq(0.7,1.82,by=0.02),xaxt="n",yaxt="n",ylim=c(0,3.5),
col="blue", xlab=substitute(paste(~Delta,"T wrt preindustrial")), ylab="probability", 
main=substitute(paste(2016~Delta,"T wrt preindustrial")))
axis(1, at= seq(0.7,1.82, by=0.1))
axis(2, at= seq(0,3.5,by=0.5))
mtext("d",3, line=1.5, adj=-0.17, cex=1.2)
#dev.off()


##2016 estimates
avg_2016<- ann_avg_lk[167,]
YY2<- as.numeric(avg_2016-clm_o)    #T_2016 - T_[1850-1900]+T_c
yy3a<- XX_1[c(xb)]   			#median(yy3a)
			   			#yy3a<- rep(clm_1,10000)
yy3b<- rep(YY2,each=n)
			   			#yy3b<- rep(YY2, each=length(clm_1))
YY4<- yy3a+yy3b
YY4<- round(YY4,2)
median(YY4)
quantile(YY4, prob=0.025)
quantile(YY4, prob=0.975)
sum(YY4>=1.5)/length(YY4)
sum(YY4>=1.1)/length(YY4)


####supp_figure
X11()
par(mfrow=c(2,2))
hist(as.numeric(ann_hdclm[12,]),breaks=c(seq(-1.05,0.05,by=0.05)), prob=TRUE,
 main="1861 wrt 1961-1990", xlab=paste("L=",LL_hdclm[12],"M=",med_hdclm[12],"U=",round(UL_hdclm[12],2)))
abline(v=LL_hdclm[12],col=2)
abline(v=UL_hdclm[12],col=2)
abline(v=med_hdclm[12],col=3)

hist(as.numeric(ann_preind[12,]),breaks=c(seq(-0.7,0.7,by=0.05)), prob=TRUE, 
main="1861 wrt preindustrial", xlab=paste("L=",LL_ann_preind[12],"M=",med_ann_preind[12],"U=",round(UL_ann_preind[12],2)))
abline(v=LL_ann_preind[12],col=2)
abline(v=UL_ann_preind[12],col=2)
abline(v=med_ann_preind[12],col=3)

hist(as.numeric(ann_hdclm[167,]), prob=TRUE,#ylim=c(0,16),
main="2016 wrt 1961-1990", xlab=paste("L=",LL_hdclm[167],"M=",med_hdclm[167],"U=",round(UL_hdclm[167],2)))
abline(v=LL_hdclm[167],col=2)
abline(v=UL_hdclm[167],col=2)
abline(v=med_hdclm[167],col=3)

hist(as.numeric(ann_preind[167,]), prob=TRUE,#ylim=c(0,16),
main="2016 wrt preindustrial", xlab=paste("L=",LL_ann_preind[167],"M=",med_ann_preind[167],"U=",round(UL_ann_preind[167],2)))
abline(v=LL_ann_preind[167],col=2)
abline(v=UL_ann_preind[167],col=2)
abline(v=med_ann_preind[167],col=3)



#X11()
#plot(T[,2])
##it is read as for example F(-0.39)=P(X<= -0.39)= 0 & F(-0.09)=P(X<=-0.09)= 0.0006666667 so on
##so F(-0.38)-F(-0.39)=P(-0.39<=X<=-0.38)=0 & 
##F(-0.09)-F(-0.1)= P(-0.1<=X<=-0.09)= 0-0.0006666667=0.0006666667 so on
##https://stackoverflow.com/questions/37151571/calculate-derivative-of-cumulative-distribution-cdf-to-get-probability-density



























##A cumulative density function of the difference between a "pre-industrial" 1400-1800 baseline and 1850-1900.
T<- read.table("C:/Users/maryam/Dropbox/PhD UCL/kriging_and_FEOFs/post_upgrade/paper1.5oC_r_files/data_andrew/cdf_All_Forcings.txt")

		X11()
		plot(T[,2])

		##it is read as for example F(-0.39)=P(X<= -0.39)= 0 & F(-0.09)=P(X<=-0.09)= 0.0006666667 so on
		##so F(-0.38)-F(-0.39)=P(-0.39<=X<=-0.38)=0 & 
		##F(-0.09)-F(-0.1)= P(-0.1<=X<=-0.09)= 0-0.0006666667=0.0006666667 so on

	##A probability density function of the difference between a "pre-industrial" 1400-1800 baseline and 1850-1900.
	pdf_a<- diff(T[,2])
	freq<- 30000*pdf_a    #multiply with 10,002 instead of 10000 so that after rounding we have 10,000
		#sum(round(freq))
		#range(T[,1])

	##Frequency distribution
	breaks<- seq(-0.39,0.4, by=0.01)
	Cl<- cbind(seq(breaks[1],breaks[length(breaks)-1],by=0.01), seq(breaks[2],breaks[length(breaks)],by=0.01))
	XX<- apply(Cl,1,mean)
	YY<- cbind(XX,round(freq))
	
##samples of the difference of climatology 1400-1800 to 1850-1900 
set.seed(1002)
CLM_A<- rep(XX,round(freq))
clm_a<- sample(CLM_A,size=10000, replace=FALSE) 

		X11()
		par(mfrow=c(1,2))
		plot(T[,2])
		hist(clm_a)
		abline(v=median(clm_a), col="red")

CLM_O<- read.table("C:/Users/maryam/Dropbox/PhD UCL/kriging_and_FEOFs/post_upgrade/hd_lk_ts_data/1850_1900_clm.txt")
clm_o<- apply(CLM_O,2, mean)

AVG_2016<- read.table("C:/Users/maryam/Dropbox/PhD UCL/kriging_and_FEOFs/post_upgrade/hd_lk_ts_data/2016_mon_area_avg.txt")
avg_2016<- apply(AVG_2016,2, mean)


YY<- avg_2016-clm_o+sample(clm_a,10000,replace=TRUE)         ##used went sent before YY<- mean(avg_2016)-clm_o+clm_a

			median(YY)			#1.24874
			##5-95% range
			quantile(YY,0.05)		#1.07549
			quantile(YY,0.95)		#1.458502 

			#95% confidence interval
			quantile(YY,0.025)	#1.050423
			quantile(YY,0.975)	#1.477733

			#99% confidence interval
			quantile(YY, 0.005)     #0.9987
			quantile(YY,0.995)	#1.5623

			#empirical probability
			length(which(YY>1.5))/length(YY)

II<- seq(1:100)
ann_avg_lk<- read.table(paste("C:/Users/maryam/Dropbox/PhD UCL/kriging_and_FEOFs/post_upgrade/hd_lk_ts_data/ann_avg_lk_fr_hdens",II[1],".txt", sep=""))

for(i in 2:length(II)){
HH<- read.table(paste("C:/Users/maryam/Dropbox/PhD UCL/kriging_and_FEOFs/post_upgrade/hd_lk_ts_data/ann_avg_lk_fr_hdens",II[i],".txt", sep=""))
ann_avg_lk<- cbind(HH, ann_avg_lk)
}
med_ann_lk<- apply(ann_avg_lk,1,median)
X<- c(1:nrow(ann_avg_lk))

	#check
	#avg_15<- ann_avg_lk[15,]
	#YYY<- avg_15-clm_o+clm_a  

X11()
##e<- paste("C:/Users/maryam/Dropbox/PhD UCL/kriging_and_FEOFs/post_upgrade/paper1.5oC_r_files/figures/fig_1.png",sep="")
##png(e,width = 8.5, height = 6, units = 'in', res = 300)
	grid<- matrix(c(1,1,1,0,0,
                	    1,1,1,3,3,
		    	    1,1,1,3,3,
		    	    2,2,2,3,3,
		    	    2,2,2,3,3,
		    	    2,2,2,0,0),nrow=6, ncol=5, byrow=TRUE)
layout(grid)
par(mar=c(4,4,1,0), oma=c(0.5,0.5,0,0))  

##panel(a)
	CLM_HD<- read.table("C:/Users/maryam/Dropbox/PhD UCL/kriging_and_FEOFs/post_upgrade/hd_lk_ts_data/1961_1990_clm.txt")
	clm_hd<- apply(CLM_HD,2, mean)
	ann_hdclm<- ann_avg_lk- matrix(rep(clm_hd, each=nrow(ann_avg_lk)),ncol=10000)

	med_hdclm<- apply(ann_hdclm,1,median)
	LL_hdclm<- apply(ann_hdclm,1,function(x)quantile(x, prob=0.025))
	UL_hdclm<- apply(ann_hdclm,1,function(x)quantile(x, prob=0.975))

plot(X,med_hdclm, type='n',ylab=substitute(paste(~Delta,"T wrt 1961-90")),
	 xlab="", ylim=c(-1, 2),xaxt="n",yaxt="n")                          
	 axis(side = 1,at = seq(1,length(X), by=10), labels = seq(1850,2016, by=10),tck= -.02)
	 axis(side = 2,at = seq(-1,2,by=0.2))
	 abline(h=0, lty=2, col='gray')
	 #abline(h=1.5, col='red')
	 #abline(h=2, col='red')
polygon(c(X,rev(X)),c(LL_hdclm,rev(UL_hdclm)),col = "grey75", border = FALSE)
lines(X,med_hdclm, col="black")
legend(-7,2.1, bty='n',col= c("black","grey75"),c("Median","Uncertainty"),lty=c(1,1,1), lwd=c(1,4,4))

##panel(b)
	ann_preind<- ann_avg_lk- matrix(rep(clm_o, each=nrow(ann_avg_lk)),ncol=10000)+ matrix(rep(clm_a, each=nrow(ann_avg_lk)),ncol=10000)
	med_ann_preind<- apply(ann_preind,1,median)
	LL_ann_preind<- apply(ann_preind,1,function(x)quantile(x, prob=0.025))
	UL_ann_preind<- apply(ann_preind,1,function(x)quantile(x, prob=0.975))

plot(X,med_ann_preind, type='n',ylab=substitute(paste(~Delta,"T wrt preindustrial")),
	 xlab="Time series of global annual mean temperatures", 
       ylim=c(-1, 2),xaxt="n",yaxt="n")                          
	 axis(side = 1,at = seq(1,length(X), by=10), labels = seq(1850,2016, by=10),tck= -.02)
	 axis(side = 2,at = seq(-1,2,by=0.2))
	 abline(h=0, lty=2, col='gray')
	 abline(h=1.5, col='red')
	 abline(h=2, col='red')
polygon(c(X,rev(X)),c(LL_ann_preind,rev(UL_ann_preind)),col = "grey75", border = FALSE)
lines(X,med_ann_preind, col="black")

##panel (c)

##dev.off()










































##other options
##https://stackoverflow.com/questions/42941091/sample-from-custom-distribution-in-r
##http://pj.freefaculty.org/guides/stat/Distributions/DrawingRandomSamples/DrawingSamples.pdf


library(sROC)
CLM_O<- read.table("C:/Users/maryam/Dropbox/PhD UCL/kriging_and_FEOFs/post_upgrade/hd_lk_ts_data/1850_1900_clm.txt")
clm_o<- apply(CLM_O,2, mean)

AVG_2016<- read.table("C:/Users/maryam/Dropbox/PhD UCL/kriging_and_FEOFs/post_upgrade/hd_lk_ts_data/2016_mon_area_avg.txt")
avg_2016<- apply(AVG_2016,2, mean)

YY2<- avg_2016-clm_o

##~~CDF
##A cumulative density function of the difference between a "pre-industrial" 1400-1800 baseline and 1850-1900.
T<- read.table("C:/Users/maryam/Dropbox/PhD UCL/kriging_and_FEOFs/post_upgrade/paper1.5oC_r_files/data_andrew/cdf_All_Forcings.txt")

		X11()
		plot(T[,2])
		##it is read as for example F(-0.39)=P(X<= -0.39)= 0 & F(-0.09)=P(X<=-0.09)= 0.0006666667 so on
		##so F(-0.38)-F(-0.39)=P(-0.39<=X<=-0.38)=0 & 
		##F(-0.09)-F(-0.1)= P(-0.1<=X<=-0.09)= 0-0.0006666667=0.0006666667 so on

##~~Obtain samples having the CDF above
		##A probability density function of the difference between a "pre-industrial" 1400-1800 baseline and 1850-1900.
	pdf_a<- diff(T[,2])
	freq<- 30*pdf_a    #multiply with 10,002 instead of 10000 so that after rounding we have 10,000
		#sum(round(freq))
		#range(T[,1])

	##Frequency distribution
	breaks<- seq(min(T[,1]),max(T[,1]), by=0.01)
	Cl<- cbind(seq(breaks[1],breaks[length(breaks)-1],by=(T[,1][2]-T[,1][1])), seq(breaks[2],breaks[length(breaks)],by=(T[,1][2]-T[,1][1])))
	XX<- apply(Cl,1,mean)
	YY<- cbind(XX,round(freq))

##samples of the difference of climatology 1400-1800 to 1850-1900 
clm_a<- rep(XX,round(freq))   #median(clm_a)=0.075 as is the median of CDF above 

##Kernel density estimate and confidence interval of the empirical cdf above
x<- clm_a
x.CDF <- kCDF(x)
x.CDF
ll<- CI.CDF(x.CDF)$Fhat.lower
uu<- CI.CDF(x.CDF)$Fhat.upper
xhat<- x.CDF$x
fhat_uu<- sort(uu)   #  plot(sort(uu),type="n")
			   #  lines(sort(ll), col=1)
                     #  lines((ll),col=2)
			   #  lines(sort(uu), col=1)
			   #  lines(uu, col=2)
fhat_ll<- sort(ll)
fhat<- x.CDF$Fhat

		#X11()
		#par(mfrow=c(1,2))
		#plot(T[,2])
		#plot(x.CDF, alpha=0.05, main="Kernel estimate of distribution function")
	
	pdf_1<- diff(fhat)
	freq_1<- 60*pdf_1 

	breaks_1<- seq(min(xhat),max(xhat), by=c(xhat[2]-xhat[1]))
	Cl_1<- cbind(seq(breaks_1[1],breaks_1[length(breaks_1)-1],by=c(xhat[2]-xhat[1])), seq(breaks_1[2],breaks_1[length(breaks_1)],by=c(xhat[2]-xhat[1])))
	XX_1<- apply(Cl_1,1,mean)
	YY_1<- cbind(XX_1,round(freq_1))

##samples of the difference of climatology 1400-1800 to 1850-1900 
clm_1<- rep(XX_1,round(freq_1))

set.seed(102)
sam_sz<- length(clm_1) #100

		clm_sam_1<- sample(clm_1,sam_sz,replace=FALSE)
		clm_sam_full<- rep(c(clm_sam_1),10000)
		#clm_sam_full<- matrix(0, nrow=100,ncol=10000)
		#for(i in 1:10000){
		#clm_sam_full[,i]<- sample(clm_1,100,replace=FALSE)
		#}


YY3<- rep(YY2,each=sam_sz)
YY4<- YY3+clm_sam_full

median(YY4)
quantile(YY4, probs=0.025)
quantile(YY4, probs=0.975)

sum((YY4>=1.5))/length(YY4)

X11()
par(mfrow=c(2,2))
plot(T[,1],T[,2],col=1,xlab="",ylab="", pch=1)
lines(xhat,fhat,col=2)
legend("topleft", bty="n",c("Empirical CDF", "Kernel estimate"),col=c(1,2), pch=c(1,NA),lty=c(NA,1))
hist(clm_a,prob=TRUE, main="PDF from empirical CDF", xlab="median=0.075")
abline(v=median(clm_a),col=2)
hist(clm_1,prob=TRUE, main="PDF from kernel estimate",xlab="median=0.077")
abline(v=median(clm_1),col=2)
hist(YY4, freq=FALSE,breaks=seq(0.7,1.85,by=0.025),xaxt="n",
col="blue", xlab=substitute(paste(~Delta,"T wrt preindustrial")), ylab="probability", 
main=substitute(paste(2016~Delta,"T wrt preindustrial")))
axis(1, at= seq(0.7,1.85, by=0.1))

X11()
hist(clm_sam_1, freq=FALSE)


#median(YY4)
#1.25705
#quantile(YY4, probs=0.025)
#    2.5% 
#1.046185 
#quantile(YY4, probs=0.975)
#   97.5% 
#1.511734 
#sum((YY4>=1.5))/length(YY4)
#0.0308


https://cran.ms.unimelb.edu.au/web/packages/IPSUR/vignettes/IPSUR.pdf

library(logspline)