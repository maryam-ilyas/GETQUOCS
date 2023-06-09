rm(list=ls())
library(ncdf4)
library(fields)
library(LatticeKrig)
library(lhs)

##~~HadCRUT4 median 
source("C:/Users/maryam/Dropbox/Project_NCAR/abc_lk/abc_lk_functions.R")
##HadCRUT4 median ensemble member
      a1<- "D:/hd_sprs_dat/HadCRUT.4.5.0.0.median.nc"  
     	fid.t<- nc_open(a1,verbose=TRUE,write=FALSE)
	time.t<- ncvar_get(fid.t,'time')     
	origin.t <- fid.t$dim$time$units
	tunits.t <- strsplit(fid.t$dim$time$units,split=" ") 
	if(identical(tunits.t[[1]][1],"minutes")) tunits.t[[1]][1]<-"mins"
	TIME.t<- strptime(paste(tunits.t[[1]][3],tunits.t[[1]][4]),"%Y-%m-%d %H:%M:%S",tz="GMT")+
                           as.numeric(as.difftime(time.t,units=tunits.t[[1]][1]),units="secs")
				   ##1850-01-16 to 2015-05-16
	lat.t <- ncvar_get(fid.t,'latitude');lat.t                 		#length(lat.t)         #36
	lon.t <- ncvar_get(fid.t,'longitude'); lon.t               		#length(lon.t)         #72
	tempr.a <- ncvar_get(fid.t,'temperature_anomaly')          		#dim(tempr)            #72 36 1985

	##Spatial Coverage
	##Percentage of total available spatial points at each time point
	tt1<- rep(0,length(TIME.t))
	for(i in 1:length(TIME.t)){
	tt1[i]<-(sum(!is.na(tempr.a[,,i])))
	}
	tt<- (100* (tt1/(length(lat.t)*length(lon.t))))
      ##tt[which(round(tt)==50)]    #1297 (50.03858%) spatial points are at indices 986,1018, 1102
	##cbind( order(tt[which(round(tt)==50)]),sort(tt[which(round(tt)==50)]))

      #id<- c(which.min(tt1),986,which.max(tt1))
      #aval.n<- tt1[id]
      #aval.perc<- (100* (tt1/(length(lat.t)*length(lon.t))))[id] 
	id<- 137
            TIME.t[id]                          #"1861-05-16 12:00:00 GMT"
							#"1931-03-16 12:00:00 GMT"
                                          #"1988-02-15 12:00:00 GMT"
      date_min<- format(TIME.t[id[1]], "%B-%Y")    #"May-1861"             
	date_50<- format(TIME.t[id[2]], "%B-%Y")     #"February-1932" 
	date_max<- format(TIME.t[id[3]], "%B-%Y")    #"February-1988"

##50% coverage, February-1932
##data
Y<- c(tempr.a[,,id])
ind_s<- !is.na(Y)
y<- Y[ind_s]  

##locations
X<- as.matrix((expand.grid(lon.t,lat.t))) 
x<- X[ind_s,]

    ##~~~empirical variogram 
    ##d<- rdist.earth(x,x2,miles=FALSE)
    DD<- rdist.earth(x,miles=FALSE)[upper.tri(rdist.earth(x,miles=FALSE))];	round(max(DD)/2)
    D.max<- 6300 	
    l<- 13	

    vg<- vgm(x,y,D.max,l)

	#~~observed data and empirical variogram
	X11()
	#e<-  paste("C:/Users/maryam/Dropbox/Talks/upgrade_presentation/figs_ties/abc_1.png")
  	#png(e,width = 10, height = 5, units = 'in', res = 300)
    	par(oma=c(0,0.5,0,0),mar=c(2.5,2,2,1),mfrow=c(1,2))
	
	##observed data
	zr<- range(y, na.rm=TRUE);zr    
	brks = round(seq(-12,11, by=0.5),1);brks ;length(brks)
	cool = rainbow(sum(brks<0), start=rgb2hsv(col2rgb('cyan'))[1], end=rgb2hsv(col2rgb('blue'))[1])
	warm = rainbow(sum(brks>0), start=rgb2hsv(col2rgb('red'))[1], end=rgb2hsv(col2rgb('yellow'))[1])
	lut = c(rev(cool),rev(warm))
	NL=length(brks)-1
	mapping = list(x=lon.t, y=lat.t,z = matrix(Y,nrow=length(lon.t),ncol=length(lat.t)))
	image.plot(mapping, col=lut, main= paste("Data (",date_min,"), n=",length(y) ), 
	zlim=zr,xlim=range(lon.t), ylim=range(lat.t),ylab="",xlab="",horizontal=TRUE,breaks=brks, lab.breaks=brks)
	
	##Empirical variogram
	plot(vg[,1], vg[,3], xlab="distance (in km)", ylab="",main="Empirical variogram", 
	ylim=c(0,5),xlim=range(vg[,1]), col=2, pch=2) 		
	#dev.off()

xr=cbind( c(-180, 180), c(-90,90))
nc=3
nl=3
ALPHA<- c(0.2450707, 0.01606482, 0.7388645)
n<- nrow(x)
N<- nrow(X)
M<- 1

##1) simulate from the prior 
I<- 250     	#no. of iterations 

set.seed(1224)
##Draw a Latin Hypercube sample from uniform distribution for lambda and a.wght
##l_aw<- 1.0001 ; u_aw<- 3; l_l<- 0.001; u_l<- 4

LKinfo1<- LKrigSetup( xr, LKGeometry="LKSphere",
                      NC=nc, a.wght=1.1,lambda=0.6,
                      nlevel=nl, alpha=ALPHA, choleskyMemory=list(nnzR= 4E6))

l.omg<- -4.5	; omega2Awght(l.omg, LKinfo1)
u.omg<- 0.55	; omega2Awght(u.omg, LKinfo1)

l.lnl<- -6.9	; exp(l.lnl)  
u.lnl<- 1.4       ; exp(u.lnl)  0

L<- randomLHS(I,2)       #draw a Latin hypercube first from a set of standard uniform distribution
LL<- cbind(qunif(L[,1],l.lnl,u.lnl), qunif(L[,2], l.omg,u.omg))  #transform the sample from U(0,1) to U(0.01,6) and U(0.01,16)

theta_l<- exp(LL[,1])     #for lambda
theta_aw<- omega2Awght(LL[,2], LKinfo1)
summary(theta_l)
summary(theta_aw)

		#theta_aw<- omega2Awght(LL[,2], LKinfo1)      	#for a.wght	
		#not converting awght to omega as it does not work well for very small values of a.wght
		#and covariance changes drastically over small values of a.wght
		#check for a.wght=1.0001, Awght2Omega(1.0001,LKinfo1)=-4.60517 whereas omega2Awght(-4.60517, LKinfo1)=1.02
		#omega2Awght(-10, LKinfo1)=1.000091
      	#identify reason later????

#one sample for each parametric value
#system.time({
	rss<- rep(0,I)
	for(i in 1:I){   			   		#3 sec for one iteration
	LKinfo<- LKrigSetup( xr,LKGeometry="LKSphere",
                      NC=nc, a.wght=theta_aw[i],
                      nlevel=nl, alpha=ALPHA, lambda=(theta_l[i]),
                      choleskyMemory=list(nnzR= 4E6))                       
	ystar<- LKrig.sim(x, LKinfo=LKinfo, M=M) + sqrt(theta_l[i])*matrix(rnorm(n*M),n,M)	
	vg.t<- vgm(x,ystar, D.max,l)
	rss[i]<- sum((vg.t[,3]-vg[,3])^2)  #residual sums of square	
	#rss[i]<- sum(( (vg.t[,3]-vg[,3])/ sd(vg.t[,3]) )^2)	#weighted distance function https://projecteuclid.org/download/pdfview_1/euclid.ba/1460641065
										#weighted distance function would make sense of we consider more than one replication as
										#with that we would be able to estimate the sd of summary
										#so consider sqrt below in ref to Euclidean distance
	}
#})
ind<- (order(sqrt(rss)))[1:10]	
(length(ind)/I)*100                      #epsilon=3.5, 10 (without sqrt)
(sqrt(rss[ind]))[10]			     #epsilon=1.8,10 (with sqrt)
						     ##ep=1.43,
# user  system elapsed 
# 473.25   10.17  500.15 

aw_abc<- theta_aw[ind] 
lam_abc<- theta_l[ind]

mean(lam_abc); sd(lam_abc)
mean(aw_abc); sd(aw_abc) 
summary(lam_abc)
summary(aw_abc)

 #Figure-3  
X11()
e1<- paste("C:/Users/maryam/Dropbox/AMT_paper/1861_may_abc_post_r.png")
png(e1,width = 10, height = 5, units = 'in', res = 300)
par(oma=c(0,1,0,0),mar=c(2.5,2,2,1),mfrow=c(1,2))

 hist(lam_abc, freq=FALSE,breaks=seq(0.4,1.8, by=0.22), main=expression(paste("(a)", ~~lambda)), xlab="", xlim=c(-3,6), ylim=c(0,3))
 abline(v=mean(lam_abc), lty=2)
 curve(dunif(x,min(theta_l),max(theta_l)), col="darkgrey", add=TRUE) 
 legend("topleft", bty="n",c("ABC posterior","posterior mean", "prior"),lty=c(1,2,1), col=c(1,1,"darkgrey"))

 hist(aw_abc, breaks=seq(1,3.5, by=0.5),main=expression(paste("(b)", ~~a.wght)), xlab="",xlim=c(-3,6), ylim=c(0,4))
 abline(v=mean(aw_abc), lty=2)
 curve(dunif(x,min(theta_aw-0.05),max(theta_aw)), col="darkgrey", add=TRUE) 
 dev.off()

#50% coverage
#id<- 986
#format(TIME.t[986], "%B-%Y") 	#"February-1932"
#tt[986]#50.03858

#max coverage
#id<- 1658
#format(TIME.t[1658], "%B-%Y") 	#"February-1988"
#tt[1658]  #77.77%

filename=paste("D:/parameters10from250.nc",sep="")
		ncnew<- nc_open(filename, write=TRUE) 
		AA<- ncvar_get(ncnew,"parameters")
      	dim(AA)
aw_abcc<- AA[,,986][,2][-11]  	#11th value is tolerance threshold
lam_abcc<- AA[,,986][,1][-11]	#11th value is tolerance threshold

aw_abccc<- AA[,,1658][,2][-11]  	#11th value is tolerance threshold
lam_abccc<- AA[,,1658][,1][-11]	#11th value is tolerance threshold

#Priors and Posteriors
X11()
e1<- paste("C:/Users/maryam/Dropbox/AMT_paper/1932_feb_abc_post_r.png") 
png(e1,width = 10, height = 5, units = 'in', res = 300)
par(oma=c(0,1,0,0),mar=c(2.5,2,2,1),mfrow=c(1,2))

 #breaks=seq(0.4,1.8, by=0.22),
 hist(lam_abcc, freq=FALSE, main=expression(paste("(a)", ~~lambda)), xlab="", xlim=c(-3,6), ylim=c(0,3))
 abline(v=mean(lam_abcc), lty=2)
 curve(dunif(x,min(theta_l),max(theta_l)), col="darkgrey", add=TRUE) 
 legend("topleft", bty="n",c("ABC posterior","posterior mean", "prior"),lty=c(1,2,1), col=c(1,1,"darkgrey"))

 hist(aw_abcc,breaks=seq(1,2.6, length=5),main=expression(paste("(b)", ~~a.wght)), xlab="",xlim=c(-3,6), ylim=c(0,9))
 abline(v=mean(aw_abcc), lty=2)
 curve(dunif(x,min(theta_aw-0.05),max(theta_aw)), col="darkgrey", add=TRUE) 
 dev.off()

#Priors and Posteriors
X11()
e1<- paste("C:/Users/maryam/Dropbox/AMT_paper/1988_feb_abc_post_r.png") 
png(e1,width = 10, height = 5, units = 'in', res = 300)
par(oma=c(0,1,0,0),mar=c(2.5,2,2,1),mfrow=c(1,2))

 #breaks=seq(0.4,1.8, by=0.22),
 hist(lam_abccc, freq=FALSE, main=expression(paste("(a)", ~~lambda)), xlab="", xlim=c(-3,6), ylim=c(0,3))
 abline(v=mean(lam_abccc), lty=2)
 curve(dunif(x,min(theta_l),max(theta_l)), col="darkgrey", add=TRUE) 
 legend("topleft", bty="n",c("ABC posterior","posterior mean", "prior"),lty=c(1,2,1), col=c(1,1,"darkgrey"))
 
 #breaks=seq(1,3.5, by=0.5)
 hist(aw_abccc,main=expression(paste("(b)", ~~a.wght)), xlab="",xlim=c(-3,6), ylim=c(0,9))
 abline(v=mean(aw_abccc), lty=2)
 curve(dunif(x,min(theta_aw-0.05),max(theta_aw)), col="darkgrey", add=TRUE) 
 dev.off()

###########################################################################################################################################


#Figure-3  
#Priors and Posteriors
X11()
#e1<- paste("C:/Users/maryam/Dropbox/thesis/latex_thesis/figures/1861_may_abc_post_r.png")
e1<- paste("C:/Users/maryam/Dropbox/AMT_paper/1932_feb_abc_post_r.png") 
png(e1,width = 10, height = 5, units = 'in', res = 300)
par(oma=c(0,1,0,0),mar=c(2.5,2,2,1),mfrow=c(1,2))

 hist(lam_abc, freq=FALSE,breaks=seq(0.4,1.8, by=0.22), main=expression(paste("(a)", ~~lambda)), xlab="", xlim=c(-3,6), ylim=c(0,3))
 abline(v=mean(lam_abc), lty=2)
 curve(dunif(x,min(theta_l),max(theta_l)), col="darkgrey", add=TRUE) 
 legend("topleft", bty="n",c("ABC posterior","posterior mean", "prior"),lty=c(1,2,1), col=c(1,1,"darkgrey"))

 hist(aw_abc, breaks=seq(1,3.5, by=0.5),main=expression(paste("(b)", ~~a.wght)), xlab="",xlim=c(-3,6), ylim=c(0,4))
 abline(v=mean(aw_abc), lty=2)
 curve(dunif(x,min(theta_aw-0.05),max(theta_aw)), col="darkgrey", add=TRUE) 
 dev.off()



 #Figure-3  
 #Priors and Posteriors
 X11()
 #e1<- paste("C:/Users/maryam/Dropbox/thesis/latex_thesis/figures/1861_may_abc_post.png")
 #png(e1,width = 10, height = 5, units = 'in', res = 300)
 par(oma=c(0,1,0,0),mar=c(2.5,2,2,1),mfrow=c(1,2))
 d_a<- density(lam_abc)
 plot(d_a, main=expression(paste("(a)", ~~lambda)), xlab="", xlim=c(-3,6), ylim=c(0,1.4))
 abline(v=mean(lam_abc), lty=2)
 curve(dunif(x,min(theta_l),max(theta_l)), col="darkgrey", add=TRUE) 
 legend("topleft", bty="n",c("ABC posterior","posterior mean", "prior"),lty=c(1,2,1), col=c(1,1,"darkgrey"))

 d_b<- density(aw_abc, adjust=1)
 plot(d_b, main=expression(paste("(b)", ~~aw)), xlab="",xlim=c(-3,6), ylim=c(0,1.4))
 abline(v=mean(aw_abc), lty=2)
 curve(dunif(x,min(theta_aw),max(theta_aw)), col="darkgrey", add=TRUE) 
 #dev.off()

###Approximate Bayesian Inference
###Consider likely values of lambda and a.wght 
theta_l_hat<- theta_l[ind]     #likely values of lambda
theta_aw_hat<- theta_aw[ind]   #likely values of a.wght

Y_hat<- matrix(0,ncol=length(theta_l_hat),nrow=N)
Y_hat_se<- matrix(0,ncol=length(theta_l_hat),nrow=N)

for(k in 1:length(theta_l_hat)){
LKinfo3<- LKrigSetup( xr,LKGeometry="LKSphere",
                      NC=nc, a.wght=theta_aw_hat[k], lambda=theta_l_hat[k],
                      nlevel=nl, alpha=ALPHA,choleskyMemory=list(nnzR= 4E6))
obj<- LKrig(x,y,LKinfo=LKinfo3)
Y_hat[,k]<- predict(obj, xnew=X)
Y_hat_se[,k]<- predictSE(obj, xnew=X)
}

yhat<- apply(Y_hat,1,mean) 
sd_yhat<- apply(Y_hat_se,1,mean) 


##~~1~~centres of MR basis on sphere 
LKinfo11<- LKrigSetup( xr,
                      NC=nc, LKGeometry="LKSphere", a.wght=1.4,
                      nlevel=nl, alpha=ALPHA,
                      choleskyMemory=list(nnzR= 4E6))

hold0<- LKrigLatticeCenters( LKinfo11, Level=1)        #organize grid of lattice centres (162)
hold_1<- LKrigLatticeCenters( LKinfo11, Level=2)       #organize grid of lattice centres (642)
hold_2<- LKrigLatticeCenters( LKinfo11, Level=3)       #organize grid of lattice centres (2562)

hold<- LKrig.basis( X, LKinfo11)      #2592 x nb (3366=162+642+2562).....PHI
                                      #values of nb basis functions at 2592 locations
xg<- make.surface.grid( list( abcissa= lon.t, ordinate= lat.t) )

range(yhat)
range(sd_yhat)

##Extra calculations
Y.na<- tempr.a[,,id[1]]      
          	Y3<- as.vector(Y.na)  
		good3 <- !is.na(Y3)
g<- (yhat[good3]-y)/sd_yhat[good3]
m<- mean(g)
std<- sqrt(var(g))

############################################################
##spatial predictions of Ilyas et al data and ABC ensemble##
############################################################

SS<- read.table("C:/files_moved/hadcrut_v4_5_0_0_gapfilled_data/lk_hd_median.txt")
id<- 137
yhat_old<- as.numeric(SS[id[1],])


	X11()
	par(op) 		
	e<- paste("C:/Users/maryam/Dropbox/thesis/latex_thesis/figures/spatial_estimate_abc_may_1861_new_2.png")
      png(e,width = 13, height = 13, units = 'in', res = 300)
	par(oma=c(1,1,1,0),mar=c(7,3,3,2),mfrow=c(2,2))

	zr<- range(Y,yhat, na.rm=TRUE);zr  
	brks = seq(-4.4,6.6, by=0.4);brks ;length(brks) 	 #make sure the sequence has zero
      pal_full<- length(brks)+ ((sum(brks>0)-(sum(brks<0)))) #full palette
	pal_ret<- 1:(length(brks)-1)                           #retained palatte
	lut = rev(colorRampPalette(brewer.pal(11,"RdBu"))(pal_full)[pal_ret]) ##c(rev(cool),rev(warm))
	NL=length(brks)-1
      	
##(a)
  	mapping = list(x=lon.t, y=lat.t,z = matrix(yhat,nrow=length(lon.t),ncol=length(lat.t)))
      image.plot(mapping, nlevel=NL, col=lut, horizontal=TRUE,yaxt="n", xaxt="n", main="(a)",cex.main=1.5, 
	zlim=zr,xlim=range(lon.t), ylim=range(lat.t),ylab="", xlab="",breaks=brks, lab.breaks=brks,axis.args = list(cex.axis = 1.5))
	 xxx3<- X[!is.na(tempr.a[,,id[1]] ),]
	points(xxx3,cex=1,  pch=4, col="purple")
	axis(1, at=seq(-150,150, by=30),cex.axis=1.5)
	axis(2, at=seq(-90,90,30) ,cex.axis=1.5)
 
      world(add=TRUE)
##(b)
  	mapping = list(x=lon.t, y=lat.t,z = matrix(yhat_old,nrow=length(lon.t),ncol=length(lat.t)))
      image.plot(mapping, nlevel=NL, col=lut, horizontal=TRUE,yaxt="n", xaxt="n", main="(b)", cex.main=1.5,
	zlim=zr,xlim=range(lon.t), ylim=range(lat.t),ylab="", xlab="",breaks=brks, lab.breaks=brks, axis.args = list(cex.axis = 1.5))
	 xxx3<- X[!is.na(tempr.a[,,id[1]] ),]
	points(xxx3,cex=1,  pch=4, col="purple")
	axis(1, at=seq(-150,150, by=30), cex.axis=1.5)
	axis(2, at=seq(-90,90,30), cex.axis=1.5)  
      world(add=TRUE)

##(a)-(b)
	#zr2<- range(yhat-yhat_old);zr2
	#brks2 =  seq(-2.4,1.4,by=0.06);brks2 ;length(brks2)
	#warm2 = rainbow(length(brks2)-1, start=rgb2hsv(col2rgb('red'))[1], end=rgb2hsv(col2rgb('yellow'))[1])
	#lut2 = c(rev(warm2))
	#NL2=length(lut2)

	zr2<- range(yhat-yhat_old);zr2 
	brks2 = seq(-2.3,1.3, by=0.1);brks2 ;length(brks2) 	     #make sure the sequence has zero
      pal_full2<- length(brks2)+ ((sum(brks2<0)-(sum(brks2>0)))) #full palette
	pal_ret2<- 1:(length(brks2)-1)                             #retained palatte
	lut2 = rev(colorRampPalette(brewer.pal(11,"RdBu"))(pal_full2))[pal_ret2] ##c(rev(cool),rev(warm))
	NL2=length(brks2)-1 

	mapping = list(x=lon.t, y=lat.t,z = matrix(yhat-yhat_old,nrow=length(lon.t),ncol=length(lat.t)))
	image.plot(mapping, nlevel=NL2, col=lut2, horizontal=TRUE,yaxt="n", xaxt="n", main="(c)",cex.main=1.5, 
	zlim=zr2,xlim=range(lon.t), ylim=range(lat.t),ylab="", xlab="",breaks=brks2, lab.breaks=brks2,  axis.args = list(cex.axis = 1.5))
 		xx3<- X[!is.na(tempr.a[,,id[1]] ),]
      points(xx3,cex=1,  pch=4, col="purple")
	axis(1, at=seq(-150,150, by=30), cex.axis=1.5)
	axis(2, at=seq(-90,90,30), cex.axis=1.5) 
      world(add=TRUE)
	dev.off()

############################################
##SEs of Ilyas et al data and ABC ensemble##
############################################
SSS<- read.table("C:/files_moved/hadcrut_v4_5_0_0_gapfilled_data/lk_se_hd_median.txt")
SE.MIN<- as.numeric(SSS[id[1],])

	X11()
	par(op) 		
	e<- paste("C:/Users/maryam/Dropbox/thesis/latex_thesis/figures/spatial_estimate_abc_may_1861_new_1.png")
  	png(e,width = 13, height = 13, units = 'in', res = 300)
	par(oma=c(1,1,1,0),mar=c(7,3,3,2),mfrow=c(2,2))

	zr1<- range(sd_yhat,SE.MIN);zr1
	brks1 =  seq(0.26,1.32,by=0.02);brks1 ;length(brks1)
	warm1 = rainbow(length(brks1)-1, start=rgb2hsv(col2rgb('red'))[1], end=rgb2hsv(col2rgb('yellow'))[1])
	lut1 = c(rev(warm1))
	NL1=length(lut1) 
##(a)
  	mapping = list(x=lon.t, y=lat.t,z = matrix(sd_yhat,nrow=length(lon.t),ncol=length(lat.t)))
	image.plot(mapping, nlevel=NL1, col=lut1, horizontal=TRUE,yaxt="n", xaxt="n", main="(a)",cex.main=1.5, 
	zlim=zr1,xlim=range(lon.t), ylim=range(lat.t),ylab="", xlab="",breaks=brks1, lab.breaks=brks1,axis.args = list(cex.axis = 1.5))
 		xx3<- X[!is.na(tempr.a[,,id[1]] ),]
      points(xx3,cex=1,  pch=4, col="purple")
	axis(1, at=seq(-150,150, by=30), cex.axis=1.5)
	axis(2, at=seq(-90,90,30), cex.axis=1.5) 
      world(add=TRUE)
##(b)
  	mapping = list(x=lon.t, y=lat.t,z = matrix(SE.MIN,nrow=length(lon.t),ncol=length(lat.t)))
	image.plot(mapping, nlevel=NL1, col=lut1, horizontal=TRUE,yaxt="n", xaxt="n", main="(b)", 
	zlim=zr1,xlim=range(lon.t), ylim=range(lat.t),ylab="", xlab="",breaks=brks1, lab.breaks=brks1, cex.main=1.5, axis.args = list(cex.axis = 1.5))
 	xx3<- X[!is.na(tempr.a[,,id[1]] ),]
      points(xx3,cex=1,  pch=4, col="purple")
	axis(1, at=seq(-150,150, by=30), cex.axis=1.5)
	axis(2, at=seq(-90,90,30), cex.axis=1.5) 
      world(add=TRUE)
##(a)-(b)
	zr2<- range(sd_yhat-SE.MIN);zr2
	brks2 =  seq(-0.52,0.08,by=0.02);brks2 ;length(brks2)
	warm2 = rainbow(length(brks2)-1, start=rgb2hsv(col2rgb('red'))[1], end=rgb2hsv(col2rgb('yellow'))[1])
	lut2 = c(rev(warm2))
	NL2=length(lut2) 

	mapping = list(x=lon.t, y=lat.t,z = matrix(sd_yhat-SE.MIN,nrow=length(lon.t),ncol=length(lat.t)))
	image.plot(mapping, nlevel=NL2, col=lut2, horizontal=TRUE,yaxt="n", xaxt="n", main="(c)", 
	zlim=zr2,xlim=range(lon.t), ylim=range(lat.t),ylab="", xlab="",breaks=brks2, lab.breaks=brks2,cex.main=1.5, axis.args = list(cex.axis = 1.5))
 		xx3<- X[!is.na(tempr.a[,,id[1]] ),]
      points(xx3,cex=1,  pch=4, col="purple")
	axis(1, at=seq(-150,150, by=30), cex.axis=1.5)
	axis(2, at=seq(-90,90,30), cex.axis=1.5) 
      world(add=TRUE)
	dev.off()





