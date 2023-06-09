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
u.lnl<- 1.4       ; exp(u.lnl)  

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

	X11()
	par(op) 		
	e<- paste("C:/Users/maryam/Dropbox/thesis/latex_thesis/figures/spatial_estimate_abc_may_1861.png")
      png(e,width = 11, height = 10, units = 'in', res = 300)
	par(oma=c(0.7,0.5,0,0),mar=c(5,2,2,2),mfrow=c(3,2))
	zr<- range(Y,yhat, na.rm=TRUE);zr  
	
	brks = seq(-4.4,6.6, by=0.4);brks ;length(brks) 	 #make sure the sequence has zero
			
      pal_full<- length(brks)+ ((sum(brks>0)-(sum(brks<0)))) #full palette
	pal_ret<- 1:(length(brks)-1)                           #retained palatte
	lut = rev(colorRampPalette(brewer.pal(11,"RdBu"))(pal_full)[pal_ret]) ##c(rev(cool),rev(warm))
	NL=length(brks)-1
      	mapping = list(x=lon.t, y=lat.t,z = matrix(Y,nrow=length(lon.t),ncol=length(lat.t)))
		image(mapping, col=lut, main=paste("(a)"),xaxt="n", yaxt="n",
		zlim=zr,xlim=range(lon.t), ylim=range(lat.t),ylab="",xlab="")
	
	#set background color
	rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "lightgray")
		image(mapping, col=lut, main=paste("(a)"),xaxt="n", yaxt="n",
		zlim=zr,xlim=range(lon.t), ylim=range(lat.t),ylab="",xlab="",add=TRUE)
		axis(1, at=seq(-150,150, by=30))
		axis(2, at=seq(-90,90,30))  
      	world(add=TRUE)

 	plot( hold0,main="(b)",			
      	xlim=c(-167,168), ylim=c(-84,84),xaxt="n", yaxt="n",
      	ylab="", xlab="", pch=19, cex=0.4, type="n")
			axis(1, at=seq(-150,150, by=30))
			axis(2, at=seq(-60,60,30))   
	points(hold_2, col="pink", pch=19, cex=0.4)
	points(hold_1, col="cyan", pch=19, cex=0.4)
	points(hold0, col="black", pch=19, cex=0.4)
	contour( as.surface(xg, hold[,1]), add=TRUE, col="black")  #First basis function
      xxx3<- X[!is.na(tempr.a[,,id[1]] ),]
	points(xxx3, col="purple",cex=0.5,  pch=4)

 	mapping = list(x=lon.t, y=lat.t,z = matrix(yhat,nrow=length(lon.t),ncol=length(lat.t)))
      image.plot(mapping, nlevel=NL, col=lut, horizontal=TRUE,yaxt="n", xaxt="n", main="(c)", #main="(c) Multi-resolution lattice kriging",
	zlim=zr,xlim=range(lon.t), ylim=range(lat.t),ylab="", xlab="",breaks=brks, lab.breaks=brks)
	points(xxx3,cex=0.5,  pch=4, col="purple")
	axis(1, at=seq(-150,150, by=30))
	axis(2, at=seq(-90,90,30))  
      world(add=TRUE)

	zr1<- range(sd_yhat);zr1
	brks1 =  seq(0.26,1.32,by=0.02);brks1 ;length(brks1)
	warm1 = rainbow(length(brks1)-1, start=rgb2hsv(col2rgb('red'))[1], end=rgb2hsv(col2rgb('yellow'))[1])
	lut1 = c(rev(warm1))
	NL1=length(lut1) 

  	mapping = list(x=lon.t, y=lat.t,z = matrix(sd_yhat,nrow=length(lon.t),ncol=length(lat.t)))
	image.plot(mapping, nlevel=NL1, col=lut1, horizontal=TRUE,yaxt="n", xaxt="n", main="(d)", 
	zlim=zr1,xlim=range(lon.t), ylim=range(lat.t),ylab="", xlab="",breaks=brks1, lab.breaks=brks1)
 		xx3<- X[!is.na(tempr.a[,,id[1]] ),]
      points(xx3,cex=0.5,  pch=4, col="purple")
	axis(1, at=seq(-150,150, by=30))
	axis(2, at=seq(-90,90,30)) 
      world(add=TRUE)

	HH<- hist(g, prob=TRUE, xlab="(predictions-observed)/standard error", main="(e)", breaks=seq(-10,6, by=1))
	curve(dnorm(x, mean=m, sd=std), col="blue", lwd=1.4, add=TRUE, yaxt="n")

	#read SE presented in paper ilyas et al 2017
	SSS<- read.table("C:/files_moved/hadcrut_v4_5_0_0_gapfilled_data/lk_se_hd_median.txt")
	SE.MIN<- as.numeric(SSS[id[1],])

	zr1<- range(sd_yhat-SE.MIN);zr1
	brks1 =  seq(-0.52,0.08,by=0.02);brks1 ;length(brks1)
	warm1 = rainbow(length(brks1)-1, start=rgb2hsv(col2rgb('red'))[1], end=rgb2hsv(col2rgb('yellow'))[1])
	lut1 = c(rev(warm1))
	NL1=length(lut1) 

  	mapping = list(x=lon.t, y=lat.t,z = matrix(sd_yhat-SE.MIN,nrow=length(lon.t),ncol=length(lat.t)))
	image.plot(mapping, nlevel=NL1, col=lut1, horizontal=TRUE,yaxt="n", xaxt="n", main="(f)", 
	zlim=zr1,xlim=range(lon.t), ylim=range(lat.t),ylab="", xlab="",breaks=brks1, lab.breaks=brks1)
 		xx3<- X[!is.na(tempr.a[,,id[1]] ),]
      points(xx3,cex=0.5,  pch=4, col="purple")
	axis(1, at=seq(-150,150, by=30))
	axis(2, at=seq(-90,90,30)) 
      world(add=TRUE)
	dev.off()






















> cbind(lam_abc, aw_abc)
        lam_abc   aw_abc
 [1,] 0.8322848 3.071688
 [2,] 0.6491180 1.604935
 [3,] 0.5843478 2.820865
 [4,] 1.1665418 1.003463
 [5,] 0.6271231 2.083420
 [6,] 1.3257647 1.006819
 [7,] 0.7899294 2.213436
 [8,] 0.7129534 1.063388
 [9,] 1.0836773 1.742086
[10,] 1.6707843 1.057278
##Profile MLE
#system.time({
#LKinfo.search<- LKrigSetup(xr,LKGeometry="LKSphere", NC=nc, nlevel=nl, a.wght=1.01, alpha=ALPHA)
#search.results<- LKrigFindLambdaAwght(x,y,LKinfo=LKinfo.search)
#search.results$summary
#})
#sort(round(search.results$lnLike.eval[,1],4))
#sort(round(search.results$lnLike.eval[,3],4))
#lam.ml<- 0.1462012 #search.results$summary[6]    #0.1462012
#awght.ml<- 1.1387  #search.results$summary[7]    #1.1387 







id=138
> cbind(lam_abc, aw_abc)
        lam_abc   aw_abc
 [1,] 0.4218572 3.828947
 [2,] 0.5284317 1.133848
 [3,] 0.8487807 1.084028
 [4,] 0.7730788 1.259355
 [5,] 0.4455003 1.880290
 [6,] 0.2475698 3.548007
 [7,] 0.6335096 1.167367
 [8,] 0.7201140 1.057440
 [9,] 1.1323322 1.017218
[10,] 0.9668538 1.055551

> id
[1] 139
> cbind(lam_abc, aw_abc)
        lam_abc   aw_abc
 [1,] 0.3030249 2.288415
 [2,] 0.3207931 1.549745
 [3,] 0.6271231 2.083420
 [4,] 0.7129534 1.063388
 [5,] 0.7670692 1.049362
 [6,] 0.4847806 1.802487
 [7,] 0.5843478 2.820865
 [8,] 0.3739696 3.722133
 [9,] 0.3120511 1.363045
[10,] 0.9671312 1.081973

> id
[1] 140
> cbind(lam_abc, aw_abc)
          lam_abc   aw_abc
 [1,] 0.320793101 1.549745
 [2,] 0.021674544 2.630716
 [3,] 0.064799016 3.829781
 [4,] 0.199239799 1.175774
 [5,] 0.176228056 1.673479
 [6,] 0.008874635 2.258914
 [7,] 0.010122625 1.700338
 [8,] 0.908423862 1.000350
 [9,] 0.001523823 3.171011
[10,] 0.164999034 3.605068

> id
[1] 133
> cbind(lam_abc, aw_abc)
       lam_abc   aw_abc
 [1,] 2.192790 1.048626
 [2,] 2.089398 1.104972
 [3,] 2.517501 1.027476
 [4,] 2.854751 1.112348
 [5,] 1.816702 2.516519
 [6,] 2.718094 1.042475
 [7,] 2.313164 1.316367
 [8,] 1.083677 1.742086
 [9,] 2.277785 1.002048
[10,] 2.454556 2.101642




> lam_abc
  [1] 0.599964832 0.553912217 0.824500025 0.641751927 0.773728438 0.288202100 0.940264380 0.929719646 0.647552050 0.707014814 0.723780872 0.300263754
 [13] 0.803017210 0.881113613 0.765832647 0.655747843 0.714540588 1.102393230 0.956295061 0.588288102 0.819171398 0.172483445 1.472666484 0.950890489
 [25] 1.492813375 1.218613214 0.742084505 0.637112865 1.527893141 0.839858146 1.135526462 0.281089283 1.588047677 0.404441802 0.458153296 1.124979667
 [37] 0.263527126 0.610349091 1.208969686 1.018161427 0.889672731 1.409069155 0.601803870 1.111720894 0.632226050 0.581600845 1.952862801 0.726385248
 [49] 1.762511187 0.241750843 1.083649704 1.434935533 1.687036960 0.815689426 1.331723525 1.271123044 0.437689927 0.986630020 0.471465143 1.503932141
 [61] 1.330808940 1.512191717 0.361484689 1.455030838 0.061058935 0.375491189 0.620459939 0.191738154 0.679438749 1.788753854 1.834262371 0.749078273
 [73] 0.007240987 1.697822065 1.198162174 0.615397095 0.577712977 0.243859977 2.100637911 1.316067216 1.291356375 0.894951842 0.970862853 1.537289585
 [85] 1.610223509 1.239268379 1.384306663 1.915729943 1.934734806 1.553062221 0.017274708 0.479222970 0.522473716 0.207007722 0.029108726 0.995251495
 [97] 0.220079833 1.031753385 1.412703235 1.432191457

> aw_abc
  [1] 1.756137 1.083761 1.551016 1.165750 1.222255 2.058141 1.051030 1.293353 1.858923 1.460654 2.380847 1.437923 1.580634 1.236520
 [15] 1.262078 1.118375 2.290626 1.158431 1.402468 1.214204 1.808650 2.141360 1.000773 1.087617 1.008800 1.059876 1.115788 2.810840
 [29] 1.002894 1.699078 1.458061 2.544297 1.000524 1.515317 1.526950 1.003705 1.337502 1.117863 1.003430 1.033091 1.663303 1.008367
 [43] 2.080385 1.011722 1.299523 1.304034 1.005134 1.846491 1.000684 1.624592 1.546414 1.054718 1.001663 1.061333 1.021899 1.007367
 [57] 1.385392 1.016941 1.021180 1.019596 1.002943 1.000751 1.089853 1.003777 1.129796 1.060274 1.601423 1.491653 1.020008 1.004013
 [71] 1.000229 1.059399 1.171172 1.001009 1.397080 2.796381 1.734702 2.494200 1.000183 1.000432 1.006194 1.839790 1.005035 1.070056
 [85] 1.010485 1.000426 1.075276 1.000193 1.000260 1.006217 1.313430 1.022083 1.013131 2.330835 1.904411 1.225061 1.134951 1.000549
 [99] 1.000178 1.000737


###l_aw<- 1.0001 ; u_aw<- 1.5; l_l<- 0.001; u_l<- 2.5
> rbind(aw_abc,lam_abc)
             [,1]      [,2]     [,3]      [,4]      [,5]     [,6]     [,7]     [,8]     [,9]     [,10]    [,11]     [,12]
aw_abc  1.2141042 1.0756139 1.218152 1.1706533 1.1392095 1.232316 1.027751 1.211690 1.033784 1.2857404 1.123135 1.4073325
lam_abc 0.4175948 0.6731552 1.039760 0.5497753 0.7623193 1.203917 1.014485 1.099722 1.133694 0.2091019 0.801927 0.9096012
            [,13]     [,14]     [,15]     [,16]    [,17]    [,18]     [,19]    [,20]     [,21]     [,22]     [,23]
aw_abc  1.0748035 1.1619783 1.2342957 1.3388872 1.041489 1.319363 1.0295178 1.000296 1.1190786 1.0796789 1.0808808
lam_abc 0.5308388 0.7294215 0.6089678 0.4924454 1.421802 1.064906 0.9189674 1.883305 0.7951755 0.7643606 0.7139264
            [,24]    [,25]   [,26]    [,27]     [,28]    [,29]    [,30]    [,31]    [,32]    [,33]    [,34]    [,35]
aw_abc  1.2908917 1.020440 1.00725 1.093770 1.0168111 1.005829 1.219847 1.004571 1.001788 1.003695 1.111031 1.048393
lam_abc 0.4188022 1.709672 1.37139 1.074895 0.8434086 1.797435 1.122538 1.405083 1.349868 1.380532 1.178677 1.563556
           [,36]    [,37]    [,38]    [,39]     [,40]     [,41]    [,42]     [,43]    [,44]    [,45]    [,46]    [,47]
aw_abc  1.001656 1.000399 1.000325 1.038750 1.3049373 1.1479335 1.010463 1.0559318 1.078909 1.181692 1.016700 1.353212
lam_abc 1.285921 1.645982 1.320660 1.090282 0.1288348 0.1423612 1.622099 0.7413178 0.600081 0.177154 1.159837 1.301620
            [,48]    [,49]    [,50]    [,51]    [,52]    [,53]    [,54]    [,55]    [,56]    [,57]     [,58]     [,59]
aw_abc  1.0214433 1.004606 1.000600 1.000105 1.007632 1.163149 1.003051 1.018296 1.008347 1.058579 1.0059484 1.3591360
lam_abc 0.5711825 1.574420 2.092218 1.411694 1.676548 0.563640 1.504496 1.488266 1.616737 1.436298 0.8562962 0.4984508
            [,60]     [,61]     [,62]     [,63]    [,64]     [,65]    [,66]   [,67]     [,68]   [,69]    [,70]    [,71]
aw_abc  1.0323535 1.1245680 1.0242499 1.0949401 1.347935 1.0067132 1.489317 1.18366 1.0211387 1.00012 1.010222 1.144038
lam_abc 0.3872882 0.6482704 0.8148335 0.6312616 0.153033 0.9243034 1.319890 1.39144 0.6380137 1.63485 0.886001 1.665900
           [,72]    [,73]    [,74]     [,75]     [,76]    [,77]    [,78]    [,79]    [,80]    [,81]     [,82]    [,83]
aw_abc  1.034825 1.000928 1.000149 1.0435455 1.2268142 1.066204 1.001298 1.045412 1.395355 1.002380 1.2001518 1.453291
lam_abc 1.981784 1.593915 1.745016 0.4541295 0.5602029 1.479795 1.237457 1.256032 0.204230 1.170347 0.5735393 0.441945
            [,84]    [,85]     [,86]      [,87]   [,88]     [,89]    [,90]     [,91]    [,92]    [,93]    [,94]    [,95]
aw_abc  1.1398860 1.031643 1.1339364 1.09215465 1.00543 1.0224446 1.000210 1.0581438 1.000415 1.410859 1.000837 1.000685
lam_abc 0.4875511 1.334712 0.2173479 0.04188291 1.05983 0.9506416 1.455531 0.4449801 1.046055 1.140783 1.212290 1.225065
           [,96]    [,97]    [,98]    [,99]   [,100]
aw_abc  1.001464 1.005458 1.000388 1.050839 1.000175
lam_abc 1.034849 1.080461 1.847578 1.520367 1.248567

##############
THETA<- cbind(theta_aw, theta_l)
system.time({	 
apply(THETA,1,function(Z){  #x is a mtrix of dim= length(prior)x2=Ix2

	   LKinfo<- LKrigSetup( xr,LKGeometry="LKSphere",
                      NC=nc, a.wght=Z[1],
                      nlevel=nl, alpha=ALPHA, lambda=Z[2],
                      choleskyMemory=list(nnzR= 4E6)) 

	   ystar<- LKrig.sim(x, LKinfo=LKinfo, M=M) + sqrt(Z[2])*matrix(rnorm(n*M),n,M)
	   vg.t<- vgm(x,ystar, D.max,l)
	   rss<- sum((vg.t[,3]-vg[,3])^2)
	   rss
}
)
})
AAA<- matrix(c(1:12), ncol=2)

apply(AAA,1,function(Z)rnorm(n=2,mean=Z[1], sd=Z[2]))

###ensemble-median....I=150
####l_aw<- 1.0001 ; u_aw<- 1.5; l_l<- 0.001; u_l<- 2.5
> cbind(aw_abc,lam_abc)
        aw_abc   lam_abc
 [1,] 1.381304 0.9392599
 [2,] 1.048364 0.9927910
 [3,] 1.003872 1.2900582
 [4,] 1.028496 1.5556191
 [5,] 1.018221 1.1417306
 [6,] 1.246339 0.7742875
 [7,] 1.218930 0.5905058
 [8,] 1.022000 1.2079026
 [9,] 1.000920 1.2579534
[10,] 1.339887 0.7419879

####ensemble-1
####l_aw<- 1.0001 ; u_aw<- 1.5; l_l<- 0.001; u_l<- 2.5
       aw_abc   lam_abc
 [1,] 1.381304 0.9392599
 [2,] 1.048364 0.9927910
 [3,] 1.003872 1.2900582
 [4,] 1.028496 1.5556191
 [5,] 1.018221 1.1417306
 [6,] 1.022000 1.2079026
 [7,] 1.246339 0.7742875
 [8,] 1.218930 0.5905058
 [9,] 1.000920 1.2579534
[10,] 1.010396 1.6192053

####ensemble-15
####l_aw<- 1.0001 ; u_aw<- 1.5; l_l<- 0.001; u_l<- 2.5
        aw_abc   lam_abc
 [1,] 1.381304 0.9392599
 [2,] 1.048364 0.9927910
 [3,] 1.028496 1.5556191
 [4,] 1.003872 1.2900582
 [5,] 1.018221 1.1417306
 [6,] 1.022000 1.2079026
 [7,] 1.010396 1.6192053
 [8,] 1.246339 0.7742875
 [9,] 1.218930 0.5905058
[10,] 1.000920 1.2579534

####ensemble-80
####l_aw<- 1.0001 ; u_aw<- 1.5; l_l<- 0.001; u_l<- 2.5
        aw_abc   lam_abc
 [1,] 1.381304 0.9392599
 [2,] 1.028496 1.5556191
 [3,] 1.010396 1.6192053
 [4,] 1.003872 1.2900582
 [5,] 1.048364 0.9927910
 [6,] 1.022000 1.2079026
 [7,] 1.018221 1.1417306
 [8,] 1.246339 0.7742875
 [9,] 1.292762 1.4564964** is different than (1.218930 0.5905058) the above two
[10,] 1.000920 1.2579534


################################################################################
#l_aw<- 1.01 ; u_aw<- 2; l_l<- 0.01; u_l<- 2
cbind(aw_abc,lam_abc)
> cbind(aw_abc,lam_abc)
        aw_abc   lam_abc
 [1,] 1.081655 0.8367572
 [2,] 1.182263 0.4565120
 [3,] 1.684889 0.8995253
 [4,] 1.058007 1.1136761
 [5,] 1.642630 0.7487882
 [6,] 1.072755 1.2708178
 [7,] 1.016936 0.6547119
 [8,] 1.042396 0.9779666
 [9,] 1.079297 0.9178756
[10,] 1.814891 0.8739516


l_aw<- 1.0001 ; u_aw<- 1.5; l_l<- 0.01; u_l<- 2.5
> cbind(aw_abc,lam_abc)
        aw_abc   lam_abc
 [1,] 1.218930 0.9032129
 [2,] 1.246339 1.0937076
 [3,] 1.339887 1.0612840
 [4,] 1.048364 1.3036041
 [5,] 1.455612 0.6264803
 [6,] 1.000920 1.5408266
 [7,] 1.002545 1.3666948
 [8,] 1.159960 0.7582929
 [9,] 1.262670 0.3163778
[10,] 1.000255 1.6376954

l_aw<- 1.001 ; u_aw<- 1.5; l_l<- 0.01; u_l<- 2.5
> cbind(aw_abc,lam_abc)
        aw_abc   lam_abc
 [1,] 1.273627 0.9032129
 [2,] 1.298152 1.0937076
 [3,] 1.376853 1.0612840
 [4,] 1.050066 0.5389609
 [5,] 1.466422 0.6264803
 [6,] 1.016964 1.0142010
 [7,] 1.016307 1.1169976
 [8,] 1.007013 1.1934043
 [9,] 1.010700 1.3666948
[10,] 1.217761 0.7582929

##l_aw<- 1.001 ; u_aw<- 1.5; l_l<- 0.1; u_l<- 2.5
> cbind(aw_abc,lam_abc)
        aw_abc   lam_abc
 [1,] 1.273627 1.3840716
 [2,] 1.177225 0.6887232
 [3,] 1.312410 0.7514182
 [4,] 1.047904 1.0632804
 [5,] 1.466422 1.1185221
 [6,] 1.265473 0.9238390
 [7,] 1.100277 0.7670129
 [8,] 1.083368 0.8610302
 [9,] 1.122156 0.7137883
[10,] 1.217761 1.2500653

#l_aw<- 1.00001 ; u_aw<- 1.5; l_l<- 0.1; u_l<- 2.5**
cbind(aw_abc,lam_abc)
        aw_abc   lam_abc
 [1,] 1.175503 1.3840716
 [2,] 1.117815 1.2500653
 [3,] 1.445143 1.1185221
 [4,] 1.221180 0.7514182
 [5,] 1.166478 0.9238390
 [6,] 1.022055 0.8610302
 [7,] 1.001279 1.5663212
 [8,] 1.001370 1.4806995
 [9,] 1.042959 0.7137883
[10,] 1.005165 1.1429100

l_aw<- 1.00001 ; u_aw<- 1.5; l_l<- 0.01; u_l<- 2.5
> cbind(aw_abc,lam_abc)
        aw_abc   lam_abc
 [1,] 1.175503 0.9032129
 [2,] 1.306821 1.0612840
 [3,] 1.203865 1.0937076
 [4,] 1.025791 1.3036041
 [5,] 1.445143 0.6264803
 [6,] 1.000168 1.5408266
 [7,] 1.001044 1.5684967
 [8,] 1.221180 0.3163778
 [9,] 1.117815 0.7582929
[10,] 1.007466 1.4388637

l_aw<- 1.0001 ; u_aw<- 1.5; l_l<- 0.1; u_l<- 2.5
cbind(aw_abc,lam_abc)
        aw_abc   lam_abc
 [1,] 1.218930 1.3840716
 [2,] 1.021217 1.0246902
 [3,] 1.455612 1.1185221
 [4,] 1.262670 0.7514182
 [5,] 1.159960 1.2500653
 [6,] 1.210015 0.9238390
 [7,] 1.042756 0.8610302
 [8,] 1.120527 0.6887232
 [9,] 1.004795 1.4806995
[10,] 1.019968 1.0632804









