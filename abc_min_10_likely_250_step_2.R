rm(list=ls())
library(ncdf4)
library(fields)
library(LatticeKrig)
library(lhs)

	##~~HadCRUT4 median 
	source("C:/Users/maryam/Dropbox/Project_NCAR/abc_lk/abc_lk_functions.R")
	##HadCRUT4 median ensemble member
      ###--del##a1<- "D:/hd_sprs_dat/HadCRUT.4.5.0.0.median.nc"  
	a1<- paste("D:/4_4_19/HadCRUT.4.6.0.0.anomalies.1.nc", sep="")
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

	##locations
	X<- as.matrix((expand.grid(lon.t,lat.t))) 
##***
	id<- c(1:24)  
  	D.max<- 6300 	
    	l<- 13
	xr=cbind( c(-180, 180), c(-90,90))
	nc=3
	nl=3
	ALPHA<- c(0.2450707, 0.01606482, 0.7388645)
	N<- nrow(X)
	M<- 1	
      
	##simulate from the prior
	I<- 250    		#No. of iterations
	LKinfo1<- LKrigSetup( xr, LKGeometry="LKSphere",
                      NC=nc, a.wght=1.1,lambda=0.6,
                      nlevel=nl, alpha=ALPHA, choleskyMemory=list(nnzR= 4E6))
	l.omg<- -4.5	; omega2Awght(l.omg, LKinfo1) #1.000123
	u.omg<- 0.55	; omega2Awght(u.omg, LKinfo1) #4.004166

	l.lnl<- -6.9	; exp(l.lnl)  
	u.lnl<- 1.4       ; exp(u.lnl)  
	set.seed(1224)
	L<- randomLHS(I,2)      #draw a Latin hypercube first from a set of standard uniform distribution
	LL<- cbind(qunif(L[,1],l.lnl,u.lnl), qunif(L[,2], l.omg,u.omg))  #transform the sample from U(0,1) to U(0.01,6) and U(0.01,16)

	theta_l<- exp(LL[,1])     			#summary(theta_l)
	theta_aw<- omega2Awght(LL[,2], LKinfo1)   #summary(theta_aw)

		 ##For writing .nc files
                xvals <- seq(1, 11, 1)  #1-10 parameters values and 11 is the sqrt(rss)=epsilon
                yvals <- seq(1, 2, 1)
                nx <- length(xvals)
                ny <- length(yvals)
                lam <- ncdim_def("lambda", "usual", xvals)
                aw <- ncdim_def("a.wght", "usual", yvals)
                time <- ncdim_def("time","months",c(1:24) , unlim=TRUE)
                var_temp <- ncvar_def("parameters", "usual", list(lam, aw, time), longname="lk parameters")

for(k in 1:length(id)){
                
		Y<- c(tempr.a[,,id[k]])
		ind_s<- !is.na(Y)
		y<- Y[ind_s]  
		x<- X[ind_s,]

    		##~~~empirical variogram 
    		DD<- rdist.earth(x,miles=FALSE)[upper.tri(rdist.earth(x,miles=FALSE))];	round(max(DD)/2)
    		vg<- vgm(x,y,D.max,l)
    		n<- nrow(x)

		#one sample for each parametric value
		rss<- rep(0,I)
			for(i in 1:I){   			   		#3 sec for one iteration
				LKinfo<- LKrigSetup( xr,LKGeometry="LKSphere",
                      	NC=nc, a.wght=theta_aw[i],
                      	nlevel=nl, alpha=ALPHA, lambda=(theta_l[i]),
                      	choleskyMemory=list(nnzR= 4E6))                       
				ystar<- LKrig.sim(x, LKinfo=LKinfo, M=M) + sqrt(theta_l[i])*matrix(rnorm(n*M),n,M)	
				vg.t<- vgm(x,ystar, D.max,l)
				rss[i]<- sum((vg.t[,3]-vg[,3])^2)  #residual sums of square	
				}
 		 ind<- (order(sqrt(rss)))[1:10]
  		 ep<- (sqrt(rss[ind]))[10]
  		 aw_abc<- c(theta_aw[ind],ep)
  		 lam_abc<- c(theta_l[ind],ep)

 		    ##write in .nc file
		    f_name=paste("C:/Users/maryam/Dropbox/h_ths_ens/parameters10from250.nc",sep="")
                ncnew<- nc_open(f_name, write=TRUE)
                data_nc<- c(cbind(lam_abc, aw_abc))
                ncvar_put(ncnew, var_temp, data_nc,start=c(1,1,id[k]),count=c(-1,-1,1))
                nc_close(ncnew)

		rm(list=ls())
		library(ncdf4)
		library(fields)
		library(LatticeKrig)
		library(lhs)
 
		source("C:/Users/maryam/Dropbox/Project_NCAR/abc_lk/abc_lk_functions.R")
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

		##locations
		X<- as.matrix((expand.grid(lon.t,lat.t))) 
##*****
		id<- c(1:24)             ##TIME.t[id]	
		D.max<- 6300 	
    		l<- 13
		xr=cbind( c(-180, 180), c(-90,90))
		nc=3
		nl=3
		ALPHA<- c(0.2450707, 0.01606482, 0.7388645)
		N<- nrow(X)
		M<- 1	
      
		##simulate from the prior
		I<- 250    #No. of iterations
		LKinfo1<- LKrigSetup( xr, LKGeometry="LKSphere",
                      NC=nc, a.wght=1.1,lambda=0.6,
                      nlevel=nl, alpha=ALPHA, choleskyMemory=list(nnzR= 4E6))
		l.omg<- -4.5	; omega2Awght(l.omg, LKinfo1)
		u.omg<- 0.55	; omega2Awght(u.omg, LKinfo1)

		l.lnl<- -6.9	; exp(l.lnl)  
		u.lnl<- 1.4       ; exp(u.lnl)  
		set.seed(1224)
		L<- randomLHS(I,2)       #draw a Latin hypercube first from a set of standard uniform distribution
		LL<- cbind(qunif(L[,1],l.lnl,u.lnl), qunif(L[,2], l.omg,u.omg))  #transform the sample from U(0,1) to U(0.01,6) and U(0.01,16)

		theta_l<- exp(LL[,1])     			#summary(theta_l)
		theta_aw<- omega2Awght(LL[,2], LKinfo1)   #summary(theta_aw)

		 ##For writing .nc files
                xvals <- seq(1, 11, 1)  #1-10 parameters values and 11 is the sqrt(rss)=epsilon
                yvals <- seq(1, 2, 1)
                nx <- length(xvals)
                ny <- length(yvals)
                lam <- ncdim_def("lambda", "usual", xvals)
                aw <- ncdim_def("a.wght", "usual", yvals)
                time <- ncdim_def("time","months",c(1:24) , unlim=TRUE)
                var_temp <- ncvar_def("parameters", "usual", list(lam, aw, time), longname="lk parameters")
}

rm(list=ls())

#f_name=paste("D:/parameters10from250.nc")
#f_name=paste("D:/lk_par_ens_new_data/lam_aw_par.nc")

f_name=paste("C:/Users/maryam/Dropbox/h_ths_ens/parameters10from250.nc",sep="")
ncnew<- nc_open(f_name, write=TRUE)
FF<- ncvar_get(ncnew,"parameters")
nc_close(ncnew)






