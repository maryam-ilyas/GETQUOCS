rm(list=ls())
library(ncdf4)

##create file with NA's (NaN)
		xvals <- seq(1, 11, 1)  #1-11 parameters values and 101 is the sqrt(rss)=epsilon
		yvals <- seq(1, 2, 1) 
		nx <- length(xvals)
		ny <- length(yvals)
		lam <- ncdim_def("lambda", "usual", xvals)
		aw <- ncdim_def("a.wght", "usual", yvals)

		time <- ncdim_def("time","months",c(1:24) , unlim=TRUE)
		var_temp <- ncvar_def("parameters", "usual", list(lam, aw, time), longname="lk parameters") 

		filename=paste("C:/Users/maryam/Dropbox/h_ths_ens/parameters10from250.nc",sep="")
		
		#for first months 
		ncnew<- nc_create(filename, list(var_temp))
		data<- c(rep(NA,(24*22)))
		ncvar_put(ncnew, var_temp, data)
		nc_close(ncnew)	

		##Check
		ncnew<- nc_open(filename, write=TRUE) 
		AA<- ncvar_get(ncnew,"parameters")
      	dim(AA)
            nc_close(ncnew)
####################################################################################################
##Enter data
rm(list=ls())

		##create file with NA's (NaN)
		xvals <- seq(1, 11, 1)  #1-100 parameters values and 101 is the sqrt(rss)=epsilon
		yvals <- seq(1, 2, 1) 
		nx <- length(xvals)
		ny <- length(yvals)
		lam <- ncdim_def("lambda", "usual", xvals)
		aw <- ncdim_def("a.wght", "usual", yvals)

		time <- ncdim_def("time","months",c(1:2016) , unlim=TRUE)
		var_temp <- ncvar_def("parameters", "usual", list(lam, aw, time), longname="lk parameters") 

id<- 133:144
H=1
		filename=paste("D:/lk_par_ens/parameters10from250.nc",sep="")
		ncnew<- nc_open(filename, write=TRUE) 
		lam_abc<- 1:11
            aw_abc<-  5:15	 
		data<- c(cbind(lam_abc, aw_abc))
		ncvar_put(ncnew, var_temp, data,start=c(1,1,id[H]),count=c(-1,-1,1))
		nc_close(ncnew)

