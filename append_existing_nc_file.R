rm(list=ls())
library(ncdf4)  

   ##~~for writing NC file
        xvals <- seq(-177.5, 177.5, 5)
        yvals <- seq(-87.5, 87.5, 5)
        nx <- length(xvals)
        ny <- length(yvals)
        lon1 <- ncdim_def("longitude", "degrees", xvals)
        lat2 <- ncdim_def("latitude", "degrees", yvals)
        time <- ncdim_def("time","months",c(1:24) , unlim=TRUE)       ##24 as it is for two years 2017 and 2018
        var_temp <- ncvar_def("temperature", "celsius", list(lon1, lat2, time), longname="Temperature_anomaly")
	 
		    ##write nc-files
		    for(i in 1:100){
                for(j in 1:100){   
		    filename=paste("D:/old_ens_17_18/hd_ens",i,"/ensemble",j,".nc",sep="")
       	   
                filename<- paste("C:/Users/maryam/Desktop/test.nc")
		    ncnew<- nc_create(filename, list(var_temp))
                data<- c(rep(NA,(24*2592)))
		    ncvar_put(ncnew, var_temp, data)
                nc_close(ncnew)     	#Don't forget to close the file
                }
		    }


rm(list=ls())
id=1
filename<- paste("C:/Users/maryam/Desktop/test.nc")
xvals <- seq(-177.5, 177.5, 5)
yvals <- seq(-87.5, 87.5, 5)
lon1 <- ncdim_def("longitude", "degrees", xvals)
lat2 <- ncdim_def("latitude", "degrees", yvals)
time <- ncdim_def("time","months",c(1:24) , unlim=TRUE)    
ncnew<- nc_open(filename, write=TRUE)
var_temp <- ncvar_def("temperature", "celsius", list(lon1, lat2, time), longname="Temperature_anomaly")
data<- c(1:2592)
ncvar_put(ncnew, var_temp, data,start=c(1,1,id),count=c(-1,-1,1))
nc_close(ncnew) 

rm(list=ls())
id=25
filename<- paste("C:/Users/maryam/Desktop/test.nc")
xvals <- seq(-177.5, 177.5, 5)
yvals <- seq(-87.5, 87.5, 5)
lon1 <- ncdim_def("longitude", "degrees", xvals)
lat2 <- ncdim_def("latitude", "degrees", yvals)
time <- ncdim_def("time","months",c(1:36) , unlim=TRUE)    
ncnew<- nc_open(filename, write=TRUE)
var_temp <- ncvar_def("temperature", "celsius", list(lon1, lat2, time), longname="Temperature_anomaly")
data<- c(1:2592)
ncvar_put(ncnew, var_temp, data,start=c(1,1,id),count=c(-1,-1,1))
nc_close(ncnew) 





 ncnew<- nc_open(filename, write=TRUE)
                        #ncnew<- nc_create(filename, list(var_temp)) #for first months
                        data<- c(LK_SIM[,j])
                        ncvar_put(ncnew, var_temp, data,start=c(1,1,id-id2),count=c(-1,-1,1))
                        nc_close(ncnew)     #Don't forget to close the file







  