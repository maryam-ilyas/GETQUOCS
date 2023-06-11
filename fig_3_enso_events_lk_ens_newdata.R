rm(list=ls())
library(ncdf4)
library(fields)

##ENSO events over lattice kriging ensembles taken from hd median
##R-files of input files computation:sst_nin_3_4_hd_and_lk_ens_(from_hdmed).R

			##FUNCTIONS
			##Function of monthly global average
			##sum(x_i_gridpoint * areaof_x_i_gridpoint)/sum(areaof_x_i_gridpoint)
			##and W= cos(lat*0.01745329)      
     
            	MON_AVG<- function(A, lon, lat){  #A-3D matrix (lon, lat, time) 		
			A1<- rep(0,dim(A)[3])
                  w<- cos(lat*(pi/180)) 	    #weights  
			W<- matrix(rep(w,each=length(lon)),nrow=length(lon)) #weights matrix

			for (i in 1:dim(A)[3]){
                        A2<- A[,,i]
                        W1_temp<- c(W)
                        W1_temp[is.na(c(A2))]<- NA
                        W1<- matrix(W1_temp, nrow=nrow(A2), ncol=ncol(A2)) 
                        A3<- A2*W1
                        A1[i]<- sum(A3,2, na.rm=TRUE)/sum(W1, na.rm=TRUE)                    
			}
			A1
			}

			op<- par(no.readonly=TRUE)
			par(oma=c(2,2,2,2))
			plot(1,1,type="n",xlab="",ylab="",xaxt="n",yaxt="n")
			for(side in 1:4){
  			inner<-round(par()$mar[side],0)-1
  			for(line in 0:inner){
  			mtext(text=paste0("Inner line ",line),side=side,line=line)
  			}
  			outer<-round(par()$oma[side],0)-1
   			for(line in 0:inner){
   			mtext(text=paste0("Outer line ",line),side=side,line=line,outer=TRUE)
  			}
			}

##HadCRUT4 median ensemble member
      a1<- "D:/HadCRUT.4.5.0.0.median.nc"  #downloaded 16th Feb. 2017
     	fid.t<- nc_open(a1,verbose=TRUE,write=FALSE)
	time.t<- ncvar_get(fid.t,'time')     
	origin.t <- fid.t$dim$time$units
	tunits.t <- strsplit(fid.t$dim$time$units,split=" ") 
	if(identical(tunits.t[[1]][1],"minutes")) tunits.t[[1]][1]<-"mins"
	TIME.t<- strptime(paste(tunits.t[[1]][3],tunits.t[[1]][4]),"%Y-%m-%d %H:%M:%S",tz="GMT")+
                           as.numeric(as.difftime(time.t,units=tunits.t[[1]][1]),units="secs")
				   #1850-01-16 to 2015-05-16
	lat.t <- ncvar_get(fid.t,'latitude');lat.t                 		#length(lat.t)         #36
	lon.t <- ncvar_get(fid.t,'longitude'); lon.t 
      time.t<- format(TIME.t, "%b %Y")
      n<- length(TIME.t)
      p<- 2592

BP<- c(1961,1990)   	   #30 year base-period span
clm.n.ind<- which(TIME.t>= paste(BP[1],'-01-16',sep="") & TIME.t<= paste(BP[2],'-12-17',sep="")) #climatological period

##Raw SST anomalies in the Nino 3.4 region over lattice kriging ensemble
sst_n34_lk<- read.table("D:/lk_par_ens_new_data/nin_3_4_par_ens.txt")

n.ens<- ncol(sst_n34_lk)

##extract winter PDF
yrs<- (1851:2016)

avg_DJF_lk<- matrix(0, nrow=length(yrs), ncol=ncol(sst_n34_lk))

for(i in 1:ncol(sst_n34_lk)){
K<- matrix(c(sst_n34_lk[,i]),ncol=12,byrow=TRUE)    #take an ensemble member
	
	dec<- K[-nrow(K),12]                          #remove last
	jan_feb<- K[-1,c(1,2)]			          #remove 1st
	DJF<- cbind(dec,jan_feb)
	avg_DJF_lk[,i]<- apply(DJF,1,mean)
	}

##Spread of each PDF
SD_LK<- as.vector(apply(avg_DJF_lk,1,sd))

W<- rep(NA, length(yrs)); C<- rep(NA, length(yrs)); G<- rep(NA,length(yrs))
               
		for(i in 1:length(yrs)){
		sst_tmp<- as.numeric(avg_DJF_lk[i,])
		W[i]<- sum(sst_tmp>= 0.5)/n.ens
		C[i]<- sum(sst_tmp<= -0.5)/n.ens  	     
		G[i]<- sum(sst_tmp>= -0.25 & sst_tmp<= 0.25)/n.ens			                            
		}

SSS<- rep(0,length(yrs))    
SSS[which(SD_LK>0.5)]<- "Ambiguous" ; length(SSS[which(SD_LK>0.5)])
SSS[which(G>=0.66)]<- "Neutral"  ;length(SSS[which(G>=0.66)])
SSS[which(C>=0.66)]<- "La Niña"  ;length(SSS[which(C>=0.66)])
SSS[which(W>=0.66)]<- "El Niño"  ;length(SSS[which(W>=0.66)])
SSS[which(SSS==0)]<- "Transit"
DF_lk<- data.frame(yrs,SSS)
DF_lk

#write.table(DF_lk,"D:/lk_par_ens_new_data/lk_enso_ev.txt")

new<- read.table("D:/lk_par_ens_new_data/lk_enso_ev.txt")
########################################################################################################################################
##check
AA<- read.table("C:/Users/maryam/Dropbox/PhD UCL/kriging_and_FEOFs/post_upgrade/paper_draft/submission/supp_files/lk_enso_ev.txt")

all(DF_lk[,2]==AA[,2])

el<- which(DF_lk[,2]=="El Niño")
ln<-  which(DF_lk[,2]=="La Niña")
nu<- which(DF_lk[,2]=="Neutral")
tr<- which(DF_lk[,2]=="Transit")
am<- which(DF_lk[,2]=="Ambiguous")

#Check respective probabilities

W[el] 
C[ln]
G[nu]


X11()
par(op)
#e<-  paste("C:/Users/maryam/Dropbox/PhD UCL/kriging_and_FEOFs/post_upgrade/paper_r_files/figs/fig3_r.png",sep="")
#png(e,width = 6 , height = 10, units = 'in', res = 300)
par(oma=c(1.5,1.5,1.5,1),mar=c(3.1,2.1,1,1),mfrow=c(3,2))

hist(avg_DJF_lk[el[3],], prob=TRUE,ylim=c(0,2),xlim=c(-3,3),breaks=c(seq(-3,3,0.125)),
main=substitute(paste("(a) El Niño (", EL ,",",~sigma,"=",ELS,~degree~C,")"), 
list(EL = yrs[el[3]],ELS=round(SD_LK[el[3]],2)) ) )
abline(v=0.5, col="red", lwd=1.5)
abline(v=-0.5, col="blue", lwd=1.5)
abline(v=0, col="black", lwd=1.5)


hist(avg_DJF_lk[ln[8],], prob=TRUE,ylim=c(0,2),xlim=c(-3,3),breaks=c(seq(-3,3,0.125)),
					#main=paste("(b) La Niña", "(yrs[ln[5]],", SD =",round(SD_LK[ln[5]],2),"°C)") )
main=substitute(paste("(b)  La Niña (", LN ,",",~sigma,"=",LNS~degree~C,")"), 
list(LN = yrs[ln[8]],LNS=round(SD_LK[ln[8]],2)) ) )
abline(v=0.5, col="red", lwd=1.5)
abline(v=-0.5, col="blue", lwd=1.5)
abline(v=0, col="black", lwd=1.5)

hist(avg_DJF_lk[nu[1],], prob=TRUE,ylim=c(0,2),xlim=c(-3,3),breaks=c(seq(-3,3,0.125)),
					#main=paste("(c) Neutral definite", "(yrs[nu[1]],", SD =",round(SD_LK[nu[1]],2),"°C)") )
main=substitute(paste("(c)  Definite-neutral (", ND ,",",~sigma,"=",NDS~degree~C,")"), 
list(ND = yrs[nu[1]],NDS=round(SD_LK[nu[1]],2)) ) )

abline(v=0.5, col="red", lwd=1.5)
abline(v=-0.5, col="blue", lwd=1.5)
abline(v=0, col="black", lwd=1.5)


hist(avg_DJF_lk[am[43],], prob=TRUE,ylim=c(0,2),xlim=c(-3,3),breaks=c(seq(-3,3,0.125)),
					#main=paste("(e) Ambiguous","(yrs[am[43]],", SD =",round(SD_LK[am[43]],2),"°C)") )
main=substitute(paste("(d)  Ambiguous (", AM ,",",~sigma,"=",AMS~degree~C,")"), 
list(AM = yrs[am[43]],AMS=round(SD_LK[am[43]],2)) ) )
abline(v=0.5, col="red", lwd=1.5)
abline(v=-0.5, col="blue", lwd=1.5)
abline(v=0, col="black", lwd=1.5)


hist(avg_DJF_lk[tr[26],], prob=TRUE,ylim=c(0,2),xlim=c(-3,3),breaks=c(seq(-3,3,0.125)),
					#main=paste("(d) Neutral residual", "(yrs[tr[26]],", SD =",round(SD_LK[tr[26]],2),"°C)") )
main=substitute(paste("(e)  Residual-neutral (", RN ,",",~sigma,"=",RNS~degree~C,")"), 
list(RN = yrs[tr[26]],RNS=round(SD_LK[tr[26]],2)) ) )
abline(v=0.5, col="red", lwd=1.5)
abline(v=-0.5, col="blue",lwd=1.5)
abline(v=0, col="black", lwd=1.5)

mtext(text="Niño 3.4 SST anomaly",side=1,line=0,outer=TRUE)
#mtext(text="Frequency",side=2,line=0,outer=TRUE)
dev.off()

