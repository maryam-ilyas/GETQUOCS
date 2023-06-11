rm(list=ls())
library(ncdf4)
library(fields)
library(LatticeKrig)
library(chron)
library(nlme)

##~~HadCRUT4 median descriptives
##~~1 GISS, HadCRUT4 median (with & without gap filling) annual averages and
##~~ linear time trends

##FUNCTIONS
		##Function of monthly global average
		##sum(x_i_gridpoint * areaof_x_i_gridpoint)/sum(areaof_x_i_gridpoint)
		##and W= cos(lat*0.01745329)      
     
            MON_AVG<- function(A, lon, lat){  #A-3D matrix (lon, lat, time) 
		
			A1<- rep(0,dim(A)[3])
                  w<- cos(lat*(pi/180)) #weights  
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

		##function to calculate annual averages of a time series
		ann_avg<- function(x,S1){    
                aa<- ts(x,start=c(S1,1),frequency=12)
		    bb<- as.vector(colMeans(matrix(window(aa,start=c(S1,1)),nrow=12)))
      	bb
		}

	##Function to calculate 95% prediction intervals for gls using the bootstrap approach
      ##See the reference below but do not sample errors as it results into weird confidence bands
      ##https://www.youtube.com/watch?v=c3gD_PwsCGM
      ##x:independent variable, y:dependent variable, B:no.of bootstrap samples, corr=autpregressive errors

	gls_pred_int<- function(x,y,B,corr){ 
	Y_hat<- matrix(0, nrow=length(x), ncol=B)
	for(i in 1:B){
		ind<- sort(sample(1:length(x),size=length(x) ,replace=TRUE))
            x_bs<- x[ind] ; y_bs<- y[ind]
		reg_bs<- gls(y_bs~x_bs, correlation=corr)   #bootstarp model

		y_hat<- predict(reg_bs,newdata=data.frame(x_bs=c(x)))#obtain predictions using bootstrap model
		Y_hat[,i]<- y_hat 
		}
	ci<- apply (Y_hat,1,quantile, probs=c(0.025,0.975)) #95% prediction interval for every point
	t(ci)
	}

##HadCRUT4 median ensemble member
	a1<- "D:/hd_sprs_dat/HadCRUT.4.5.0.0.median.nc" #downloaded 16th Feb. 2017
	fid.t<- nc_open(a1,verbose=TRUE,write=FALSE)
	time.t<- ncvar_get(fid.t,'time')     
	origin.t <- fid.t$dim$time$units
	tunits.t <- strsplit(fid.t$dim$time$units,split=" ") 
	if(identical(tunits.t[[1]][1],"minutes")) tunits.t[[1]][1]<-"mins"
	TIME.t<- strptime(paste(tunits.t[[1]][3],tunits.t[[1]][4]),"%Y-%m-%d %H:%M:%S",tz="GMT")+
                           as.numeric(as.difftime(time.t,units=tunits.t[[1]][1]),units="secs")
	rng_TIME<- c(format(TIME.t[1], "%B-%Y"),format(TIME.t[length(TIME.t)], "%b-%Y")) ;rng_TIME   
	lat.t <- ncvar_get(fid.t,'latitude');lat.t                 		#length(lat.t)         #36
	lon.t <- ncvar_get(fid.t,'longitude'); lon.t               		#length(lon.t)         #72
	tempr.a <- ncvar_get(fid.t,'temperature_anomaly')          		#dim(tempr)            #72 36 1985

##cut off at Dec 2016
tempr.a<- tempr.a[,,c(1:2004)]
TIME.t<- TIME.t[1:2004]

hd<- tempr.a
##gap-filled data
E1<- paste("D:/hadcrut_v4_5_0_0_gapfilled_olddata/lk_hd_median.txt")
E2<- paste("D:/hadcrut_v4_5_0_0_gapfilled_olddata/lk_hd_median_p3.txt")

lk_med_1<- read.table(E1,sep="\t")
lk_med_2<- read.table(E2,sep="\t")
lk_med<- rbind(lk_med_1, lk_med_2)
lk<- array(t(lk_med),dim=c(length(lon.t), length(lat.t),length(TIME.t)))

###Monthly average temperature###
med.lk<- MON_AVG(lk,lon.t,lat.t)
med.had<- MON_AVG(hd,lon.t,lat.t)

##GISS--Gridded data (http://www.esrl.noaa.gov/psd/data/gridded/data.gistemp.html)
##Combined SST/Air Temperature.....Monthly Anomaly: 250km smoothed
b1<- "D:/giss_anomalies/air.2x2.1200.mon.anom.comb.nc"
fid<- nc_open(b1,verbose=TRUE, write=FALSE)
lon <- ncvar_get(fid,'lon')             #length(lon)....180
lat <- ncvar_get(fid,'lat')             #length(lat)....90
time<- ncvar_get(fid,'time')
air<- ncvar_get(fid,'air')              #180   90 1653, 44

##cutoff time and temperature upto 2016
time<- time[1:1644]
air<- air[,,c(1:1644)]

	tunits <- ncatt_get(fid, "time", "units")
	tustr <- strsplit(tunits$value, " ")
	tdstr <- strsplit(unlist(tustr)[3], "-")
	tmonth = as.integer(unlist(tdstr)[2])
	tday = as.integer(unlist(tdstr)[3])
	tyear = as.integer(unlist(tdstr)[1])
Tm<- chron(time, origin = c(tmonth, tday, tyear))    #T[1];T[length(T)]; 01/01/80 -- 8/01/16
T<- as.POSIXlt(Tm, format="%m/%d/%Y")
rng_T<- c(format(T[1], "%B-%Y"),format(T[length(T)], "%B-%Y"));rng_T
AIR<- air
giss_mth_avg<- (MON_AVG(AIR,lon,lat))                #Jan.1880 to Aug.2016

##annual averages
hd_ind<- (12*trunc(length(TIME.t)/12))    #leave months at the end if not a full year 
giss_ind<- (12*trunc(length(T)/12))       #leave months at the end if not a full year 

	S1<- as.numeric(format(TIME.t[1],'%Y')) #beginning HD
	S2<- as.numeric(format(T[1],'%Y'))      #beginning GISS
	S3<- 1979                               #Cut
	S_end<- as.numeric(format(TIME.t[hd_ind],"%Y"))#end HD

ann_lk<- ann_avg(med.lk[1:hd_ind],S1)
ann_hd<- ann_avg(med.had[1:hd_ind],S1)
ann_giss<- ann_avg(giss_mth_avg[1:giss_ind],S2)

SSS<- seq(as.Date(TIME.t[1]),as.Date(TIME.t[hd_ind]),by = 'year')   #1850-2015...HD
YR<- SSS[-which(SSS< paste(S2,"-01-16", sep=""))];YR   #1880-2015
ANN_LK<- ann_lk[-which(SSS< paste(S2,"-01-16", sep=""))]    
ANN_HD<- ann_hd[-which(SSS< paste(S2,"-01-16", sep=""))]
ANN_GISS<- ann_giss   		

SS<- seq(as.Date(paste(S2,"-12-01", sep="")),as.Date(paste(S_end,"-12-01", sep="")),by = 'year')
SS.ind2<- which(SS>=paste(S3,"-01-1", sep=""))

##Cowtan & Way (2014)
fff<- "C:/Users/maryam/Dropbox/thesis/latex_thesis/r_files/cow_way_anomalies/had4_krig_v2_0_0.nc"
fid_cw<- nc_open(fff,verbose=TRUE, write=FALSE)
temp_cw<- ncvar_get(fid_cw,'temperature_anomaly')
cw_mth_avg<- (MON_AVG(temp_cw,lon.t,lat.t))  
ann_cw<- ann_avg(cw_mth_avg[1:hd_ind],S1)
ANN_CW<- ann_cw[-which(SSS< paste(S2,"-01-16", sep=""))]

	##Check errors in all the series over longer and shorter time period
      X11()
	set.panel(3,4)
      acf(resid(lm(ANN_GISS~time(ANN_GISS))))
	pacf(resid(lm(ANN_GISS~time(ANN_GISS)))) #geometric decay of acf and pacf significant at p=1 so AR(1)

      acf(resid(lm(ANN_GISS~time(ANN_GISS))))
	pacf(resid(lm(ANN_GISS[SS.ind2]~time(ANN_GISS[SS.ind2]))))

      acf(resid(lm(ANN_CW~time(ANN_CW))))
	pacf(resid(lm(ANN_CW~time(ANN_CW)))) #geometric decay of acf and pacf significant at p=1 so AR(1)

      acf(resid(lm(ANN_CW~time(ANN_CW))))
	pacf(resid(lm(ANN_CW[SS.ind2]~time(ANN_CW[SS.ind2]))))

      acf(resid(lm(ANN_LK~time(ANN_LK))))
	pacf(resid(lm(ANN_LK~time(ANN_LK)))) #geometric decay of acf and pacf significant at p=1 so AR(1)

      acf(resid(lm(ANN_LK~time(ANN_LK))))
	pacf(resid(lm(ANN_HD[SS.ind2]~time(ANN_HD[SS.ind2]))))

	#library(forecast)
	#auto.arima(as.ts(resid(lm(ANN_GISS~time(ANN_GISS)))))
	#auto.arima(as.ts(resid(lm(ANN_GISS[SS.ind2]~time(ANN_GISS[SS.ind2]))))) #no autocorrelation confirm uisng the test below
      #library(lmtest)
      #dwtest(lm(ANN_GISS[SS.ind2]~time(ANN_GISS[SS.ind2])))
      #p-value = 0.1115, alternative hypothesis: true autocorrelation is greater than 0
     
	#auto.arima(as.ts(resid(lm(ANN_HD~time(ANN_HD)))))
	#auto.arima(as.ts(resid(lm(ANN_HD[SS.ind2]~time(ANN_HD[SS.ind2]))))) #no autocorrelation confirm uisng the test below
      #dwtest(lm(ANN_HD[SS.ind2]~time(ANN_HD[SS.ind2])))
	#DW = 1.5886, p-value = 0.07262


	#auto.arima(as.ts(resid(lm(ANN_LK~time(ANN_LK)))))
	#auto.arima(as.ts(resid(lm(ANN_LK[SS.ind2]~time(ANN_LK[SS.ind2]))))) #no autocorrelation confirm uisng the test below
      #dwtest(lm(ANN_LK[SS.ind2]~time(ANN_LK[SS.ind2])))
	#DW = 1.7531, p-value = 0.1738


#X_LAB_3<- paste("Time (",S2,"-",S_end,")")

##Important check
##		X11()
##		plot(time(ANN_HD),ANN_HD,type='l',
##		ylab="Global average temperature anomaly (deg C)",xlab= X_LAB_3,
##		col="red",ylim=c(-0.6,0.7))
##		axis(side = 3,at = seq(1,length(ANN_HD), by=10), 
##         	labels = seq(S2,S_end, by=10),tck= -.02)
##		cbind(seq(1,length(ANN_HD),by=10),seq(S2,S_end, by=10))

X11()
e<- paste("C:/Users/maryam/Dropbox/thesis/latex_thesis/figures/figure_3_new_r.png",sep="")
png(e,width = 7, height = 6, units = 'in', res = 300)
layout(matrix(c(1,1,2),nc=3, byrow = TRUE))
par(mar=c(4.1,4.1,0.5,0.5), oma=c(0,0,0,0)) 
##par(oma=c(0,1,0,0))
ld<- 1.2 
plot(time(ANN_CW),ANN_CW,xaxt="n",type='n',yaxt="n",
ylab="Annual, global mean temperature anomaly (°C)",xlab= "Year",
col="green",ylim=c(-0.6,1))
axis(side = 1,at = seq(1,length(ANN_CW), by=10), 
          labels = seq(S2,S_end, by=10),tck= -.02)
axis(side=2, at= round(seq(-0.6,0.9, by=0.1),1))
abline(h=0, lty=2, col='gray')

##GISS     
      fit1.giss<- gls(ANN_GISS~time(ANN_GISS), correlation=corARMA(p=1))
	fit2.giss<- lm(ANN_GISS[SS.ind2]~time(ANN_GISS[SS.ind2]))
L2<- fit1.giss$fitted             
L22<- fit2.giss$fitted		     
	giss1.b<- summary(fit1.giss)$tTable[2];giss1.b
	giss1.p<- summary(fit1.giss)$tTable[8];giss1.p
giss1.ci<- as.vector(as.matrix(confint(fit1.giss, level=0.90))[2,]) #interval estimate of slope

	giss2.b<- summary(fit2.giss)$coefficients[2];giss2.b
	giss2.p<- summary(fit2.giss)$coefficients[8];giss2.p
      giss2.ci<- as.vector(as.matrix(confint(fit2.giss, level=0.90))[2,])
preds.giss1 <- gls_pred_int(x=time(ANN_GISS),y=ANN_GISS,B=1000, corr=corARMA(p=1))#prediction intervals
preds.giss2 <- predict(fit2.giss,interval = 'confidence',level=0.90)#interval estimate of predicted values

polygon(c(rev(time(ANN_GISS)),(time(ANN_GISS))),
c(rev(preds.giss1[ ,2]), preds.giss1[ ,1]),col = "grey75", border = FALSE)

polygon(c(rev(time(ANN_GISS)[SS.ind2]),(time(ANN_GISS)[SS.ind2])),
c(rev(preds.giss2[ ,3]), preds.giss2[ ,2]),col = "grey75", border = FALSE)

##HadCRUT4
	fit1.cw<- gls(ANN_CW~time(ANN_CW), correlation=corARMA(p=1))
	fit2.cw<- lm(ANN_CW[SS.ind2]~time(ANN_CW[SS.ind2]))
L4<- fit1.cw$fitted			
L44<- fit2.cw$fitted

	cw1.b<- summary(fit1.cw)$tTable[2];cw1.b
	cw1.p<- summary(fit1.cw)$tTable[8];cw1.p
      cw1.ci<- as.vector(as.matrix(confint(fit1.cw, level=0.90))[2,])

	cw2.b<- summary(fit2.cw)$coefficients[2];cw2.b
	cw2.p<- summary(fit2.cw)$coefficients[8];cw2.p
      cw2.ci<- as.vector(as.matrix(confint(fit2.cw, level=0.90))[2,])

	preds.cw1 <- gls_pred_int(x=time(ANN_CW),y=ANN_CW,B=1000, corr=corARMA(p=1))#prediction intervals
	preds.cw2 <- predict(fit2.cw,interval = 'confidence',level=0.90)

polygon(c(rev(time(ANN_CW)),(time(ANN_CW))),
c(rev(preds.cw1[ ,2]), preds.cw1[ ,1]),col = "yellow", border = FALSE)    #grey60

polygon(c(rev(time(ANN_CW)[SS.ind2]),(time(ANN_CW)[SS.ind2])),
c(rev(preds.cw2[ ,3]), preds.cw2[ ,2]),col = "yellow", border = FALSE)

##MRLK
	fit1.lk<- gls(ANN_LK~time(ANN_LK), correlation=corARMA(p=1))
	fit2.lk<- lm(ANN_LK[SS.ind2]~time(ANN_LK[SS.ind2]))
			
L1<- lm(ANN_LK~time(ANN_LK))$fitted
L11<- lm(ANN_LK[SS.ind2]~time(ANN_LK[SS.ind2]))$fitted  
 
	lk1.b<- summary(fit1.lk)$tTable[2];lk1.b
	lk1.p<- summary(fit1.lk)$tTable[8];lk1.p
      lk1.ci<- as.vector(as.matrix(confint(fit1.lk, level=0.90))[2,])

	lk2.b<- summary(fit2.lk)$coefficients[2];lk2.b
	lk2.p<- summary(fit2.lk)$coefficients[8];lk2.p
      lk2.ci<- as.vector(as.matrix(confint(fit2.lk, level=0.90))[2,])
	preds.lk1 <- gls_pred_int(x=time(ANN_LK),y=ANN_LK,B=1000, corr=corARMA(p=1))
	preds.lk2 <- predict(fit2.lk,interval = 'confidence',level=0.90)
			
polygon(c(rev(time(ANN_LK)),(time(ANN_LK))),
c(rev(preds.lk1[ ,2]), preds.lk1[ ,1]),col = "pink", border = FALSE) #grey90

polygon(c(rev(time(ANN_LK)[SS.ind2]),(time(ANN_LK)[SS.ind2])),
c(rev(preds.lk2[ ,3]), preds.lk2[ ,2]),col = "pink", border = FALSE)

lines(time(ANN_GISS),ANN_GISS,col="black",lwd=ld)
lines(time(ANN_GISS),L2,col="black",lwd=ld)
lines(time(ANN_GISS)[SS.ind2],L22,col="black",lwd=ld)
			
lines(time(ANN_HD),L4, col="green",lwd=ld)
lines(time(ANN_HD)[SS.ind2],L44, col="green",lwd=ld)
lines(time(ANN_HD),ANN_HD,col="green", lwd=ld)

lines(time(ANN_LK),ANN_LK,col=22,lwd=ld)
lines(time(ANN_LK),L1,col=22,lwd=ld)
lines(time(ANN_LK)[SS.ind2],L11,col=22,lwd=ld)

legend("topleft",bty='n', c("Cowtan & Way(2014)","Gap-filled HadCRUT4","GISS"),
col=c("green",22, "black"),lty=c(1,1,1), lwd=c(ld,ld,ld))

X_LAB_4<- paste(S2,"-",S_end,"&",S3,"-",S_end)
lim<- c(0.05,0.21)
plot(1,cw1.b*10, col='green', pch=19,xlab=X_LAB_4, ylab='Trend (°C/ decade)',ylim=lim,,xaxt='n', yaxt='n')
axis(2,at=seq(0.04, 0.2, by=0.01))
points(1,cw2.b*10, col='green', pch=3)
segments(x0=1,x1=1,y0=cw1.ci[1]*10,y1=cw1.ci[2]*10, col='green')
segments(x0=1,x1=1,y0=cw2.ci[1]*10,y1=cw2.ci[2]*10, col='green', lty=2)

points(1.2,giss1.b*10,col='black', pch=19)
points(1.2,giss2.b*10,col='black', pch=3)
segments(x0=1.2,x1=1.2,y0=giss1.ci[1]*10,y1=giss1.ci[2]*10, col='black')
segments(x0=1.2,x1=1.2,y0=giss2.ci[1]*10,y1=giss2.ci[2]*10, col='black', lty=2)

points(0.8,lk1.b*10,col=22, pch=19)
points(0.8,lk2.b*10,col=22, pch=3)
segments(x0=0.8,x1=0.8,y0=lk1.ci[1]*10,y1=lk1.ci[2]*10, col=22)
segments(x0=0.8,x1=0.8,y0=lk2.ci[1]*10,y1=lk2.ci[2]*10, col=22, lty=2)
dev.off()

##slope per decade
bhat_full<- c(giss1.b,cw1.b,lk1.b)*10 ;bhat_full
bhat_cut<-  c(giss2.b,cw2.b,lk2.b)*10 ;bhat_cut

p_full<- c(giss1.p,cw1.p,lk1.p) ;p_full
p_cut<-  c(giss2.p,cw2.p,lk2.p) ;p_cut
DF<- data.frame(bhat_full,p_full,bhat_cut,p_cut)
rownames(DF)<- c("GISS", "Cowtan & Way","HadCRUT4 (gap-filled)")
DF


##> DF
                       bhat_full       p_full  bhat_cut        p_cut
GISS                  0.07579064 1.199503e-10 0.1698860 4.886738e-14
Cowtan & Way          0.07170704 3.012728e-14 0.1868980 1.387001e-15
HadCRUT4 (gap-filled) 0.06999475 1.193274e-14 0.1779717 6.317052e-15









> DF
                       bhat_full       p_full  bhat_cut        p_cut
GISS                  0.07579064 1.199503e-10 0.1698860 4.886738e-14
HadCRUT4              0.06735899 1.189841e-12 0.1773596 6.525412e-15
HadCRUT4 (gap-filled) 0.06999475 1.193274e-14 0.1779717 6.317052e-15


##global average temperature anomalies
DF_gat<- data.frame(c(1880:2015),ANN_GISS,ANN_HD,ANN_LK)
colnames(DF_gat)<- c("Year","GISS anomalies","HadCRUT4","HadCRUT4(gap-filled using LK)")
DF_gat

#https://www.youtube.com/watch?v=-vSzKfqcTDg  (at time 5.47 similar situation of acf and pacf )
#http://stat565.cwick.co.nz/lectures/11-regression.pdf
#http://www.stats.uwo.ca/faculty/braun/ss359a/2004/notes/Chapter14/glsnotes.pdf

#############################
> DF
                       bhat_full       p_full  bhat_cut        p_cut
GISS                  0.07209890 6.037624e-11 0.1642531 5.397132e-14
HadCRUT4              0.06614553 3.844146e-13 0.1716232 5.425372e-14
HadCRUT4 (gap-filled) 0.06800860 1.125719e-15 0.1708448 4.144263e-14
