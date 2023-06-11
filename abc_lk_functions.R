#~~Functions
#######################empirical variogram###########################
####semivariance: r^(h)=1/(2Nh) sum_{i=1}^{Nh} (Y(xi+h)-Y(xi))^2 ####
####
####	X:2D spatial locations
####	y:observations at each of spatial locations in X
####    d.max: max distance
####    L: number of lag classes

vgm<- function(X,y,d.max,L){
N<- nrow(X)	

#Great circle distance matrix in kilometers
D_GC<- rdist.earth(X,miles=FALSE) 

#pairwise distances 
d_gc<- round(D_GC[t(upper.tri(D_GC))],4) 

#indices of corresponding pairs of distances computed above
pr <- rep.int(seq.int(N - 1), rev(seq.int(N - 1)))  #row index of distance matrix
pc <- rev(abs(sequence(seq.int(N - 1)) - N) + 1)    #column index of distance matrix

#Pairs and distances
P_D<- cbind(pr,pc,d_gc)

##Distance or lag classes
dist.max<- d.max          		 #max(d_gc)/2      #variogram not to be estimated beyond half 
							             #of the maximal distance (10019.15km ~ 10,000km) #take 7000km instead													
U<- seq(0,dist.max,length=L)
U1<- U[-length(U)]
U2<- U[-1]-0.0001
bin_c<- rep(0,length(U1))

		for(i in 1:length(U1)){
			bin_c[i]<- (U[i]+U[i+1])/2
		}		
	r_h<- numeric(0)  #semivariance
	l_h<- numeric(0)

	for(i in 1:length(U1)){
		ind<- which(d_gc>=U1[i] & d_gc<= U2[i])
		l_h[i]<- length(ind)

		z1<- y[P_D[ind,1]]
		z2<- y[P_D[ind,2]]	
		r_h[i]<- (1/(2*length(ind)))*sum((z1-z2)^2)
		}
DF.vgm<- data.frame(bin_c, l_h, r_h)
return(DF.vgm)
vgm
}

##check variogram function
#	library(LatticeKrig)
#	lon<- seq(-177.5, 177.5, by=5)
#	lat<- seq(-87.5, 87.5, by=5)
#	N<- length(lon)*length(lat)            #spatial points
#	X<- as.matrix((expand.grid(lon,lat)))  #2592...lat:-90 to 90
#	M<- 1
#	sigma<- sqrt(.01)
#	xr=cbind( c(-180, 180), c(-90,90))
#	LKinfo1<- LKrigSetup( xr,
#                      NC=3, LKGeometry="LKSphere", a.wght=6.5,
#                      nlevel=1, alpha=1, sigma=sigma, 
#                      choleskyMemory=list(nnzR= 4E6))
#	set.seed(224)
#	truef<- LKrig.sim(X, LKinfo=LKinfo1, M=M)
#	error<- sigma*matrix(rnorm(N*M),N,M)
#	y<- truef+error  #Simulated data with known	parameters sigma=sqrt(0.01), a.wght=6.5
#	vg<- vgm(X,y)
#	X11()
#	plot(vg[,1],vg[,3], xlab="distance (in km)", ylab="semivariance")

#LKrig.cov.plot function using distance in kms
lk_cov_sp<- function(LKinfo){
	NP<- 200
	grid.info <- LKinfo$latticeInfo$rangeLocations
	xlim <- grid.info[, 1]
	ylim <- grid.info[, 2]
	ux <- seq(xlim[1], xlim[2], , NP)
	uy <- seq(ylim[1], ylim[2], , NP)
	center <- rbind(c(ux[NP/2], uy[NP/2]))

	x1 <- cbind(ux, rep(center[2], NP))
	x2 <- rbind(center)
	d <- c(rdist.earth(x1, x2,miles=FALSE))
	y <- c(LKrig.cov(x1, x2, LKinfo))
	x1 <- cbind(rep(center[1], NP), uy)
	d2 <- c(rdist.earth(x1, x2,miles=FALSE))
	y2 <- c(LKrig.cov(x1, x2, LKinfo))
	return(list(d = cbind(d, d2), u = cbind(ux, uy), cov = cbind(y, y2), center = center, LKinfo = LKinfo))
}

#LKrig.cov.plot function using distance in radians

lk_cov_sp_r<- function(LKinfo){
	NP<- 200
	grid.info <- LKinfo$latticeInfo$rangeLocations
	xlim <- grid.info[, 1]
	ylim <- grid.info[, 2]
	ux <- seq(xlim[1], xlim[2], , NP)
	uy <- seq(ylim[1], ylim[2], , NP)
	center <- rbind(c(ux[NP/2], uy[NP/2]))

	x1 <- cbind(ux, rep(center[2], NP))
	x2 <- rbind(center)
	d <- c(rdist.earth(x1, x2, R=1))
	y <- c(LKrig.cov(x1, x2, LKinfo))
	x1 <- cbind(rep(center[1], NP), uy)
	d2 <- c(rdist.earth(x1, x2, R=1))
	y2 <- c(LKrig.cov(x1, x2, LKinfo))
	return(list(d = cbind(d, d2), u = cbind(ux, uy), cov = cbind(y, y2), center = center, LKinfo = LKinfo))
}

#LKrig.cov.plot function using distance in radians and with Y-transect only

lk_cov_sp_r<- function(LKinfo){
	NP<- 200
	grid.info <- LKinfo$latticeInfo$rangeLocations
	xlim <- grid.info[, 1]
	ylim <- grid.info[, 2]
	ux <- seq(xlim[1], xlim[2], , NP)
	uy <- seq(ylim[1], ylim[2], , NP)
	center <- rbind(c(ux[NP/2], uy[NP/2]))

	x1 <- cbind(ux, rep(center[2], NP))
	x2 <- rbind(center)
	d <- c(rdist.earth(x1, x2, R=1))
	y <- c(LKrig.cov(x1, x2, LKinfo))
	x1 <- cbind(rep(center[1], NP), uy)
	d2 <- c(rdist.earth(x1, x2, R=1))
	y2 <- c(LKrig.cov(x1, x2, LKinfo))
	return(list(d = cbind(d, d2), u = cbind(ux, uy), cov = cbind(y, y2), center = center, LKinfo = LKinfo))
}



##check lk_cov function
#	library(LatticeKrig)
#	lon<- seq(-177.5, 177.5, by=5)
#	lat<- seq(-87.5, 87.5, by=5)

#	N<- length(lon)*length(lat)            #spatial points
#	X<- as.matrix((expand.grid(lon,lat)))  #2592...lat:-90 to 90
#	M<- 1
#	sigma<- sqrt(.01)
#	lam<- sigma^2
#	xr=cbind( c(-180, 180), c(-90,90))
#	LKinfo1<- LKrigSetup( xr,
#                      NC=3, LKGeometry="LKSphere", a.wght=6.5,
#                      nlevel=1, alpha=1, lambda=lam, 
#                      choleskyMemory=list(nnzR= 4E6))
#   lkcv<- lk_cov_sp(LKinfo=LKinfo1)
#	matplot(lkcv$d,lkcv$cov,type="l",lty=c(1,2))
# X11()
# mapping = list(x=lon, y=lat,z = matrix(0,nrow=length(lon),ncol=length(lat)))
# image(mapping)
# points(rbind(lkcv$center),cex=3.5,  pch=4, col="red")
# points( cbind(lkcv$u[,1],rep(lkcv$center[2],length(lkcv$u[,1])) ))
# points(cbind(rep(lkcv$center[2],length(lkcv$u[,1])),lkcv$u[,2]))


#edited version of rdist.earth.vec function of fields package
#It at times gives error generating very small negative values (eg -1.54x10^-34) in 
#sqrt(1 - a) so round it to avoid warning 

rdist.earth.vec_ed<- function (x1, x2, miles = TRUE, R = NULL) 
{
    if (is.null(R)) {
        if (miles) 
            R = 3963.34
        else R = 6378.388
    }
    x1 = x1 * (pi/180)
    x2 = x2 * (pi/180)
    lonDist2 = (x2[, 1] - x1[, 1]) * (1/2)
    latDist2 = (x2[, 2] - x1[, 2]) * (1/2)
    a = sin(latDist2) * sin(latDist2) + cos(x1[, 2]) * cos(x2[, 
        2]) * sin(lonDist2) * sin(lonDist2)
    return(2 * atan2(sqrt(a), sqrt(round(1 - a,4))) * R)    ###.....use sqrt(round(1 - a,4))) instead of sqrt(1 - a)
}

#edited version of vgram function of fields package
# only difference is the use of rdist.earth.vec_ed function

vgram_ed<- function (loc, y, id = NULL, d = NULL, lon.lat = FALSE, dmax = NULL, 
    N = NULL, breaks = NULL, type = c("variogram", "covariogram", 
        "correlogram")) 
{
    type = match.arg(type)
    y <- cbind(y)
    if (is.null(id)) {
        n <- nrow(loc)
        is = rep(1:n, n)
        js = rep(1:n, rep(n, n))
        ind <- is > js
        id <- cbind(is, js)[ind, ]
    }
    if (is.null(d)) {
        loc <- as.matrix(loc)
        if (lon.lat) {
            d <- rdist.earth.vec_ed(loc[id[, 1], ], loc[id[, 2], 
                ])
        }
        else {
            d <- rdist.vec_ed(loc[id[, 1], ], loc[id[, 2], ])
        }
    }
    if (type == "correlogram") {
        sigma = apply(y, 2, sd, na.rm = TRUE)
        y = sweep(y, 2, (1/sigma), FUN = "*")
    }
    colMeans <- apply(y, 2, mean, na.rm = TRUE)
    yCntr = sweep(y, 2, colMeans)
    y1Cntr = yCntr[id[, 1], ]
    y2Cntr = yCntr[id[, 2], ]
    if (type == "variogram") {
        vg <- 0.5 * rowMeans(cbind((y1Cntr - y2Cntr)^2), na.rm = TRUE)
    }
    else {
        vg <- rowMeans(cbind(y1Cntr * y2Cntr), na.rm = TRUE)
    }
    call <- match.call()
    if (is.null(dmax)) {
        dmax <- max(d)
    }
    od <- order(d)
    d <- d[od]
    vg <- vg[od]
    ind <- d <= dmax & !is.na(vg)
    out <- list(d = d[ind], vgram = vg[ind], call = call, type = type)
    if (!is.null(breaks) | !is.null(N)) {
        out <- c(out, stats.bin(d[ind], vg[ind], N = N, breaks = breaks))
    }
    class(out) = c("vgram", class(out))
    out
}


