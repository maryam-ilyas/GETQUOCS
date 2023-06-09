##figure
rm(list=ls())
ff<- read.table("C:/Users/maryam/Dropbox/AMT_paper/mn.txt")$x
gg<- read.table("C:/Users/maryam/Dropbox/AMT_paper/ar5_com_sam.txt")
GG<- as.numeric(gg[137,])

X11()
hist(ff,col=rgb(1,0,0,0.5),prob=TRUE)
hist(ff[GG], col=rgb(0,0,1,0.5), add=TRUE,prob=TRUE)

##Distribution of Lahore
rm(list=ls())

op<- par(no.readonly=TRUE)
			par(oma=c(2,2,2,2))
			plot(1,1,type="n",xaxt="n",yaxt="n")
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

ff<- read.table("C:/Users/maryam/Dropbox/AMT_paper/lhr_full.txt")$x
gg<- read.table("C:/Users/maryam/Dropbox/AMT_paper/lhr_sampl.txt")$x

#median(ff);median(gg)
#quantile(ff, c(0.025,0.975))
#quantile(gg, c(0.025,0.975))

X11()
#par(op)
e<- paste("C:/Users/maryam/Dropbox/AMT_paper/lahore.png")
png(e,width = 5, height = 7, units = 'in', res = 300)
#par(oma=c(3,3,0.5,0.5),mar=c(1,1,0,0))
par(mar=c(4,4,0,0))

hist(ff,col=rgb(1,0,0,0.5),prob=TRUE,main="", xlab="temperature anomaly (°C)", ylab="probability", ylim=c(0,0.35))
hist(gg, col=rgb(0,0,1,0.5), add=TRUE,prob=TRUE)

legend("topleft",bty="n", legend=c("Full ensemble", "Subsampled ensemble"), 
       col=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)),lty=c(1,1),lwd=c(3,3))
dev.off()