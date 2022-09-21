####EXPLORE DATA ANALYSIS II:
##EXLPORE THE SAMLL -SCALE DEPENDENCE

install.packages("geoR")
install.packages("gstat")
library(geoR)
library(gstat)

#1. Coal data

#The coal seam on the Rodena Mine Property in Greence County, Pennsylvania. 
#The data frame contains 208 coal 
#core samples collectes on a grid given by x and y planar coordinates. 

data(coalash)
summary(coalash)
head(coalash)

#Converts an Object to the Class "geodata"#
geocoal <- as.geodata(coalash)
points(geocoal, cex.min = 2, cex.max = 2, col = "gray")
plot(geocoal)


#Remove outlier, data in row number 50.
geocoal2 <- as.geodata(coalash[-50,])
plot(geocoal2)

#Proposed trend 
lm1<-lm(coalash~x,data=coalash[-50,])
summary(lm1)
#Remove the trend
residus<-round(residuals(lm1),digits=3)
coalash_res<-cbind(coalash[-50,],residus) 

geocoal_res<- as.geodata(coalash_res,data.col=4)
plot(geocoal_res)



###
##Variog: Computes sample (empirical) variograms with options 
#for the classical or robust estimators. 
#Output can be returned as a binned variogram, a variogram cloud 
#or a smoothed variogram.
# option = c("bin", "cloud", "smooth")
# estimator.type = c("classical", "modulus")
#Trend: option default cte  
  #trend="1st" (trend~coords[,1]+coords[,2]) 
  #or trend=2nd"  trend~poly(coords[,1]+coords[,2],2)
  #quadratic form on the x-coordinate trend~ coords[,2]+poly(coords[,1],2)
# uvec= a vector with values used to define variogram binning
  #default uvec= seq(0, maxdist, l = 13)

maxdist<-max(dist(cbind(coalash_res$x,coalash_res$y)))

vario.c.classical<- variog(geodata=geocoal_res, option= "cloud", estimator.type="classical")
vario.bc.classical<-variog(geodata=geocoal_res, option = "bin", bin.cloud =TRUE, pairs.min=30, max.dist=maxdist/2, estimator.type = "classical");
vario.bc.classical$bins.lim
vario.bc.classical$ind.bin
vario.bc.classical$u  
vario.bc.classical$v  # the distances 
vario.b.classical <- variog(geodata=geocoal_res, option = "bin", pairs.min=30, max.dist=maxdist/2,estimator.type = "classical")

par(mfrow=c(1, 3))
plot(vario.c.classical, main = "CLOUD", cex.main=1, cex.lab=1);
plot(vario.bc.classical, bin.cloud=TRUE,cex.lab=1,main = "\nBINNED BOXPLOTS", cex.main=1);
plot(vario.b.classical, main = "EMPIRICAL VARIOGRAMS (Classical)\nBINNED",cex.main=1, cex.lab=1, cex=1, pch=16);

#Robust variogram
vario.c.robust<-variog(geodata=geocoal_res,    option= "cloud", estimator.type="modulus")
vario.bc.robust<- variog(geodata=geocoal_res,  option = "bin", bin.cloud = T, pairs.min=30, max.dist=maxdist/2, estimator.type = "modulus");
vario.b.robust <- variog(geodata=geocoal_res, option = "bin", pairs.min=30, max.dist=maxdist/2,estimator.type = "modulus")

par(mfrow=c(1, 3))
plot(vario.c.robust, main = "\nCLOUD", cex.main=1, cex.lab=1);
plot(vario.bc.robust, bin.cloud=T, cex.lab=1);
title(main = "\nBINNED BOXPLOTS", cex.main=1 );
plot(vario.b.robust, main = "EMPIRICAL VARIOGRAMS (MODULUS)\nBINNED",cex.main=1, cex.lab=1, cex=1, pch=16);

#Variogrames direccionals:
###direction = c(0, pi/4, pi/2, 3*pi/4), tolerance = pi/8
vario.d <- variog4(geocoal_res, max.dist=maxdist/2,estimator.type="modulus")  # computes the variogram in 4 driections or 4 degrees.
par(mfrow=c(1, 1))
windows()
plot(vario.d,lwd =2,legend=FALSE)
legend(x="bottomright", inset=0.01, lty=c(1,2,3,4), col=c("black", "red", "green","blue"),
       legend=c("0º", "45º", "90º","135º"), cex=1);

title(main="DIRECTIONAL EMPIRICAL VARIOGRAM", cex.main=1);


#variog.mc.env: Computes envelops for empirical variograms by permutation 
#of the data values on the spatial locations.

#variog.mc.env(geodata, coords = geodata$coords, data = geodata$data,
#              obj.variog, nsim = 99, save.sim = FALSE, messages) #
par(mfrow=c(1, 1))
indep.env<-variog.mc.env(geocoal_res,obj.variog=vario.b.robust)
indep.env<-variog.mc.env(geocoal_res,obj.variog=vario.b.robust,save.sim = TRUE)

plot(vario.b.robust, envelope = indep.env, main="CONFIDENCE BANDS FOR INDEPENDENT MODEL", lwd=2,   pch=16,cex.main=0.75)

###########################################################################
########SCALLOPS
###############################
setwd("D:/Master d'Estadistica i IO/Geoestadística/Diapositives curs/Capitol 2")

scallops<- read.table("Scallops_R.txt",head=TRUE,sep=" ",dec=".");
scallops$lgcatch <- log(scallops$tcatch + 1);
scallops$long <- -scallops$long;

# axes rotation
rotate.axis <- function(xy,theta=0){
  # xy is a nx2 matrix of coordinates
  # theta is the angle of rotation desireed in degrees
  # theta can be negative
  # XY is the new coordinates
  pimult <- 2*pi*theta/360
  newx <- c(cos(pimult),sin(pimult))
  newy <- c(-sin(pimult),cos(pimult))
  XY <- as.matrix(xy) %*% cbind(newx,newy)
  as.data.frame(XY)
}
xy <- scallops[c("long","lat")]
scall.rot <- cbind(rotate.axis(xy,52),lgcatch=scallops$lgcatch)
summary(as.geodata(scall.rot))
geoscal.rot<-as.geodata(scall.rot)
windows()
plot(geoscal.rot)
with(geoscal.rot , plot(coords[, 1], data, xlab = "W-E",ylab = "elevation data", pch = 20, cex = 0.7))
lines(with(geoscal.rot, lowess(data ~ coords[, 1])))
with(geoscal.rot , plot(coords[, 2], data, xlab = "S-N",ylab = "elevation data", pch = 20, cex = 0.7))
lines(with(geoscal.rot, lowess(data ~ coords[, 2])))

lm.scp.rot<-lm(lgcatch ~poly(newx,1)+poly(newy,3),data=scall.rot)
summary(lm.scp.rot)


#Sampling variogram  lgcath original data. This is not correct, not run because of the trend
#of the data. And here, it is not remove.

maxdist<-max(dist(cbind(scall.rot$newx,scall.rot$newy)))
maxdist
cloud <- variog(geoscal.rot, option = "cloud", estimator.type = "modulus");
bp <- variog(geoscal.rot, option = "bin", bin.cloud = T,
             pairs.min=30, max.dist=maxdist/2,
             estimator.type = "modulus");
bin <- variog(geoscal.rot, option = "bin", pairs.min=30, max.dist=maxdist,
              estimator.type = "modulus");

windows()
par(mfrow=c(1,3))

plot(cloud,main="CLOUD", cex.main=1, cex.lab=1, cex=2);
plot(bp, bin.cloud=T, cex.lab=1);
title(main = "BINNED BOXPLOTS", cex.main=1);
plot(bin,main="EMPIRICAL VARIOGRAMS (MODULUS)\nBINNED",cex.main=1, cex.lab=1, cex=2, pch=16);



# Samppling variogram: Residus

residus.scallops<-round(residuals(lm.scp.rot),digits=3)
scall.rot.residus <- cbind(scall.rot,residus.scallops)
geoscal.rot.residus<-as.geodata(scall.rot.residus,data.col=4)
plot(geoscal.rot.residus)

cloud.2 <- variog(geoscal.rot.residus, option = "cloud", estimator.type = "modulus");
bp.2 <- variog(geoscal.rot.residus, option = "bin", bin.cloud = T,
               pairs.min=30, max.dist=maxdist/2,
               estimator.type = "modulus");
bin.2 <- variog(geoscal.rot.residus, option = "bin", pairs.min=25, max.dist=maxdist/2,
                estimator.type = "modulus");

windows()
par(mfrow=c(1,3))

plot(cloud.2,main="CLOUD", cex.main=1, cex.lab=1, cex=2);
plot(bp.2, bin.cloud=T, cex.lab=1);
title(main = "BINNED BOXPLOTS", cex.main=1);
plot(bin.2,main="EMPIRICAL VARIOGRAMS (MODULUS)\nBINNED",cex.main=1, cex.lab=1, cex=2, pch=16);


par(mfrow=c(1, 2))
plot(bin,main="EMPIRICAL VARIOGRAMS (MODULUS)\nBINNED",cex.main=1, cex.lab=1, cex=2, pch=16);
plot(bin.2,main="EMPIRICAL VARIOGRAMS",cex.main=1, cex.lab=1, cex=2, pch=16);


# check isotropy: 
##Residus
par(mfrow=c(1,1))
variod <- variog4(geoscal.rot.residus,max.dist=maxdist/2, pairs.min=30,estimator.type = "modulus");
plot(variod,lyt=2,legend=FALSE)
legend(x="bottomright", inset=0.01, lty=c(1,2,3,4), col=c("black", "red", "green","blue"),
       legend=c("0º", "45º", "90º","135º"), cex=0.5);

#variog.mc.env: Computes envelops for empirical variograms by permutation 
#of the data values on the spatial locations.

#variog.mc.env(geodata, coords = geodata$coords, data = geodata$data,
#              obj.variog, nsim = 99, save.sim = FALSE, messages) #
set.seed(100)
windows()
par(mfrow=c(1,1))
indep.env<-variog.mc.env(geoscal.rot.residus,coords=geoscal.rot.residus$coords, data=geoscal.rot.residus$data,obj.variog=bp.2,nsim=200)
plot(bp.2, envelope = indep.env, main="CONFIDENCE BANDS FOR INDEPENDENT MODEL", lwd=2,     pch=16)


#######################Study the Anisotropy of the process
source("RoseDiagram.R")
rose.diagram(data.var=geoscal.rot.residus$data,data.cds=geoscal.rot.residus$coord,max.dist=maxdist/2,numcases=10,numdirec=4,poly.tnd="cte",crit.val=5)

rose.diagram(data.var=geoscal.rot.residus$data,data.cds=geoscal.rot.residus$coord,max.dist=maxdist/2,numdirec=8,poly.tnd="cte",crit.val=5)


angle<-pi/2
ratio<-3

#Revissar aquesta funció
coords.aniso(geoscal.rot.residus$coords, c(angle,ratio), reverse = FALSE)
scall.ani<-cbind(coords.aniso(geoscal.rot.residus$coords, c(angle,ratio),reverse = FALSE),geoscal.rot$data)
geoscal.ani<-as.geodata(scall.ani)
variod <- variog4(geoscal.ani, max.dist=maxdist/2, pairs.min=30,estimator.type = "classical");
plot(variod,lyt=2,legend=FALSE)
legend(x="bottomright", inset=0.01, lty=c(1,2,3,4), col=c("black", "red", "green","blue"),
       legend=c("0º", "45º", "90º","135º"), cex=1)