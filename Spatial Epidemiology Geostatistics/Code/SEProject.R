###### IMPORT PACKAGES #####
#install.packages("geodata")
#install.packages("geoR")
#install.packages("gstat")
library(geoR)
library(gstat)

##### IMPORT DATA ######
setwd("C:/Users/usuari/Desktop/Barcelona_courses/SpatialEpidimiology/Project_Spatial")
cls <- c(X="numeric", Y="numeric")
poly84=read.table("poly84.txt", header = TRUE)[1:2]
elevationsIslet=read.table("elevationsIslet.txt", header = TRUE)

##### DATA: POLY #####
head(poly84)
summary(poly84)
#geoscalg <- as.geodata(poly84)
head(poly84)
plot(poly84)
#plot(geoscalg)

##### DATA: ELEVATIONS ####
head(elevationsIslet)
summary(elevationsIslet)
geoscalelev <- as.geodata(elevationsIslet)
head(elevationsIslet[1:2])
plot(elevationsIslet)
plot(geoscalelev, borders = poly84)

##### Plots #####
points(geoscalelev, cex.min = 2, cex.max = 2, col = "gray")
plot(geoscalelev, borders = poly84)

geoscalelev2 <- as.geodata(elevationsIslet[-50,])
plot(geoscalelev2)

#Proposed trend 
lm1<-lm(elevationsIslet[1:2],data=elevationsIslet[-50,])
#lm(medv ~ lstat + I(lstat^2), data = train.data)
#model <- lm(noisy.y ~ poly(q,3))
#lm3 <- lm(elevationsIslet[1:2]  ~ poly(q, 2), data = elevationsIslet[-50,])
lm2<- lm(elevationsIslet[1:2] + elevationsIslet[1:2]^2, data=elevationsIslet[-50,])
summary(lm2)
#Remove the trend
residus<-round(residuals(lm1),digits=3)
elevationsIslet_res<-cbind(elevationsIslet[-50,],residus) 

geoscalelev_res<- as.geodata(elevationsIslet_res,data.col=4)
plot(geoscalelev_res)


maxdist<-max(dist(cbind(elevationsIslet_res$x,elevationsIslet_res$y)))


vario.c.classical<- variog(geodata=geoscalelev_res, option= "cloud", estimator.type="classical")
vario.bc.classical<-variog(geodata=geoscalelev_res, option = "bin", bin.cloud =TRUE, pairs.min=30, max.dist=maxdist/2, estimator.type = "classical");
vario.bc.classical$bins.lim
vario.bc.classical$ind.bin
vario.bc.classical$u
vario.bc.classical$v
vario.b.classical <- variog(geodata=geoscalelev_res, option = "bin", pairs.min=30, max.dist=maxdist/2,estimator.type = "classical")

par(mfrow=c(1, 3))
plot(vario.c.classical, main = "CLOUD", cex.main=1, cex.lab=1);
plot(vario.bc.classical, bin.cloud=TRUE,cex.lab=1,main = "\nBINNED BOXPLOTS", cex.main=1);
plot(vario.b.classical, main = "EMPIRICAL VARIOGRAMS (Classical)\nBINNED",cex.main=1, cex.lab=1, cex=1, pch=16)



#Robust variogram
vario.c.robust<-variog(geodata=geoscalelev_res,    option= "cloud", estimator.type="modulus")
vario.bc.robust<- variog(geodata=geoscalelev_res,  option = "bin", bin.cloud = T, pairs.min=30, max.dist=maxdist/2, estimator.type = "modulus");
vario.b.robust <- variog(geodata=geoscalelev_res, option = "bin", pairs.min=30, max.dist=maxdist/2,estimator.type = "modulus")

par(mfrow=c(1, 3))
plot(vario.c.robust, main = "\nCLOUD", cex.main=1, cex.lab=1);
plot(vario.bc.robust, bin.cloud=T, cex.lab=1);
title(main = "\nBINNED BOXPLOTS", cex.main=1 );
plot(vario.b.robust, main = "EMPIRICAL VARIOGRAMS (MODULUS)\nBINNED",cex.main=1, cex.lab=1, cex=1, pch=16)


#Variogrames direccionals:
###direction = c(0, pi/4, pi/2, 3*pi/4), tolerance = pi/8
vario.d <- variog4(geoscalelev_res, max.dist=maxdist/2,estimator.type="modulus")
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
indep.env<-variog.mc.env(geoscalelev_res,obj.variog=vario.b.robust)
indep.env<-variog.mc.env(geoscalelev_res,obj.variog=vario.b.robust,save.sim = TRUE)

plot(vario.b.robust, envelope = indep.env, main="CONFIDENCE BANDS FOR INDEPENDENT MODEL", lwd=2,   pch=16,cex.main=0.75)


