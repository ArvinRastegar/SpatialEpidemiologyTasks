###### IMPORT PACKAGES #####
#install.packages("geodata")
#install.packages("geoR")
#install.packages("gstat")
library(geoR)
library(gstat)
library(sp)
 
##### IMPORT DATA ######
#setwd("G:/My Drive/UNI/UPC/Spatial Epidemiology/R")
cls <- c(X="numeric", Y="numeric")
poly84=read.table("poly84.txt", header = TRUE)[1:2]
elevationsIslet=read.table("elevationsIslet.txt", header = TRUE)

##### DATA: POLY #####
head(poly84)
summary(poly84)
x_coord <- poly84$X
y_coord <- poly84$Y
xym <- cbind(x_coord, y_coord)
xym <- rbind(xym, cbind(x_coord[1], y_coord[1])) # add first value to close polygon

# convert to polygon object
#p = Polygon(xym)
#ps = Polygons(list(p),1)
#sps = SpatialPolygons(list(ps))
#plot(sps)

# DATA: ELEVATIONS #####
# Arvin #### 
#setwd("C:/Users/usuari/Desktop/Barcelona_courses/SpatialEpidimiology/Project_Spatial")
head(elevationsIslet)
summary(elevationsIslet)
geoscalelev <- as.geodata(elevationsIslet)
head(elevationsIslet[1:2])
plot(elevationsIslet)
plot(geoscalelev, borders = poly84)
##### Plots #####
points(geoscalelev, cex.min = 2, cex.max = 2, col = "gray")
plot(geoscalelev, borders = poly84)

geoscalelev2 <- as.geodata(elevationsIslet)
plot(geoscalelev2)

#Proposed trend 
###here
lm1<-lm(formula = data ~ x+y, data=elevationsIslet)
summary(lm1)

#Remove the trend
residus<-round(residuals(lm1),digits=3)
elevationsIslet_res<-cbind(elevationsIslet, residus) 

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


# Rodrigo ###
head(elevationsIslet)
summary(elevationsIslet)

elevations.geo <- as.geodata(elevationsIslet, coords.col=1:2, data.col = 3)
plot(elevations.geo)

# convert to polar coordinates IGNORE ####
# move center
x.new <- elevationsIslet$x - 70
y.new <- elevationsIslet$y - 27
elevations.centered <- data.frame(list(x=x.new,y=y.new,data=elevationsIslet$data))
plot(as.geodata(elevations.centered))

r <- sqrt(x.new**2 + y.new**2)
theta <- atan(y.new / x.new)
elevations.polar <- data.frame(list(x=r, y=theta, data=elevationsIslet$data))
elevations.geo.polar <- as.geodata(elevations.polar)
plot(elevations.geo.polar)

lm1<-lm(data~x,data=elevations.geo.polar)
summary(lm1)

# Calculate residuals for envelopes ####
lm.elevation <- lm(data~x+y, data=elevationsIslet)
summary(lm.elevation)
residus.elevation<-round(residuals(lm.elevation),digits=3)
elevation_res <- cbind(elevationsIslet, residus.elevation)
geo_elev_res <- as.geodata(elevation_res, data.col=4)
plot(geo_elev_res)

# variogram
set.seed(1000)
maxdist<-max(dist(cbind(elevation_res$x,elevation_res$y)))
vario.b.robust <- variog(geodata=geo_elev_res, option = "bin", pairs.min=30, max.dist=maxdist/2,estimator.type = "modulus")
plot(vario.b.robust, main = "EMPIRICAL VARIOGRAMS (MODULUS)\nBINNED",cex.main=1, cex.lab=1, cex=1, pch=16)


# envelopes
indep.env<-variog.mc.env(geo_elev_res,obj.variog=vario.b.robust)
indep.env<-variog.mc.env(geo_elev_res,obj.variog=vario.b.robust,save.sim = TRUE)

plot(vario.b.robust, envelope = indep.env, main="CONFIDENCE BANDS FOR INDEPENDENT MODEL", lwd=2,   pch=16,cex.main=0.75)

# Theoretical variograms ####

maxdist <- variog(geo_elev_res, option = "cloud")$max.dist;

bin <- variog(geo_elev_res,option = "bin", pairs.min=30, max.dist=maxdist/2,
              estimator.type = "modulus");
plot(bin,main="Empirical variogram",cex.main=1, cex.lab=1, cex=2, pch=16);
par(mfrow=c(1, 1))
windows()
eyefit(bin, silent = FALSE)

###Exponential  
#Theoretical variogram nugget fixed at 0.4
#wls1 <- variofit(bin, cov.model = "exponential", ini = c(3,0.3),fix.nugget = TRUE, nugget =0.4 , weights="cressie");

#Theoretical variogram: Estimate the nugget parameter
wls1 <- variofit(bin, ini.cov.pars = c(125,24.45), cov.model = "exponential",
                 fix.nugget = F, nugget = 47 , weights="cressie");

wls1
summary(wls1)$sum.of.squares
#Gausssian variogram

wls2 <- variofit(bin, ini.cov.pars = c(104, 29.74), cov.model = "gaussian",
                 fix.nugget = F, nugget =67.88, weights="cressie"); 

wls2
summary(wls2)$sum.of.squares

#Spherical
wls3 <- variofit(bin, cov.model = "spherical", ini.cov.pars=c(130, 42),
                 fix.nugget = F, nugget =36.55, weights="cressie"); 

wls3
summary(wls3)$sum.of.squares

#Matern variogram
#Be careful  fix.kappa=TRUE (default)  kappa=0.5 (default)

wls4 <- variofit(bin, cov.model = "matern", ini.cov.pars = c(154.5, 24.45),
                 fix.nugget = F, nugget =19.32, fix.kappa = FALSE, kappa=0.54,weights="cressie"); 

wls4
summary(wls4)$sum.of.squares


windows()
plot(bin,  main="PARAMETRIC VARIOGRAMS",cex.main=1, pch=16);

lines(wls1, lwd = 2, col="red", max.dist=maxdist/2);
lines(wls2, lwd = 2, col="blue", max.dist=maxdist/2);
lines(wls3, lwd = 2, col="green3", max.dist=maxdist/2);
lines(wls4, lwd = 2, col="yellow", max.dist=maxdist/2);
legend(x="bottomright", inset=0.01, lty=c(1, 1), col=c("red", "blue", "green3","yellow"),
       legend=c("Exponetial", "Gaussian", "Spherical","Matern"), cex=1)

# the best variograms are 1) matern (wls4) and 2) exponential (wls1)
#Preditions locations
# Arvin ###
rnx <- range(elevations.geo$coords[,1]);
rny <- range(elevations.geo$coords[,2]);
newx.grid <- seq(rnx[1],rnx[2],l=51);
newy.grid <- seq(rny[1],rny[2],l=51);
dsgr.grid <- expand.grid(newx=newx.grid, newy=newy.grid);
plot(dsgr.grid)

kc1<-krige.conv(elevations.geo,coords= elevations.geo$coords,
                data=elevations.geo$data,locations= dsgr.grid,
                krige =krige.control(type.krige="OK",obj.m = wls1,trend.l= ~coords[,1]+poly(coords[,2],3), 
                                     trend.d=~coords[,1]+poly(coords[,2],3)))
attributes(kc1)

image(kc1, main="Universal kriging estimates")
image(kc1, val=sqrt(kc1$krige.var), main="Universal kriging std. errors",col=gray(seq(1,0.1,l=30)))

contour(kc1,filled=TRUE,coords.data=elevations.geo$coords,color=heat.colors)
contour(kc1,val=sqrt(kc1$krige.var),filled=TRUE,coords.data=elevations.geo$coords,color=cm.colors)

wls1.1 <- likfit(geo_elev_res, cov.model = "exponential", ini =c(3,0.3),
                 fix.nugget = F, nugget =47 ,lik.method = "REML")

xv <- xvalid(elevations.geo, model =wls1.1, reestimate = F);
names(xv)

VC1 <- mean(xv$error/sqrt(xv$krige.var));
VC2 <- sqrt(mean((xv$error/sqrt(xv$krige.var))^2));
VC3 <- sqrt(mean(xv$error^2));
VC1;VC2;VC3

#Variogram matern wls4 

kc2<-krige.conv(elevations.geo,coords= elevations.geo$coords,
                data=elevations.geo$data,locations= dsgr.grid,
                krige =krige.control(type.krige="OK",obj.m = wls4,trend.l= ~coords[,1]+poly(coords[,2],3), 
                                     trend.d=~coords[,1]+poly(coords[,2],3)))
attributes(kc2)

image(kc2, main="Universal kriging estimates")
image(kc2, val=sqrt(kc2$krige.var), main="Universal kriging std. errors",col=gray(seq(1,0.1,l=30)))

contour(kc2,filled=TRUE,coords.data=elevations.geo$coords,color=heat.colors)
contour(kc2,val=sqrt(kc2$krige.var),filled=TRUE,coords.data=elevations.geo$coords,color=cm.colors)

wls4.1 <- likfit(geo_elev_res, cov.model = "matern", ini =c(3,0.3),
                 fix.nugget = F, nugget =47, fix.kappa =F,
                 kappa = 0.54,lik.method = "REML")

xv.2 <- xvalid(elevations.geo, model =lk2, reestimate = F);
names(xv.2)

VC1.2 <- mean(xv.2$error/sqrt(xv.2$krige.var));
VC2.2 <- sqrt(mean((xv.2$error/sqrt(xv.2$krige.var))^2));
VC3.2 <- sqrt(mean(xv.2$error^2));
VC1.2;VC2.2;VC3.2
# Show the predictions and their standard errors 
## Exponential
par(mfrow=c(1, 1))
image(kc1, main="Universal kriging estimates")
image(kc1, val=sqrt(kc1$krige.var), main="Universal kriging std. errors",col=gray(seq(1,0.1,l=30)))
## matern 
image(kc2, main="Universal kriging estimates")
image(kc2, val=sqrt(kc2$krige.var), main="Universal kriging std. errors",col=gray(seq(1,0.1,l=30)))
