#install.packages("geoR")
library(geoR)

# setwd("D:/Master d'Estadistica i IO/Geoestadística")
# setwd("d:\\Docencia");getwd()


#############################################################
# Datos: Scallops.
##############################################################

scallops<- read.table("Scallops_R.txt",head=TRUE,sep=" ",dec=".");
scallops$lgcatch <- log(scallops$tcatch + 1);
scallops$long <- -scallops$long;
geoscalg <- as.geodata(scallops[,-c(1,2)],coords.col = 2:1,data.col=6)

# rotacion de ejes
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
plot(geoscal.rot)


lm.scp.rot<-glm(lgcatch ~poly(newx,1)+poly(newy,3),data=scall.rot)
scall.rot[,"residuals"]<-residuals(lm.scp.rot)
scall.rot.residus<-cbind(rotate.axis(xy,52),residus=lm.scp.rot$residuals)
geoscalresidus.rot<-as.geodata(scall.rot.residus)
plot(geoscalresidus.rot)



# Empirical variogram (residuals) 
summary(geoscalresidus.rot)
maxdist <- variog(geoscalresidus.rot, option = "cloud")$max.dist;

bin <- variog(geoscalresidus.rot,option = "bin", pairs.min=30, max.dist=maxdist/2,
              estimator.type = "modulus");


par(mfrow=c(1, 2))
plot(bin,main="Empirical variogram",cex.main=1, cex.lab=1, cex=2, pch=16);
#title(main="Empirical variogram")


###################################################################################
#Estimate variogram Weighted or ordinary least squares

#variofit(vario, ini.cov.pars, cov.model,
#         fix.nugget = FALSE, nugget = 0,
#         fix.kappa = TRUE, kappa = 0.5,
#         simul.number = NULL, max.dist = vario$max.dist,
#         weights, minimisation.function,
#         limits = pars.limits(), messages, ...)

#ini.cov.pars  
#initial values for the covariance parameters: 
#sigma^2 (partial sill) and phi (range parameter)
# covariance model, you can choose: "exponential", "matern", "gaussian", 
#  "spherical",  "wave", "powered.exponential"

#value of the nugget parameter tau^2. 
#This is an estimate if fix.nugget = FALSE  (default) otherwise, a fixed value. 

# you do not need to initizalize the nugget only the partial sill
# and range, therefore ini=c(.5,.3) where .5 is the initial value
# for partial sill and .3 is the initial value for range.
#
# if you choose the mattern then the smoothing parameter is kappa
# here kappa is .5 (the default)
#
# weights: are used to define the loss function to be 
#minimised: npairs, cressie, equal 

par(mfrow=c(1, 1))
windows()
eyefit(bin, silent = FALSE)


###Exponential  
#Theoretical variogram nugget fixed at 0.4
wls1 <- variofit(bin, cov.model = "exponential", ini = c(3,0.3),
                 fix.nugget = TRUE, nugget =0.4 , weights="cressie");
#Theoretical variogram: Estimate the nugget parameter
wls1 <- variofit(bin, cov.model = "exponential", ini = c(3,0.3),
                 fix.nugget = F, nugget =0.4 , weights="cressie");

wls1
summary(wls1)$sum.of.squares
#Gausssian variogram

wls2 <- variofit(bin, cov.model = "gaussian", ini = c(3,0.3),
                 fix.nugget = F, nugget =0.4, weights="cressie"); 

wls2
summary(wls2)$sum.of.squares

#Spherical
wls3 <- variofit(bin, cov.model = "spherical", ini = c(3,0.3),
                 fix.nugget = F, nugget =0.4, weights="cressie"); 

wls3
summary(wls3)$sum.of.squares

#Matern variogram
#Be careful  fix.kappa=TRUE (default)  kappa=0.5 (default)
 
wls4 <- variofit(bin, cov.model = "matern", ini = c(3,0.3),
                 fix.nugget = F, nugget =0.4, fix.kappa = FALSE, kappa=0.5,weights="cressie"); 

wls4
summary(wls4)$sum.of.squares


windows()
plot(bin,  main="PARAMETRIC VARIOGRAMS",cex.main=1, pch=16);

lines(wls1, lwd = 2, col="red", max.dist=maxdist/2);
lines(wls2, lwd = 2, col="blue", max.dist=maxdist/2);
lines(wls3, lwd = 2, col="green3", max.dist=maxdist/2);
lines(wls4, lwd = 2, col="yellow", max.dist=maxdist/2);
legend(x="bottomright", inset=0.01, lty=c(1, 1), col=c("red", "blue", "green3","yellow"),
       legend=c("Exponetial", "Gaussian", "Spherical","Matern"), cex=1);


### Independent variogram 
set.seed(1)
indep.env<-variog.mc.env(geoscalresidus.rot,coords=geoscalresidus.rot$coords, data=geoscalresidus.rot$data,obj.variog=bin)
plot(bin, envelope = indep.env, main="CONFIDENCE BANDS FOR INDEPENDENT MODEL", lwd=2,pch=16)

#variog.model.env:Envelops for Empirical Variograms Based on Model Parameters
#Computes envelopes for a empirical variogram by simulating data for given model parameters
set.seed(10)
par(mfrow=c(1, 4))
env <- variog.model.env(geoscalresidus.rot,model.pars = wls1,obj.variog=bin,nsim = 999);
plot(bin,main="CONFIDENCE BANDS FOR EMPIRICAL VARIOGRAM \nMODEL:exponential", lwd=2,pch=16, envelope = env);

env <- variog.model.env(geoscalresidus.rot,obj.variog=bin, model.pars = wls2,nsim = 999);
plot(bin,main="CONFIDENCE BANDS FOR EMPIRICAL VARIOGRAM\nMODEL:Gaussian", lwd=2,pch=16, envelope = env);

env <- variog.model.env(geoscalresidus.rot,obj.variog=bin, model.pars = wls3,nsim = 999);
plot(bin,main="CONFIDENCE BANDS FOR EMPIRICAL VARIOGRAM\nMODEL:Spherical", lwd=2,pch=16, envelope = env);
                        
env <- variog.model.env(geoscalresidus.rot,obj.variog=bin, model.pars = wls4,nsim = 999);
plot(bin,main="CONFIDENCE BANDS FOR EMPIRICAL VARIOGRAM\nMODEL:Matern", lwd=2,pch=16, envelope = env);



#Likelihood Based Parameter Estimation for Gaussian Random Fields
#use likfit
#ini to initialize parameters
# partial sill and range
# kappa is power for powered exponential model
# kappa is smoothing for Matern
#method is ML or REML:default ML
#trend is cte, but you can also choose a linear
# trend by saying trend=1 or quadratic trend=2


lk1 <- likfit(geoscalresidus.rot ,cov.model = "exponential", ini =c(3,0.3),
              fix.nugget = F, nugget =0.4 ,lik.method = "REML");

lk1
summary(lk1)
lk1$AIC
lk1$sigmasq
lk1$tausq
lk1$phi
lk1$practicalRange
            
                        
lk2 <- likfit(geoscalresidus.rot,cov.model = "gaussian", ini = c(3,0.3),
              fix.nugget = F, nugget =0.4 ,lik.method = "REML");

summary(lk2)
lk2$AIC
lk2$sigmasq
lk2$tausq
lk2$phi
lk2$practicalRange

lk3 <- likfit(geoscalresidus.rot, cov.model = "spherical", ini = c(3,0.3),fix.nugget = F
              , nugget =0.4,lik.method = "REML");
summary(lk3)

lk4 <- likfit(geoscalresidus.rot, cov.model = "matern", ini = c(3,0.3),
              fix.nugget = F, nugget =0.4 ,fix.kappa = FALSE, kappa=1,lik.method = "REML");

summary(lk4)
lk4

#Expontial
lk1$AIC;
#Gaussian
lk2$AIC;
#Spherical
lk3$AIC;
#matern
lk4$AIC
# Lowest AIC gives us the best results. 
bin <- variog(geoscalresidus.rot, option = "bin", pairs.min=25, max.dist=maxdist,
              estimator.type = "modulus", uvec = seq(0, maxdist/2, l = 13))
windows()
par(mfrow=c(1, 1))
plot(bin,  main="PARAMETRIC VARIOGRAMS",cex.main=1, pch=16 ,ylim=c(0,8));
lines(lk1, lwd = 2, col="red", max.dist=maxdist/2);
lines(lk2, lwd = 2, col="blue", max.dist=maxdist/2);
lines(lk3, lwd = 2, col="green3", max.dist=maxdist/2);
lines(lk4, lwd = 2, col="yellow", max.dist=maxdist/2);
legend(x="bottomright", inset=0.01, lty=c(1, 1), col=c("red", "blue", "green","yellow"),
       legend=c("Exponetial", "Gaussian", "Spherical","Matern"), cex=1);

                      
 ### Variabilidad
set.seed(10)
indep.env<-variog.mc.env(geoscalresidus.rot,coords=geoscalresidus.rot$coords, data=geoscalresidus.rot$data,obj.variog=bin)
plot(bin, envelope = indep.env, main="CONFIDENCE BANDS FOR INDEPENDENT MODEL", lwd=2,pch=16)
 

#variog.model.env:Envelops for Empirical Variograms Based on Model Parameters
par(mfrow=c(1, 4))
env <- variog.model.env(geoscalresidus.rot,coords=geoscalresidus.rot$coords,obj.variog=bin, model.pars = lk1);
plot(bin,main="CONFIDENCE BANDS FOR EMPIRICAL VARIOGRAM\nMODEL:Exponential", lwd=2,pch=16, envelope = env);
                        
env <- variog.model.env(geoscalresidus.rot,coords=geoscalresidus.rot$coords,obj.variog=bin, model.pars = lk2);
plot(bin,main="CONFIDENCE BANDS FOR EMPIRICAL VARIOGRAM\nMODEL:Gaussian", lwd=2,pch=16, envelope = env, cex.main=1);

env <- variog.model.env(geoscalresidus.rot,coords=geoscalresidus.rot$coords,obj.variog=bin, model.pars = lk3);
plot(bin,main="CONFIDENCE BANDS FOR EMPIRICAL VARIOGRAM\nMODEL:Spherical", lwd=2,pch=16, envelope = env, cex.main=1);
                        
env <- variog.model.env(geoscalresidus.rot,coords=geoscalresidus.rot$coords,obj.variog=bin, model.pars = lk4);
plot(bin,main="CONFIDENCE BANDS FOR EMPIRICAL VARIOGRAM\nMODEL:Matern", lwd=2,pch=16, envelope = env, cex.main=1);


################################################.
#Increase the number of bins to estimate the variogram  
#############################################

bin <- variog(geoscalresidus.rot,option = "bin", pairs.min=30, max.dist=maxdist/2,
              estimator.type = "modulus",uvec=seq(0,maxdist/2,l=30));

#Return to line 92


#########################################################################
#####Working original data, and estimate the estimate the trend 

bin2 <- variog(geoscal.rot,option = "bin", pairs.min=30, max.dist=maxdist/2,
              estimator.type = "modulus",uvec=seq(0,maxdist/2,l=30),
              trend = ~geoscal.rot$coords[,1]+poly(geoscal.rot$coords[,2],3));
bin2
plot(bin2,main="Empirical variogram",cex.main=1, cex.lab=1, cex=2, pch=16);

#Not estimated the trend parametres

wls1.2 <- variofit(bin2, cov.model = "exponential", ini = c(3,0.3), 
                   fix.nugget = FALSE, nugget =0.4 , weights="cressie");

summary(wls1.2)$sum.of.squares

wls2.2 <- variofit(bin2, cov.model = "gaussian", ini = c(3,0.3),
                 fix.nugget = F, nugget =0.4, weights="cressie"); 

wls2.2
summary(wls2.2)$sum.of.squares

#Spherical
wls3.2 <- variofit(bin2, cov.model = "spherical", ini = c(3,0.3),
                 fix.nugget = F, nugget =0.4, weights="cressie"); 

wls3.2
summary(wls3.2)$sum.of.squares

#Matern variogram
wls4.2 <- variofit(bin2, cov.model = "matern", ini = c(3,0.3),
                 fix.nugget = F, nugget =0.4, fix.kappa = FALSE, kappa=0.5,weights="cressie"); 

wls4.2
summary(wls4.2)$sum.of.squares


plot(bin2,  main="PARAMETRIC VARIOGRAMS",cex.main=1, pch=16);

lines(wls1.2, lwd = 2, col="red", max.dist=maxdist/2);
lines(wls2.2, lwd = 2, col="blue", max.dist=maxdist/2);
lines(wls3.2, lwd = 2, col="green3", max.dist=maxdist/2);
lines(wls4.2, lwd = 2, col="yellow", max.dist=maxdist/2);
legend(x="bottomright", inset=0.01, lty=c(1, 1), col=c("red", "blue", "green3","yellow"),
       legend=c("Exponetial", "Gaussian", "Spherical","Matern"), cex=1);

set.seed(10)
par(mfrow=c(1, 4))

env <- variog.model.env(geoscal.rot,obj.variog=bin2, model.pars = wls1.2,nsim = 999);
plot(bin2,main="CONFIDENCE BANDS FOR EMPIRICAL VARIOGRAM\nMODEL:Exponential", lwd=2,pch=16, envelope = env);

env <- variog.model.env(geoscal.rot,obj.variog=bin2, model.pars = wls2.2,nsim = 999);
plot(bin2,main="CONFIDENCE BANDS FOR EMPIRICAL VARIOGRAM\nMODEL:Gaussian", lwd=2,pch=16, envelope = env);

env <- variog.model.env(geoscal.rot,obj.variog=bin2, model.pars = wls3.2,nsim = 999);
plot(bin2,main="CONFIDENCE BANDS FOR EMPIRICAL VARIOGRAM\nMODEL:Spherical", lwd=2,pch=16, envelope = env);

env <- variog.model.env(geoscal.rot,obj.variog=bin2, model.pars = wls4.2,nsim = 999);
plot(bin2,main="CONFIDENCE BANDS FOR EMPIRICAL VARIOGRAM\nMODEL:Matern", lwd=2,pch=16, envelope = env);


##################################
#Estimate parameters for trend and variogram

bin2 <- variog(geoscal.rot,option = "bin", pairs.min=30, max.dist=maxdist/2,
               estimator.type = "modulus",uvec=seq(0,maxdist/2,l=30),
               trend = ~geoscal.rot$coords[,1]+poly(geoscal.rot$coords[,2],3));


lk1.2 <- likfit(geoscal.rot ,cov.model = "exponential", ini =c(3,0.3),
              fix.nugget = F, nugget =0.4 ,lik.method = "REML"
              ,trend = ~geoscal.rot$coords[,1]+poly(geoscal.rot$coords[,2],3));

##########################################
#Estimate the anisotropy angle and ratio. (only using REML)

lk1.2 <- likfit(geoscalresidus.rot ,cov.model = "exponential", ini =c(3,0.3),
                fix.nugget = F, nugget =1.41 ,lik.method = "REML",
                #trend = ~geoscal.rot$coords[,1]+poly(geoscal.rot$coords[,2],3),
                fix.psiA = TRUE, psiA = pi/2, fix.psiR = TRUE, psiR = 3);

summary(lk1.2)
lk1.2
plot(bin2,  main="PARAMETRIC VARIOGRAMS",cex.main=1, pch=16);
lines(lk1.2, lwd = 2, col="red", max.dist=maxdist/2)


