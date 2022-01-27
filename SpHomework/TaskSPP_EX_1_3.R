#setwd("C:/Users/johan/Google Drive/UNI/UPC/Spatial Epidemiology/Task Spatial Point Patterns")
setwd("~/Desktop/SpHomework/TASK 3 files-20220101")
##### IMPORT #####
library(spatstat)
#install.packages("RandomFields")
#install.packages("~/Downloads/RandomFieldsUtils_1.0.11.tgz", repos = NULL, type='source')
library(RandomFieldsUtils)
library(RandomFields)
#data
win23=read.table("pol23.txt", header = TRUE) #outline long 23
nucli23=read.table("nucli23.txt", header = TRUE) #coordinates of nests in 23
#win84=read.table("poly84.txt", header = TRUE) #outline round 84
#nucli84=read.table("nucli84.txt", header = TRUE) #coordinates of nests & vars in  84

#### EXERCISE 1 #####

# EX1.1: Build a ppp object using the "win23" data as a window
p23=ppp(win23$X,win23$Y, c(-10,120),c(-10,120))
pol.illa<-list(x=win23$X,y=win23$Y)
p23p=ppp(nucli23$X,nucli23$Y,poly=pol.illa)
plot(p23p,main="Pol23 Spatial Data")
axis(1,at=c(seq(0,100,by=10)))
axis(2,at=c(seq(0,120,by=10)),pos=c(-10,0))

## EX 1.2: Briefly describe the point pattern process
# point process descriptive analysis
p23p.den<-density(p23p, dimyx=c(512,512), kernel = "gaussian")
p23p.den
summary(p23p.den)
plot(p23p.den, main="Pol23 Data")
axis(1,at=c(seq(0,80,by=10)))
axis(2,at=c(seq(0,110,by=10)),pos=c(-10,0))

# Using various variances
par(mfrow=c(1,2))
plot(density(p23p, dimyx=c(256,256), sigma=50), main="Density Plot \n (Gaussian Kernel with sigma=50)")
plot(density(p23p, dimyx=c(256,256), sigma=10), main="Density Plot \n (Gaussian Kernel with sigma=10)")
plot(density(p23p, dimyx=c(256,256), sigma=5), main="Density Plot \n (Gaussian Kernel with sigma=5)")
plot(density(p23p, dimyx=c(256,256), sigma=1), main="Density Plot \n (Gaussian Kernel with sigma=1)")

## EX 1.3: Rebuild the ppp object using the "time to nest" (DATAPOS variable) as marks
nucli84=read.table("nucli84.txt", header = TRUE)
win84=read.table("poly84.txt", header = TRUE)
summary(nucli84)

nucli84$X=nucli84$X-min(win84$X)
nucli84$Y=nucli84$Y-min(win84$Y)
win84$X=win84$X-min(win84$X)
win84$Y=win84$Y-min(win84$Y)


nuc.illa<-list(x=win84$X,y=win84$Y)
p84=ppp(win84$X,win84$Y, c(-10,150),c(-10,180))
p84p=ppp(nucli84$X,nucli84$Y,poly=nuc.illa)
plot(p84p,main="Pol84 Spatial Data")
axis(1,at=c(seq(0,150,by=10)))
axis(2,at=c(seq(0,180,by=10)),pos=c(-10,0))


p84T=ppp(nucli84$X,nucli84$Y,poly=nuc.illa,marks=nucli84$DATAPOS )
summary(p84T)

plot(p84T,main="Nucli84 Spatial Data", maxsize = 10,leg.side="right")
axis(1,at=c(seq(0,150,by=25)))
axis(2,at=c(seq(-10,180,by=25)),pos=c(-10,0))

## EX1.4: Briefly describe the marked point process
p84p.den<-density(p84p, dimyx=c(512,512), kernel = "gaussian")
plot(p84p.den, main="Pol23 Data")
axis(1,at=c(seq(0,150,by=10)))
axis(2,at=c(seq(0,180,by=10)),pos=c(-10,0))

# Using various variances
par(mfrow=c(1,2))
plot(density(p84p, dimyx=c(256,256), sigma=50), main="Density Plot \n (Gaussian Kernel with sigma=50)")
plot(density(p84p, dimyx=c(256,256), sigma=10), main="Density Plot \n (Gaussian Kernel with sigma=10)")
plot(density(p84p, dimyx=c(256,256), sigma=5), main="Density Plot \n (Gaussian Kernel with sigma=5)")
plot(density(p84p, dimyx=c(256,256), sigma=1), main="Density Plot \n (Gaussian Kernel with sigma=1)")

#### EXERCISE 2 ##### 

## EX2.1: Give a point estimate and a 95% confidence interval of the intensity assuming a homogeneous Poisson process

#point estimate
int=summary(p23p)$intensity
log(int)
#CI for log(int)
model=ppm(p23p, ~1)
model
exp(-3.475206) #0.007640754
exp(-3.070896) #0.011222

## EX2.2: Assess the Completely Spatial Randomness hypothesis

Q <- quadratcount(p23p, nx=1, ny=4)
plot(p23p, cex = 0.5, main="Quadrat Count Poly23")
plot(Q, add = TRUE)
axis(1,at=c(seq(0,80,by=10)))
axis(2,at=c(seq(0,110,by=10)),pos=c(-10,0))
Q


M <- quadrat.test(p23p, nx=1, ny=4)
M
plot(p23p)
plot(M, add = TRUE)
axis(1,at=c(seq(0,100,by=25)))
axis(2,at=c(seq(0,180,by=25)),pos=c(-10,0))

# In case expected<5- P-value Montecarlo simulation
obs<-M$observed
# Compute the area of each subregion
Z<-quadratcount(p23p, nx=1, ny=4)
quadratsM <- tiles(as.tess(Z))
quadrat.areas <- unlist(lapply(quadratsM, area.owin))
# Expected proportions
p.exp<-quadrat.areas/(sum(quadrat.areas))
chisq.test(M$observed,p=p.exp)$expected
chisq.test(M$observed,p=p.exp,simulate.p.value = T)


## EX2.3: Fit an inhomogeneous Poisson model - Give the point estimates and 95% confidence interval of parameters.

model=ppm(p23p, ~x+y)
summary(model)

#x: 0.06156428  
#y: 0.01003627 
model
par(mfrow=c(1,2))
plot(model,  dimyx=c(256,256), axes=TRUE)

## EX2.4: Fit an inhomogeneous Poisson model - Check the goodness-of-fit of the model
AIC(model)
par(mfrow=c(1,2))
KS=cdf.test(model,"y")
KS
plot(KS, main="KS test of inhomogeneous Poisson Process \n (based on distribution of y coordinate)")

KS=cdf.test(model,"x")
KS
plot(KS, main="KS test of inhomogeneous Poisson Process \n (based on distribution of x coordinate)")

## EX2.5: Fit an inhomogeneous Poisson model - Make a residual validation of the model

plot(predict(model), main="Predicated Model")
plot(p23p, add = TRUE, pch = "+")

sum(eem(model))/area(p23p)

#Smoothed residuals
xx<-diagnose.ppm(model, which = "smooth",type="pearson",sigma=5)
axis(1,at=c(seq(0,100,by=25)), pos = c(-10,0))
axis(2,at=c(seq(0,180,by=25)),pos=c(-10,0))

# Null standard deviation = 1/(2*bandwidth*sqrt(pi))
3/(2*50*sqrt(pi))

lurking(model, expression(x), type = "raw", main="Lurking variable plot \n wrt. x coordinate")
lurking(model, expression(y), type = "raw",main="Lurking variable plot \n wrt. y coordinate" )


#### EXERCISE 3 #####
# Link for explanations.
# https://www.geo.fu-berlin.de/en/v/soga/Geodata-analysis/Spatial-Point-patterns/Analysis-of-spatial-point-patterns/Interactions/index.html
## EX3.1: Explore the pattern of interaction.
# Empty space function
par(mfrow = c(1,1))
?p23p
plot(p23p)
axis(1,at=c(seq(0,100,by=10)))
axis(2,at=c(seq(0,120,by=10)),pos=c(-10,0))
class(p23p)

# Empty space distance 
par(mfrow = c(1,3))
Fp <- Fest(p23p)
Fp
par(pty = "s")
plot(Fest(p23p,correction="best"),legend=F, main = "Fest")

# NN function (distances)
Gp <- Gest(p23p)
Gp
par(pty = "s")
plot(Gest(p23p,correction="best"),legend=F, main = "Gest")

# Pairwise distances
Kp <- Kest(p23p)
Kp
par(pty = "s")
plot(Kest(p23p,correction="best"),legend=F, main = "Kest")
 
# Envelope
par(mfrow = c(1,3))
Ep <- envelope(p23p, Kest, nsim = 39,correction="best")
Ep
plot(Ep, main = "pointwise envelopes",legend=F)

Ep <- envelope(p23p, Kest, nsim = 19, rank = 1, global = TRUE,correction="best")
Ep
plot(Ep, main = "global envelopes",legend=F)

Ep <- envelope(p23p, Lest, nsim = 19, rank = 1, global = TRUE,correction="best")
plot(Ep, main = "global envelopes of L(r)", legend=F)

#### EX3.2: Fit a log-Gaussian Cox process to data. Comment the results.
# Log-Gaussian Cox process
# Log Gaussian Cox
library(RandomFields)
par(mfrow = c(1,1))
grid<-read.table("grid_height_23.txt",header=T)
grid_veg<-read.table("grid_veg_poly23.txt",header=T)

mat<-as.matrix(read.table("height_23.txt"))
height<-im(mat,grid$x,grid$y)
plot(height,axis=T)
plot(p23p, add=T, cex=0.5)


mat<-as.matrix(read.table("veg_poly23.txt"))
veg<-im(mat,grid_veg$x,grid_veg$y)
plot(veg,axis=T)
plot(p23p, add=T, cex=0.5)
fitLG <- kppm(p23p, ~height+veg, "LGCP")
summary(fitLG)
plot(fitLG,what="statistic", correction='best')
plot(fitLG, main = "log Gaussian Cox process")

sLG<-fitLG$par[1]
aLG<-fitLG$par[2]

# Covariance at distance 10
sLG*exp(-(6/aLG))

#Autocorrelation function
r<-seq(0,20,1)
plot(r,exp(-(r/aLG)),type="l", main = "Autocorrelation function")

# Range with autocorrelation = 0.2
aLG*log(0.2)*(-1)

# simulation

# Trend function

trend=as.im(log(predict(fitLG)),W=p23p$window)
#Xsim=rLGCP("exp", mu =trend, var=fitLG$modelpar[1],scale=fitLG$modelpar[2],win=p23p$window)
#plot(Xsim)

# E <- envelope(fitLG, Lest, nsim = 19, rank = 1, global = TRUE,correction="none",
#               simulate=expression(Xsim=rLGCP("exp", mu = trend, var=fitLG$modelpar[1],
#                                              scale=fitLG$modelpar[2],win=p23p$window)
#               ))

plot(E, main = "global envelopes of L(r)", legend=F)



#### EX3.3: Fit a Gibbs process model to data. Comment the results ####
# Gibbs process

# Area Interaction
# It takes some time 80s
df=data.frame(r=seq(3,6,by=0.5))
date()
pfit=profilepl(df,AreaInter,p23p,~height+veg,
               correction="none")
date()
summary(pfit)
plot(pfit)

fitA<-ppm(p23p, ~height+veg, AreaInter(r=3.5))
summary(fitA)
par(mfrow = c(1,2))
plot(fitA)

E <- envelope(fitA, Lest, nsim = 19, rank = 1, global = TRUE,correction="none")
par(mfrow = c(1,1))
plot(E, main = "global envelopes of L(r)")


# Geyer
df=data.frame(r=seq(4,15,by=0.5),sat=1:23)
pfit1=profilepl(df,Geyer,p23p,~height+veg,correction="none")
pfit1
plot(pfit1)

sat<-2:7
r<-5:10

library(tidyr)

df=as.data.frame(crossing(r,sat))
pfit2=profilepl(df,Geyer,p23p,~height+veg)
pfit2
plot(pfit2)


fitG=ppm(p23p, ~height+veg, Geyer(r = 6, sat = 7))
fitG
summary(fitG)
summary(fitA)

par(mfrow = c(1,2))
plot(fitG)



# E <- envelope(fitG, Lest, nsim = 19, rank = 1, global = TRUE,correction="none",
#               simulate=expression(Xsim=rLGCP("exp", mu = trend, var=fitLG$modelpar[1],scale=fitLG$modelpar[2],win=p23p$window)
#               ))
# plot(E, main = "global envelopes of L(r)")

AIC(fitA)
AIC(fitG)


fitG2=ppm(p23p, ~height, Geyer(r = 6, sat = 7))
AIC(fitG2)
summary(fitG2)
plot(fitG2)


#Smoothed residuals
#Smoothed residuals
par(mfrow = c(1,1))
diagnose.ppm(fitA, which = "smooth",type="pearson",rbord=0, legend=T)



diagnose.ppm(fitG, which = "smooth",type="pearson",rbord=0, legend=T)
diagnose.ppm(fitG2, which = "smooth",type="pearson",rbord=0, legend=T)

diagnose.ppm(fitG2,type="pearson")
axis(1,at=c(seq(0,100,by=10)))
axis(2,at=c(seq(0,120,by=10)),pos=c(-10,0))
# Null standard deviation = 1/(2*bandwidth*sqrt(pi))
0.0158*2
diagnose.ppm(fitA, which = "smooth",type="pearson",rbord=0, legend=T)

# Lurking variable plot
par(mfrow = c(1,3))
lurking(fitG2, expression(x), type = "raw",envelope=T)
lurking(fitG2, expression(y), type = "raw",envelope=T)
lurking(fitG2, height, type = "raw",envelope=T)


#### EXERCISE 4 #####

## EX4.1: Using the Primary Biliary Cirrhosis data (PBC) make a case-control point pattern analysis - Assessment of risk variability

## EX4.2: Using the Primary Biliary Cirrhosis data (PBC) make a case-control point pattern analysis - Assessment of interaction
