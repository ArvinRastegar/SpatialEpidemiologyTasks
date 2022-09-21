library(spatstat)

# Empty space function
data(cells)
?cells
plot(cells)
class(cells)





########################
# Inhomogeneous Data ##
#######################

nucli.84=read.table(file="nucli84.txt",header=T,sep="\t")
poly.84=read.table(file="poly84.txt",header=T,sep="\t")

#Create ppp object
summary(nucli.84)

min.X=min(nucli.84$X)
min.Y=min(nucli.84$Y)

nucli.84$X=nucli.84$X-min.X
nucli.84$Y=nucli.84$Y-min.Y
summary(nucli.84)

#Create ppp object
# Poligon
poligon=poly.84
summary(poly.84)
poligon$X=poligon$X-min.X
poligon$Y=poligon$Y-min.Y
summary(poligon)


pol.illa<-list(x=poligon$X,y=poligon$Y)
plot(owin(poly=pol.illa),main="nucli84")
n84p=ppp(nucli.84$X,nucli.84$Y,poly=pol.illa)
plot(n84p,main="nucli84")
axis(1,at=c(seq(-25,100,by=25)),pos=c(-75,0))
axis(2,at=c(seq(-75,100,by=25)),pos=c(-25,0))



# Covariables: height and vegetation

# Load data
grid<-read.table("grid.txt",header=T)
grid_veg<-read.table("grid_veg.txt",header=T)

mat<-as.matrix(read.table("height.txt"))
height<-im(mat,grid$x,grid$y)
plot(height,axis=T)
plot(n84p, add=T, cex=0.5)


mat<-as.matrix(read.table("veg.txt"))
veg<-im(mat,grid_veg$x,grid_veg$y)
plot(veg,axis=T)
plot(n84p, add=T, cex=0.5)


# Inhomogeneous PP
fit <- ppm(n84p, ~height + veg)
lam <- predict(fit, locations = n84p)
Ki <- Kinhom(n84p, lam,correction="none")
plot(Ki, main = "Inhomogeneous K function",legend=F)

# Envelope
E <- envelope(fit, Lest, nsim = 39, rank = 1,correction="none")
plot(E, main = "pointwise envelopes")

E <- envelope(fit, Lest, nsim = 19, rank = 1, global = TRUE,correction="none")
plot(E, main = "global envelopes of L(r)")


# Models

#Thomas
fit <- kppm(n84p, ~height + veg, "Thomas")
fit

plot(fit, what="statistic")
plot(fit)

plot(predict(fit))
plot(n84p, add = TRUE, pch = "+")


# Envelope
E <- envelope(fit, Lest, nsim = 39, rank = 1,correction="none")
plot(E, main = "pointwise envelopes")

E <- envelope(fit, Lest, nsim = 19, rank = 1, global = TRUE,correction="none")
plot(E, main = "global envelopes of L(r)")


#Simulate data
Xsim=rThomas(fit$par[1], sqrt(fit$par[2]), fit$mu, win=n84p$window)
plot(Xsim)


# Matern
fitM <- kppm(n84p, ~height + veg, "MatClust")
fitM
plot(fitM,what="statistic")


#Simulate data
Xsim=rMatClust(fitM$par[1],fitM$par[2],fitM$mu,win=n84p$window)
plot(Xsim)

E <- envelope(fitM, Lest, nsim = 19, rank = 1, global = TRUE,correction="none")
plot(E, main = "global envelopes of L(r)")


plot(predict(fitM))
plot(n84p, add = TRUE, pch = "+")


# Log Gaussian Cox
fitLG <- kppm(n84p, ~height+veg, "LGCP")
fitLG

plot(fitLG,what="statistic")
plot(fitLG)

sLG<-fitLG$par[1]
aLG<-fitLG$par[2]

# Covariance at distance 10
sLG*exp(-(6/aLG))

#Autocorrelation function
r<-seq(0,20,1)
plot(r,exp(-(r/aLG)),type="l")

# Range with autocorrelation = 0.2
aLG*log(0.2)*(-1)

# simulation

# Trend function

trend=as.im(log(predict(fitLG)),W=n84p$window)
# Xsim=rLGCP("exp", mu =trend, var=fitLG$modelpar[1],scale=fitLG$modelpar[2],win=n84p$window)
plot(Xsim)

# E <- envelope(fitLG, Lest, nsim = 19, rank = 1, global = TRUE,correction="none",
#               simulate=expression(Xsim=rLGCP("exp", mu = trend, var=fitLG$modelpar[1],scale=fitLG$modelpar[2],win=n84p$window)
# ))

plot(E, main = "global envelopes of L(r)", legend=F)



#### Gibbs process ####

# Area Interaction
# It takes some time 80s
# df=data.frame(r=seq(3,6,by=0.5))
# date()
# pfit=profilepl(df,AreaInter,n84p,~height+veg,
#                correction="none")
# date()
# pfit
# plot(pfit)

# fitA<-ppm(n84p, ~height+veg, AreaInter(r=3.5))
summary(fitA)
plot(fitA)

E <- envelope(fitA, Lest, nsim = 19, rank = 1, global = TRUE,correction="none")
plot(E, main = "global envelopes of L(r)")


# Geyer
df=data.frame(r=seq(4,15,by=0.5),sat=1:23)
pfit=profilepl(df,Geyer,n84p,~height+veg,correction="none")
pfit
plot(pfit)

sat<-2:7
r<-5:10

library(tidyr)

df=as.data.frame(crossing(r,sat))
pfit=profilepl(df,Geyer,n84p,~height+veg,correction="none")
pfit
plot(pfit)


fitG=ppm(n84p, ~height+veg, Geyer(r = 6, sat = 7))
fitG
summary(fitG)
plot(fitG)



# E <- envelope(fitG, Lest, nsim = 19, rank = 1, global = TRUE,correction="none",
#               simulate=expression(Xsim=rLGCP("exp", mu = trend, var=fitLG$modelpar[1],scale=fitLG$modelpar[2],win=n84p$window)
#               ))
# plot(E, main = "global envelopes of L(r)")


AIC(fitA)
AIC(fitG)


fitG2=ppm(n84p, ~height, Geyer(r = 6, sat = 7))
AIC(fitG2)
summary(fitG2)
plot(fitG2)


#Smoothed residuals
#Smoothed residuals
diagnose.ppm(fitG2, which = "smooth",type="pearson",rbord=0)

# Null standard deviation = 1/(2*bandwidth*sqrt(pi))
0.0158*2


# Lurking variable plot
lurking(fitG2, expression(x), type = "raw",envelope=T)
lurking(fitG2, expression(y), type = "raw",envelope=T)
lurking(fitG2, height, type = "raw",envelope=T)

