
# EXERCISE 4 #####

library(spatstat)
library(sparr)
pbc.data <- read.table('./PBCdata2.txt',header=T)
pbc.poly <- read.table('./PBCpoly.txt', header=T)

# EX4.1: Using the Primary Biliary Cirrhosis data (PBC) make a case-control point pattern analysis ----


## Assessment of risk variability ----

m <- factor(pbc.data$marks)
pbc <- ppp(pbc.data$x, pbc.data$y,
           poly=list(x=pbc.poly$x, y=pbc.poly$y),
           marks=m)
plot(pbc$window)
#points(split(pbc)$case, pch = 16, col = 1)
#points(split(pbc)$control, pch = 16, col = 2)
### plots ----

# split
cases <- split(pbc)$case
controls <- split(pbc)$control

# intensities
chp <- risk(cases, controls, log=F, intensity=T)

# Cases
plot(chp$f, ribbon=F)
# Controls
plot(chp$g)
# Relative Risk
plot(chp$rr)
# reference value
cases$n/controls$n

# Log Intensities
chp.log <- risk(cases,controls, intensity=T)
# Relative risk
plot(chp.log$rr)
# Reference value
log(cases$n/controls$n)
# Sealevel --> yellow
max(chp.log$rr)
min(chp.log$rr)
plot(chp.log$rr,col=beachcolours(c(-8,-1), sealevel = log(cases$n/controls$n)))


### Test ----

rho0<-cases$n/controls$n  #0 # For Log-densities
chp <- risk(cases,controls)

cellsize<-chp$rr$xstep*chp$rr$ystep
ratiorho <- cellsize*sum((chp$rr$v-rho0)^2,na.rm=T)

# Permutation function

perm_rr<-function(){
  new_ch<-rlabel(pbc)
  new_cases <- split(new_ch)$case
  new_controls <- split(new_ch)$control
  new_chp <- risk(new_cases,new_controls)
  cellsize<-new_chp$rr$xstep*new_chp$rr$ystep
  ratio_perm <- cellsize*sum((new_chp$rr$v-rho0)^2,na.rm=T)
  ratio_perm
}

nsim<-500
set.seed(2022)
rperm<-sapply(1:nsim,function(i) perm_rr())

# P-value
(sum(rperm > ratiorho)+1)/(nsim+1)

plot(density(rperm))
abline(v=ratiorho)


# P-values Contour Plot
chp.contour <- risk(cases,controls,tolerate = T,adapt=T)
plot(chp.contour)
plot(chp.contour$P)
chp_p<-cut(chp.contour$P,breaks=c(0,0.05,0.1,0.9,0.95,1))
V <- tess(image = chp_p)
plot(V,valuesAreColours=FALSE)



## Assessment of interaction ----

# Assessment of general spatial clustering


s=seq(0,2100,by=100)
khcases<-Kest(cases,r=s,correction="iso")
khcontrols<-Kest(controls,r=s, correction="iso")

plot(khcases)
lines(khcontrols$r,khcontrols$iso,lty=2)

plot(khcontrols)

# Difference of k functions with envelope
Xppp=pbc
Kdif<-function(Xppp,r,cr="iso")
{
  k1<-Kest(Xppp[marks(Xppp)=="case"], r=r, correction=cr)
  k2<-Kest(Xppp[marks(Xppp)=="control"],r=r,  correction=cr)
  D=k1[[cr]]-k2[[cr]]
  res<-data.frame(r=r, D=D)
  return(fv(res, valu="D", fname="D"))
}

nsim<-39
envKdif<-envelope(pbc, Kdif, nsim=nsim, nrank=1,r=s,
                  savefuns=TRUE,simulate=expression(rlabel(pbc)))


plot(envKdif,legend=T)

# Test

# Extract the simulated functions
simfuns<-as.data.frame(attr(envKdif, "simfuns"))[,-1]
#Compute diagonal of var-cov matrix of dif. K-functions
khcovdiag<-apply(simfuns, 1, var)


#Test
T0<-sum( ((khcases$iso-khcontrols$iso)/sqrt(khcovdiag))[-1])
T_pm<-apply(simfuns, 2, function(X){
  sum((X/sqrt(khcovdiag))[-1])
})

plot(density(T_pm), xlim=c(-100,260))
abline(v=T0, col='red')

pvalue<-2*(sum(abs(T_pm)>abs(T0))+1)/(nsim+1)
pvalue


# Alternative option for limits
plot(envKdif, legend=F)
lines(s, -1.96*sqrt(khcovdiag), lty=2)
lines(s, +1.96*sqrt(khcovdiag), lty=2)
