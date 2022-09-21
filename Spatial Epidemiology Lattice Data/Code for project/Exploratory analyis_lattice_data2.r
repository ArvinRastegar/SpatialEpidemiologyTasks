

setwd("C:/Users/rabellana/OneDrive - Universitat de Barcelona/Docencia/Master d'Estadistica i IO/Lattice data/Sessió 1 i 2")

#Load packages

install.packages("maptools")
install.packages("sf")
install.packages("spdep")

library(sf)
library(spdep)
library(maptools)


###North Carolina SIDS data
####################################################################.
#DATA: 100 COUNTIES OF NORTH CAROLINA INCLUDES COUNTS OF NUMBERS OF LIVER BIRTHS
#(ALSO NON-WHITE LIVE BIRTHS) AND NUMBERS OF SUDDEN INFANT DEATHS, FOR 1974-1978 AND 1979-1984.
#THE COUNTY SEAT LOCATION COORDINATES ARE GIVEN IN MILES.


#The nc.sids data frame has 100 rows and 21 columns. It contains data given in Cressie (1991,
#pp. 386-9), Cressie and Read (1985) and Cressie and Chan (1989) on sudden infant deaths in North
#Carolina for 1974-78 and 1979-84. The data set also contains the neighbour list given by Cressie
#and Chan (1989) omitting self-neighbours (ncCC89.nb), and the neighbour list given by Cressie
#and Read (1985) for contiguities (ncCR85.nb). The data are ordered by county ID number, not
#alphabetically as in the source tables sidspolys is a "polylist" object of polygon boundaries, and
#sidscents is a matrix of their centroids.
#Usage

#Format
#This data frame contains the following columns:
#SP_ID SpatialPolygons ID
#CNTY_ID county ID
#east eastings, county seat, miles, local projection
#north northings, county seat, miles, local projection
#L_id Cressie and Read (1985) L index
#M_id Cressie and Read (1985) M index
#names County names
#AREA County polygon areas in degree units
#PERIMETER County polygon perimeters in degree units
#CNTY_ Internal county ID
#NAME County names
#FIPS County ID
#FIPSNO County ID
#CRESS_ID Cressie papers ID
#BIR74 births, 1974-78
#SID74 SID deaths, 1974-78
#NWBIR74 non-white births, 1974-78
#BIR79 births, 1979-84
#SID79 SID deaths, 1979-84
#NWBIR79 non-white births, 1979-84

#*******************************************.


load("data_nc_sids.Rdata")

nc.sids
slot(nc.sids,"data")

#Compute sudden infant deaths rates between 1974-1978
rate<-(1000)*nc.sids$SID74/nc.sids$BIR74
hist(rate)
#Freeman-Tukey Transfomation 
ft.SID74 <- sqrt(1000)*(sqrt(nc.sids$SID74/nc.sids$BIR74) + sqrt((nc.sids$SID74+1)/nc.sids$BIR74))

#Add rate and transfomation rate  to SpatialPolygonsDataFrame

nc.sids$ft.SID74<-ft.SID74   
nc.sids$rate<-rate
slot(nc.sids,"data")

#Creates the adjacency matrix

# Contiguity based neighbors

#Areas sharing any boundary point (QUEEN) are taken as neighbors

#Poly2nb:Construct neighbours list from polygon list based on contiguous boundaries

xxnb <- poly2nb(nc.sids)
class(xxnb)
plot(nc.sids)#, border="grey")   
plot(xxnb, coordinates(nc.sids), add=TRUE, col="blue")


#ASSIGN WEIGHTS TO THE AREAS THAT ARE LINKED  

#nb2listw:The function supplements a neighbours list with spatial weights 
#for the chosen coding scheme.

#a.#B:Binary (1: neighbour, 0:otherwise)
#b.#W: standardised to sum unity row.

b.sids<-nb2listw(xxnb,  style="B",  zero.policy=TRUE)   #Spatial weights for neighbours lists
w.sids<-nb2listw(xxnb,  style="W",  zero.policy=TRUE)   #Spatial weights for neighbours lists


#Moran's I test for spatial autocorrelation

# 
sids.moran.B<-moran.test(nc.sids$ft.SID74,b.sids)
sids.moran.B

sids.moran.W<-moran.test(nc.sids$ft.SID74 ,w.sids)
sids.moran.W 


#Permutation test for Moran's I statistic

nsim <- 1000
set.seed(1234)

sids.moran.mc.B<-moran.mc(nc.sids$ft.SID74, listw=b.sids, nsim=nsim)
sids.moran.mc.B

sids.moran.mc.B$res
hist(sids.moran.mc.B$res)
max(sids.moran.mc.B$res[1:nsim])

sids.moran.mc.W<-moran.mc(nc.sids$ft.SID74, listw=w.sids, nsim=nsim)
sids.moran.mc.W

sids.moran.mc$res
mean(sids.moran.mc$res[1:nsim])
var(sids.moran.mc$res[1:nsim])
hist(sids.moran.mc$res)
max(sids.moran.mc$res[1:nsim])



#geary.test Geary's C test for spatial autocorrelation

sids.geary<-geary.test(nc.sids$ft.SID74, w.sids)
sids.geary

sids.geary.mc<-geary.mc(nc.sids$rate, w.sids, nsim)
sids.geary.mc
sids.geary.mc$res

mean(sids.geary.mc$res[1:nsim])
var(sids.geary.mc$res[1:nsim])
hist(sids.geary.mc$res)
min(sids.geary.mc$res[1:nsim])



#Scatterplot moran Tests

moran.plot(ft.SID74,w.sids ,labels=as.factor(nc.sids$NAME), pch=19)


#Plots to show the results

label.colors<-c(rep("black",100))
label.colors[nc.sids$NAME=="Northampton"]<-"red"
plot(nc.sids)
text(coordinates(nc.sids),label=as.factor(nc.sids$NAME),cex=0.5,col=label.colors)
plot(nc.sids)
text(coordinates(nc.sids),label=round(ft.SID74,1),cex=0.7,col=label.colors)

moran.plot(ft.SID74,w.sids ,labels=as.factor(nc.sids$NAME), pch=19)

label.colors<-c(rep("black",100))
label.colors[nc.sids$NAME=="Swain"]<-"red"
plot(nc.sids)
text(coordinates(nc.sids),label=as.factor(nc.sids$NAME),cex=0.5,col=label.colors)
plot(nc.sids)
text(coordinates(nc.sids),label=round(ft.SID74,1),cex=0.7,col=label.colors)

moran.plot(ft.SID74,w.sids ,labels=as.factor(nc.sids$NAME), pch=19)

label.colors<-c(rep("black",100))
label.colors[nc.sids$NAME=="Richmond"]<-"red"
plot(nc.sids)
text(coordinates(nc.sids),label=as.factor(nc.sids$NAME),cex=0.5,col=label.colors)
plot(nc.sids)
text(coordinates(nc.sids),label=round(ft.SID74,1),cex=0.7,col=label.colors)




#Local indicators.

sids.local.moran<-localmoran(nc.sids$ft.SID74, w.sids)
#Northampton is the region 5
sids.local.moran[5,]
#Swain is the region 58
sids.local.moran[58,]
#Richmond is the region number 89
sids.local.moran[89,]

head(sids.local.moran)

#Only to look the results
data.frame(nc.sids$NAME,sids.local.moran)  

#Add I Moran results to spatial polygon dataframe

slot(nc.sids,"data")<-cbind(slot(nc.sids,"data"),sids.local.moran)

#Plot  I local Moran
#Intervals

l.inf<-round(min(nc.sids$Ii),digits=2)
l.sup<-round(max(nc.sids$Ii),digits=2)
library(classInt)


I.q4<-findInterval(nc.sids$Ii,c((l.inf-0.1),-0.5,0,0.5,(l.sup+0.1)))
breakI<-c((l.inf-0.1),-0.5,0,0.5,(l.sup+0.1))

#Colours
#pal.I<-gray.colors(5,0.95,0.55)
pal.I<- c("yellow", "orange2", "red3", "brown")

#Color assignment based on intervals
plot(nc.sids, col=pal.I[I.q4])
legend("topright", legend=leglabs(breakI),fill=pal.I, bty="n")
text(coordinates(nc.sids),label=as.factor(nc.sids$NAME),cex=0.5,col=label.colors)

#Plot  p value
Ip.2<-findInterval(nc.sids$"Pr(z != E(Ii))",c(0,0.05, 1))
breaksp<-c(0,0.05, 1)
pal.p<- gray.colors(2,0.4,0.95)

plot(nc.sids, col=pal.p[Ip.2])
legend("topright", legend=leglabs(breaksp),fill=pal.p, bty="n")
text(coordinates(nc.sids),label=as.factor(nc.sids$NAME),cex=0.5,col=label.colors)


#Join  the plots

par(mfrow=c(2,1),mai=c(0.5, 0.5, 0.5, 0))
plot(nc.sids, col=pal.I[I.q4])
legend("bottomright", legend=leglabs(breakI),fill=pal.I, bty="n")
plot(nc.sids, col=pal.p[Ip.2])
legend("bottomright", legend=leglabs(breaksp),fill=pal.p, bty="n")



#plot of the Rates

par(mfrow=c(1,1))
windows()
brks <- round(quantile(nc.sids$ft.SID74, probs=seq(0,1,0.2)), digits=2)   
colours <- c("yellow", "orange2", "red3", "brown", "black")

plot(nc.sids, col=colours[findInterval(nc.sids$ft.SID74, brks,all.inside=TRUE)])
legend(x=c(-84, -80), y=c(33, 34.5), legend=leglabs(brks),fill=colours, bty="n")
text(coordinates(nc.sids),label=as.factor(nc.sids$names),cex=0.55)
invisible(title(main=paste("Rate (Transformed FT) in North Carolina","For the 1974-1978 period")))