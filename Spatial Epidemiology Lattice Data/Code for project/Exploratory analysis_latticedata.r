#Load Workspace nc_sids data
#load("data_nc_sids.Rdata")

#Load packages

install.packages("maptools")
install.packages("sf")
install.packages("spdep")

library(sf)
library(spdep)
library(maptools)



#legend:colour palletes rainbow, grey.colours,heat.colours,terrain.colours topo.colours,cm.colours( cyan-magenta default)
  library(class)
  install.packages("RColorBrewer")
  library(RColorBrewer)     #Creates nice looking color palettes especially for thematic maps
#class intervals.
  library(classInt)


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
#Import data form the Noth Carolina 
#Three files different extensions: 
#.shp (geometric data) 
#.shx  (spatial indexs)  
#.dbf (spatial data information: variables )  
  
  
#st_read::Read simple features or layers from file or database

nc <-sf::st_read(system.file("shape/nc.shp", package="sf"))
class(nc)
#convert a sf to spatial*Dataframe
spd <- sf::as_Spatial(nc)

#interesting web page: https://cran.r-project.org/web/packages/sf/vignettes/sf5.html 

class(spd)
# [1] "SpatialPolygons"
# attr(,"package")
# [1] "sp"

## grab the data from the sf object
df <- nc
df$geometry <- NULL
df <- as.data.frame(df)

## create the SpatialPolygonsDataFrame
nc.sids <- sp::SpatialPolygonsDataFrame(spd, data = df )

######## 

summary(nc.sids)
slot(nc.sids,"data")


#Compute sudden infant deaths rates between 1974-1978
rate<-(1000)*nc.sids$SID74/nc.sids$BIR74
hist(rate)
#Freeman-Tukey Transfomation
ft.SID74 <- sqrt(1000)*(sqrt(nc.sids$SID74/nc.sids$BIR74) + sqrt((nc.sids$SID74+1)/nc.sids$BIR74))
summary(ft.SID74)
shapiro.test(ft.SID74)
hist(ft.SID74)

#Add rate and transfomation rate  to SpatialPolygonsDataFrame

  nc.sids$ft.SID74<-ft.SID74   
  nc.sids$rate<-rate
  slot(nc.sids,"data")

  plot(nc.sids, border="blue", axes=TRUE, las=1)
  text(coordinates(nc.sids),label=nc.sids$NAME,cex=0.5)


#Plot of the ft.SID74

  brks <- round(quantile(nc.sids$ft.SID74, probs=seq(0,1,0.2)), digits=2)   
  colours <- c("yellow", "orange2", "red3", "brown", "black")
  
  plot(nc.sids, col=colours[findInterval(nc.sids$ft.SID74, brks,all.inside=TRUE)])
  legend(x=c(-84, -80), y=c(33, 34.5), legend=leglabs(brks),fill=colours, bty="n",cex=0.5)
  title(main=paste("Rate (Transformed FT) in North Carolina","For the 1974-1978 period"))
  text(coordinates(nc.sids),label=as.factor(nc.sids$NAME),cex=0.5)

# Different color Palette

  colours<-gray.colors(5,0.95,0.2)
  colours<- terrain.colors(5)
  colours<- rev(terrain.colors(5))
  
  plot(nc.sids, col=colours[findInterval(nc.sids$ft.SID74, brks,all.inside=TRUE)])
  legend(x=c(-84, -80), y=c(33, 34.5), legend=leglabs(brks),fill=colours, bty="n",cex=0.5)
  invisible(title(main=paste("Rate (Transformed FT) in North Carolina","For the 1974-1978 period")))
  text(coordinates(nc.sids),label=as.factor(nc.sids$NAME),cex=0.5)

# CHOOSE A NEIGHBORHOOD CRITERION USING R

#a. Areas sharing any boundary point (QUEEN) are taken as neighbors

#Poly2nb:Construct neighbours list from polygon list based on contiguous boundaries

  xxnb <- poly2nb(nc.sids)
  class(xxnb)
  plot(nc.sids)#, border="grey")   
  plot(xxnb, coordinates(nc.sids), add=TRUE, col="blue")

  summary.nb(xxnb)

  nc.sids$NAME[1]
  xxnb[[1]]
  nc.sids$NAME[xxnb[[1]]]
  
  #Show in the map the regions
  plot(nc.sids)#, border="grey")   
  label.colors<-c("red",rep("black",99))
  label.colors[xxnb[[1]]]<-"green"
  text(coordinates(nc.sids),label=as.factor(nc.sids$NAME),cex=0.5,col=label.colors)


#Region with largest number of contiguities

  cards <- card(xxnb)  #Cardinalities for neighbours lists
  sort(cards)
  which(cards == max(cards))
  maxconts <- which(cards == max(cards))[1]
  nc.sids$NAME[maxconts]
  nc.sids$FIPSNO[maxconts]
  #if(length(maxconts) > 1) maxconts <- maxconts[1]
  fg <- rep("grey", length(cards))
  fg[maxconts] <- "red"
  fg[xxnb[[maxconts]]] <- "green"
  plot(nc.sids, col=fg)
  text(coordinates(nc.sids), label=nc.sids$NAME,cex=0.5)
  title(main="Region with largest number of contiguities")

  maxconts <- which(cards == max(cards))[2]
  nc.sids$NAME[maxconts]
  #if(length(maxconts) > 1) maxconts <- maxconts[1]
  fg <- rep("grey", length(cards))
  fg[maxconts] <- "red"
  fg[xxnb[[maxconts]]] <- "green"
  plot(nc.sids, col=fg)
  text(coordinates(nc.sids), label=nc.sids$NAME,cex=0.5)
  title(main="Region with largest number of contiguities")

#b.Distance based neighbors k nearest neighbors

  coords<-coordinates(nc.sids)
  IDs<-row.names(as(nc.sids, "data.frame"))
  
  #knearneigh:K nearest neighbours for spatial weights
  #knn2nb:Neighbours list from knn object
  
  #k=1
  sids_kn1<-knn2nb(knearneigh(coords, k=1), row.names=IDs)
  class(sids_kn1)
  plot(nc.sids)
  plot(sids_kn1, coords, add=T)
  #k=2
  sids_kn2<-knn2nb(knearneigh(coords, k=2), row.names=IDs)
  plot(nc.sids)
  plot(sids_kn2, coords, add=T)
  #k=4
  sids_kn4<-knn2nb(knearneigh(coords, k=4), row.names=IDs)
  plot(nc.sids)
  plot(sids_kn4, coords, add=T)
  
#c. Distance based neighbors: Can also assign neighbors based on a specified distance
    #nbdists:spatial link distance measures
    #dnearneigh:neighbourhood contiguiy by distance
  
  dist<-unlist(nbdists(xxnb, coords))
  summary(dist)
  max_k1<-max(dist)
  
  sids_kd1<-dnearneigh(coords, d1=0, d2=0.75*max_k1, row.names=IDs)
  class(sids_kd1)
  plot(nc.sids)
  plot(sids_kd1, coords, add=T)
  
  sids_kd2<-dnearneigh(coords, d1=0, d2=1.5*max_k1, row.names=IDs)
  plot(nc.sids)
  plot(sids_kd2, coords, add=T)
  
#ASSIGN WEIGHTS TO THE AREAS THAT ARE LINKED  

#nb2listw:The function supplements a neighbours list with spatial weights 
#for the chosen coding scheme.

#a.#B:Binary (1: neighbour, 0:otherwise)
#b.#W: standardised to sum unity row.
  
  b.sids<-nb2listw(xxnb, glist=NULL, style="B",  zero.policy=TRUE)   #Spatial weights for neighbours lists
  w.sids<-nb2listw(xxnb, glist=NULL, style="W",  zero.policy=TRUE)   #Spatial weights for neighbours lists

  b.sids$weights
  summary(unlist(b.sids$weights))
  summary(sapply(b.sids$weights,sum))


  w.sids$weights
  summary(unlist(w.sids$weights))
  summary(sapply(w.sids$weights,sum))

#c. #General spatial weights.
  #glist:  argument (belive that the strenght of neighbour relationship attenuates with the distance.
  #weights to be proportional to the inverse distance between points representing the areas

  dsts<-nbdists(xxnb,coordinates(nc.sids))  #Distances between neighbors 
  dsts
  dsts[[1]]
  #idw<-lapply(dsts,function(x) 1/(x/10))    #Inversa de la distancia 
  idw<-lapply(dsts,function(x) 1/(x))    #Inversa de la distancia 
  dsts[[1]]
  idw[[1]]
  
  d.sids<-nb2listw(xxnb, glist=idw,style="B")
  attributes(d.sids)
  unlist(d.sids$weights)[1:6]
  unlist(w.sids$weights)[1:6]

#Moran's I test for spatial autocorrelation

  sids.moran<-moran.test(nc.sids$ft.SID74 ,w.sids)
  sids.moran 

  sids.moran2<-moran.test(nc.sids$ft.SID74 ,b.sids)
  sids.moran2 


  sids.moran3<-moran.test(nc.sids$rate ,w.sids)
  sids.moran3 


#Permutation test for Moran's I statistic

nsim <- 1000
set.seed(1234)

sids.moran.mc<-moran.mc(nc.sids$rate, listw=w.sids, nsim=nsim)
sids.moran.mc

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

  st.ft.SID74<-(nc.sids$ft.SID74-mean(nc.sids$ft.SID74))/sd(nc.sids$ft.SID74)
  xx<-moran.plot(st.ft.SID74,w.sids ,labels=as.factor(nc.sids$NAME), pch=19)
  xx<-moran.plot(nc.sids$ft.SID74,w.sids ,labels=as.factor(nc.sids$NAME), pch=19)

#Plots to show the results

  label.colors<-c(rep("black",100))
  label.colors[nc.sids$NAME=="Northampton"]<-"red"
  plot(nc.sids)
  text(coordinates(nc.sids),label=as.factor(nc.sids$NAME),cex=0.5,col=label.colors)
  plot(nc.sids)
  text(coordinates(nc.sids),label=round(st.ft.SID74,1),cex=0.7,col=label.colors)


  label.colors<-c(rep("black",100))
  label.colors[nc.sids$NAME=="Anson"]<-"red"
  plot(nc.sids)
  text(coordinates(nc.sids),label=as.factor(nc.sids$NAME),cex=0.5,col=label.colors)
  plot(nc.sids)
  text(coordinates(nc.sids),label=round(st.ft.SID74,1),cex=0.7,col=label.colors)


  label.colors<-c(rep("black",100))
  label.colors[nc.sids$NAME=="Richmond"]<-"red"
  plot(nc.sids)
  text(coordinates(nc.sids),label=as.factor(nc.sids$NAME),cex=0.5,col=label.colors)
  plot(nc.sids)
  text(coordinates(nc.sids),label=round(st.ft.SID74,1),cex=0.7,col=label.colors)




#Local indicators.

  sids.local.moran<-localmoran(nc.sids$ft.SID74, w.sids)
  #Northampton is the region 5
  #Anson is the region 85
  #Richmond is the region number 89
  
  head(sids.local.moran)
  #Only to look the results
  data.frame(nc.sids$NAME,sids.local.moran)  

#Add I Moran results to spatial polygon dataframe

  slot(nc.sids,"data")<-cbind(slot(nc.sids,"data"),sids.local.moran)

  #Plot  I local Moral
  #Intervals

  l.inf<-round(min(nc.sids$Ii),digits=2)
  l.sup<-round(max(nc.sids$Ii),digits=2)
  library(classInt)


  I.q4<-findInterval(nc.sids$Ii,c((l.inf-0.1),-0.5,0,1.6,(l.sup+0.1)))
  breakI<-c((l.inf-0.1),-0.5,0,1.6,(l.sup+0.1))

#Colours
#pal.I<-gray.colors(5,0.95,0.55)
  pal.I<- c("yellow", "orange2", "red3", "brown")

#Color assignment based on intervals
  plot(nc.sids, col=pal.I[I.q4])
  legend("topright", legend=leglabs(breakI),fill=pal.I, bty="n")

#Plot  p value
  Ip.2<-findInterval(nc.sids$"Pr(z > 0)",c(0,0.05, 1))
  breaksp<-c(0,0.05, 1)
  pal.p<- gray.colors(2,0.4,0.95)

  plot(nc.sids, col=pal.p[Ip.2])
  legend("topright", legend=leglabs(breaksp),fill=pal.p, bty="n")
  text(coordinates(nc.sids),label=as.factor(nc.sids$NAME),cex=0.5,col=label.colors)


#Els dos gràfics junts.

  par(mfrow=c(2,1),mai=c(0.5, 0.5, 0.5, 0))
  plot(nc.sids, col=pal.I[I.q4])
  legend("bottomright", legend=leglabs(breakI),fill=pal.I, bty="n")
  plot(nc.sids, col=pal.p[Ip.2])
  legend("bottomright", legend=leglabs(breaksp),fill=pal.p, bty="n")



#plot of the Rates
par(mfrow=c(1,1))
windows()
plot(nc.sids, col=colours[findInterval(nc.sids$ft.SID74, brks,all.inside=TRUE)])
legend(x=c(-84, -80), y=c(33, 34.5), legend=leglabs(brks),fill=colours, bty="n")
text(coordinates(nc.sids),label=as.factor(nc.sids$names),cex=0.55)
invisible(title(main=paste("Rate (Transformed FT) in North Carolina","For the 1974-1978 period")))




