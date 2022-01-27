setwd("C:/Users/johan/Google Drive/UNI/UPC/Spatial Epidemiology/Task Lattice Data")


##### IMPORT #####
install.packages("maptools")
install.packages("sf")
install.packages("spdep")

library(sf)
library(spdep)
library(maptools)

#### EXERCISE 1 ####
larynx=read.table("larynx_data.txt")
ne <-sf::st_read("NWEngland.shp")

#calculate and plot SMR
larynx$SMR <- larynx$O / larynx$E
hist(larynx$SMR)

#convert a sf to spatial Dataframe
spd <- sf::as_Spatial(ne)
spd$SMR <- larynx$SMR
df <- as.data.frame(ne)
df$geometry <- NULL

## create the SpatialPolygonsDataFrame
ne.lar <- sp::SpatialPolygonsDataFrame(spd, data = df )
summary(ne.lar)
slot(ne.lar,"data")

#add SMR
ne.lar$SMR <-larynx$SMR
slot(ne.lar, "data")

#plot SMR
brks <- round(quantile(ne.lar$SMR, probs=seq(0,1,0.2)), digits=2)   
colours <- c("yellow", "orange2", "red3", "brown", "black")

plot(ne.lar,axes=TRUE,col=colours[findInterval(ne.lar$SMR, brks,all.inside=TRUE)])
legend(290000,400000, legend=leglabs(brks),fill=colours, bty="n",cex=0.5)
title(main=paste("Standardized Morbility Rate (SMR)","\n 1982-1991"))

#find regions with most/ least number of neighbors
#neighbors=sharing geographic border
xxnb <- poly2nb(ne.lar)
plot(ne.lar) 
plot(xxnb, coordinates(ne.lar), add=TRUE, col="blue")
summary.nb(xxnb)

#regions with largest number of neighbors: region 88 with 16 links
cards <- card(xxnb)
sort(cards)
which(cards == max(cards))
maxconts <- which(cards == max(cards))
ne.lar$SMR[maxconts]
fg <- rep("grey", length(cards))
fg[maxconts] <- "red"
fg[xxnb[[maxconts]]] <- "green"
plot(ne.lar, col=fg, axes=TRUE)
title(main="Region with largest number of Neighbors")

#region with lowest number of neighors

#10 least connected regions:
#50 64 65 76 78 94 103 110 128 144 with 2 neighors

minconts <- which(cards == min(cards))
fg <- rep("grey", length(cards))
for(el in minconts){
  fg[el] <- "red"
  fg[xxnb[[el]]] <- "green"
}
fg[xxnb[[maxconts]]] <- "green"
plot(ne.lar, col=fg, axes=TRUE)
title(main="Region with smallest number of neighbors")

#adjacency matrix 
#standardised to sum unity row


b.sids<-nb2listw(xxnb, glist=NULL, style="B",  zero.policy=TRUE)   #Spatial weights for neighbours lists
w.sids<-nb2listw(xxnb, glist=NULL, style="W",  zero.policy=TRUE)   #Spatial weights for neighbours lists

b.sids$weights
summary(unlist(b.sids$weights))
summary(sapply(b.sids$weights,sum))


w.sids$weights
summary(unlist(w.sids$weights))
summary(sapply(w.sids$weights,sum))



