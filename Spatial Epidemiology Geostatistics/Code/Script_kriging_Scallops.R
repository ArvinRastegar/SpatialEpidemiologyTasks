library(geoR)

# setwd("F:\\Master d'Estadistica i IO\\Geoestad?stica")
# setwd("D:/Master d'Estadistica i IO/Geoestadística/Diapositives curs/Capitol 5")

load("scallops.Rdata")

#############################################################
# Datos: Scallops.
##############################################################


###Tendency
#Remember 
# lm.scp.rot<-lm(lgcatch ~newx+poly(newy,3),data=scall.rot)

#Preditions locations
rnx <- range(geoscal.rot$coords[,1]);
rny <- range(geoscal.rot$coords[,2]);
newx.grid <- seq(rnx[1],rnx[2],l=51);
newy.grid <- seq(rny[1],rny[2],l=51);
dsgr.grid <- expand.grid(newx=newx.grid, newy=newy.grid);
plot(dsgr.grid)


########################################################
#Options available implement the following types of kriging:
#"SK" (simple kriging), "OK" (ordinary kriging).

#Kriging with external trend and universal kriging can be defined 
#setting type.krige = "OK" 
#and specifying the trend model using the arguments trend.d and trend.l.

#krige.conv(geodata, coords=geodata$coords, data=geodata$data,
#           locations, borders, krige, output)
#
#krige.control(type.krige = "ok", trend.d = "cte", trend.l = "cte",
#              obj.model = NULL, beta, cov.model, cov.pars, kappa,
#              nugget, micro.scale = 0, dist.epsilon = 1e-10, 
#              aniso.pars, lambda)

#Maximum estimation of variogram parameters
#Attention we don't take into account the trend of the data
#Not use
lk1.1 <- likfit(geoscalresidus.rot ,cov.model = "exponential", ini =c(3,0.3),
                fix.nugget = F, nugget =1.41 ,lik.method = "REML");

#Take into account the trend of tha data
#likfit(trend = "cte) 
#or "1st"the mean is assumed to be a first order polynomial on the coordinates:
#mu(x) = beta0 + beta1*x1 + beta2*x2.

lk1 <- likfit(geoscal.rot ,cov.model = "exponential", ini =c(3,0.3),
              fix.nugget = F, nugget =1.41 ,lik.method = "REML",trend=~coords[,1]+poly(coords[,2],3));



kc1<-krige.conv(geoscal.rot,coords= geoscal.rot$coords,
                data=geoscal.rot$data,locations= dsgr.grid,
                krige =krige.control(type.krige="OK",obj.m = lk1,trend.l= ~coords[,1]+poly(coords[,2],3), 
                                     trend.d=~coords[,1]+poly(coords[,2],3)))
attributes(kc1)

image(kc1, main="Universal kriging estimates")
image(kc1, val=sqrt(kc1$krige.var), main="Universal kriging std. errors",col=gray(seq(1,0.1,l=30)))

contour(kc1,filled=TRUE,coords.data=geoscal.rot$coords,color=heat.colors)
contour(kc1,val=sqrt(kc1$krige.var),filled=TRUE,coords.data=geoscal.rot$coords,color=cm.colors)

### Convex hull
##Computes the subset of points which lie on the convex hull of the set of points specified.
dg.chull <- chull(geoscal.rot$coords[,1], geoscal.rot$coords[,2])
dg.poly <- list(x = geoscal.rot$coords[,1][dg.chull], y = geoscal.rot$coords[,2][dg.chull]);
inside <- point.in.polygon(dsgr.grid[,1], dsgr.grid[,2], dg.poly$x, dg.poly$y);

image(kc1,coords=geoscal.rot$coords,border=dg.poly )
contour(kc1,filled=TRUE,coords.data=geoscal.rot$coords,color=terrain.colors,border=dg.poly)


##cros validations.

xv <- xvalid(geoscal.rot, model =lk1.1, reestimate = F);
names(xv)

VC1 <- mean(xv$error/sqrt(xv$krige.var));
VC2 <- sqrt(mean((xv$error/sqrt(xv$krige.var))^2));
VC3 <- sqrt(mean(xv$error^2));
VC1;VC2;VC3



####variogram GAUSSIAN

lk2 <- likfit(geoscal.rot,trend=~geoscal.rot$coords[,1]+poly(geoscal.rot$coords[,2],3),cov.model = "gaussian", ini = c(3,0.3),
              fix.nugget = F, nugget =1.41 ,fix.kappa = FALSE, kappa=1,lik.method = "REML");
lk2

kc2<-krige.conv(geoscal.rot,coords= geoscal.rot$coords,
                data=geoscal.rot$data,locations= dsgr.grid,
                krige =krige.control(type.krige="ok",trend.l= ~coords[,1]+poly(coords[,2],3), 
                                     trend.d=~coords[,1]+poly(coords[,2],3),obj.m = lk2))

contour(kc2,filled=TRUE,coords.data=geoscal.rot$coords,color=heat.colors)

contour(kc2,val=sqrt(kc2$krige.var),filled=TRUE,coords.data=geoscal.rot$coords,color=cm.colors)

##cros validations variogram Gaussian

xv.2 <- xvalid(geoscal.rot, model =lk2, reestimate = F);
names(xv.2)

VC1.2 <- mean(xv.2$error/sqrt(xv.2$krige.var));
VC2.2 <- sqrt(mean((xv.2$error/sqrt(xv.2$krige.var))^2));
VC3.2 <- sqrt(mean(xv.2$error^2));
VC1.2;VC2.2;VC3.2
#Compare with the variogram exponential
VC1;VC2;VC3



