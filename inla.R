
################################################################
#################### --- FIT THE MODELS --- ####################
################################################################

#Libraries
library(sp)
library(geoR)
library(dismo)
library(hSDM)
library(rgdal)
library(spdep)
library(fields)
library(raster)
library(maptools)
library(gridExtra)
library(ggplot2)
library(rworldmap)
require(rworldxtra)
library(openxlsx)
library(INLA)
library(rgeos)
#library(plotKML)


#####################################
### --- Set working directory --- ###
#####################################
dir1<- "~/Documentos/TESIS_DOCTORAL/CURSOS/IMPARTIDOS/2018/CursoSDMBarca"
setwd(paste0(dir1, "/Day_4/Practice/Hake2"))

#####################################
### --- Read the data         --- ###
#####################################
data <- readRDS("data.rds")
str(data)


#####################################
### --- Boundaries Barcelona --- ####
#####################################
### --- Define the polygon --- ###
### --- Spain --- ###
spain <- getData('GADM',country="ESP",level=0)

### --- Catalunya --- ###
ext<-extent(-2,4,35,44)
cat <- crop(spain, ext) #Only Catalunya
plot(cat)
points(data[,2:3], pch=20)

### --- Define the polygon around the data --- ###
xym<- as.matrix(data.frame(x = c(min(data$Lon), min(data$Lon) + 0.8, max(data$Lon) + 0.2, min(data$Lon) ), 
                           y = c(39.5, 39.5, max(data$Lat), max(data$Lat))))
p = Polygon(xym)
ps = Polygons(list(p),1)
sps = SpatialPolygons(list(ps))
plot(sps)
points(data[,2:3], pch=20)

#Crop with Catalunya
cat_rec<-crop(cat, sps)
plot(cat_rec)
points(data[,2:3], pch=20)
proj4string(sps)<-proj4string(cat_rec)
plot(cat_rec)

### --- Define square to plot --- ###
### --- square --- ###
xym2<- as.matrix(data.frame(x = c(min(data$Lon), max(data$Lon) + 0.2, max(data$Lon) + 0.2, min(data$Lon) ), 
                            y = c(39.5, 39.5, max(data$Lat), max(data$Lat))))
p2 = Polygon(xym2)
ps2 = Polygons(list(p2),1)
sps2 = SpatialPolygons(list(ps2))
plot(sps2)
plot(cat_rec, add=TRUE)
points(data[,2:3], pch=20)

### --- Select the polygon which contains the data --- ###
coast <- gDifference(sps, cat_rec )
plot(coast)
points(data[,2:3], pch=20)

### --- Plot the data --- ###
#pdf("data.pdf")
plot(sps2)
points(data[,2:3][which(data$presence==1),], col="red", pch=20)
points(data[,2:3][which(data$presence==0),], col="blue", pch=20)
legend(1.5, 40, legend=c("presence", "absence"), col=c("red", "blue"), pch=20)
plot(cat_rec, add=TRUE, col="gray")
#dev.off()


#######################################
### --- Environmental variables --- ###
#######################################
#Get the file names
files<-(list.files(paste0(dir1, "/Day_4/Practice/Hake2/predictors"), full.names=T, pattern=".asc"))#change directory
predictors <- stack(files)#bioclimatic variables 
names(predictors) <- c("calcite","chlomean","nitrate","ph","phos","salinity","silicate","sstmean")

### --- Crop predictors only for our subset
ext<-extent(-2,4,35,44)
predictors<-crop(predictors,ext)

######################################
### --- Standardize predictors --- ###
######################################
predictors2<-scale(predictors)
round(apply(values(predictors2), 2, summary), 4)

#######################################################
### --- Create the datasets with the predictors --- ###
#######################################################
data<-cbind(data, extract(predictors2, as.matrix(data[,2:3])))




###########################################################################
###########################################################################
####################### --- ESTIMATION --- ################################
###########################################################################
###########################################################################
### --- Built the mesh --- ###
#### create a domain of the study area using inla.nonconvex.hull()
boundary=inla.nonconvex.hull(as.matrix(data[,2:3]))
mesh<-inla.mesh.2d(boundary=boundary, max.edge=c(0.07, 0.3),
                   cutoff=0.06,  offset=c(-0.1, -0.3))

#pdf("mesh.pdf", width = 10, height = 10)
  plot(mesh)
  plot(spain, add=TRUE, col="gray")
  plot(mesh, add=TRUE)
  points(data[,2:3][which(data$presence==1),], col="red", pch=20)
  points(data[,2:3][which(data$presence==0),], col="blue", pch=20)
#dev.off()
####################################################


### --- Definition of the spde --- ####
spde <- inla.spde2.matern(mesh)

### --- Matrix which link data with the mesh --- ###
A.est <- inla.spde.make.A(mesh, loc=cbind(data$Lon, data$Lat))

### --- inla.stack to stimate --- ###
stk.est<-inla.stack(data=list(y=data$presence),
                    A=list(A.est, 1),
                    effects=list(spatial=1:spde$n.spde,
                                 data.frame(beta0=1, data)),
                    tag='est')

################################################################
##### --- Selection of variables      --- ######################
##### --- Response variable: binomial --- ######################
##### --- Barriers: yes               --- ######################
################################################################
source("Bdiclcpomodel_stack.R")

### --- Define the variable to select the best model --- ###

variables <- c("calcite", "chlomean", "sstmean", "f(spatial, model=spde)")

# variables <- c("calcite", "chlomean", "nitrate",  "ph", "phos", "salinity", "silicate", "sstmean", 
#               "f(spatial, model=spde)")

### --- Response variable --- ###
resp=data$presence

### --- Call the function --- ###
models_bin<-Bdiclcpomodel_stack(resp=resp, variables=variables, datos=inla.stack.data(stk.est), n=20,
                                          family="binomial",
                                          control.predictor=list(compute=TRUE, A=inla.stack.A(stk.est)),
                                          control.compute = list(config=TRUE, dic=TRUE, cpo=TRUE, waic=TRUE),
                                          num.threads=3,
                                          control.inla=list(strategy="gaussian"),
                                          verbose=FALSE)


# saveRDS(models_bin, "models/models_bin.rds")
models_bin<- readRDS("models/models_bin.rds")

models_bin$`Modelos según dic`[1:10,]
models_bin$`Modelos según waic`[1:10,]
models_bin$`Modelos según lcpo`[1:10,]



### --- fitting the model --- ###
formula.1 <- y~-1 + beta0  + sstmean +  f(spatial,model=spde)
model.est <- inla(formula.1, 
                  data=inla.stack.data(stk.est), family="binomial" ,
                  control.compute=list(dic=TRUE,cpo=TRUE, waic=TRUE, return.marginals=TRUE), 
                  control.predictor=list(A=inla.stack.A(stk.est), compute=TRUE, 
                                         quantiles=c(0.025, 0.25, 0.5, 0.75, 0.975)),
                  control.family=list(quantile=c(0.025)),
                  #control.inla=list(strategy = "laplace"),
                  num.threads = 3,
                  verbose=FALSE)

### --- Save the model in a .rds ---###
#saveRDS(model.est, "models/model_est.rds")
model.est<- readRDS("models/model_est.rds")
summary(model.est)

model.est$logfile #check the logfile
summary(model.est)

1-inla.pmarginal(0, model.est$marginals.fixed$sstmean) #P(beta1>0)=0.98 So relevant

### --- Plot posteriors --- ###
plot(model.est$marginals.fixed$sstmean, type="l",main="Posterior distribution of beta1")
abline(v=0, col="red", lwd=2)

####################################################
##### --- Posterior distribution hyperpars --- #####
####################################################
### --- Postprocess --- ###
spde.result = inla.spde2.result(model.est, "spatial", spde, do.transform=TRUE)

range<-inla.emarginal(function(x) x, spde.result$marginals.range.nominal[[1]]) 

### --- Check the range is smaller than the offset --- ###
range < max(diff(range(data[,1])), diff(range(data[,2])))*0.40 #Yes!!!

### --- Plot --- ###
par(mfrow=c(2,2), mar=c(3,3,1,0.5)) 
plot(spde.result$marginals.range.nominal[[1]], type='l') 



####################################################
### --- Interpolate the posterior mean and sd --- ##
####################################################
### --- plot in a grid m X m --- ##

### --- Customize the grid to predict --- ###
bbox(coast)
(dxy <- apply(bbox(coast),1, diff))
(r <- dxy[1]/dxy[2])
m<-150
proj.grid.mat <- 
  inla.mesh.projector(mesh, 
                      xlim=bbox(coast)[1,],
                      ylim=bbox(coast)[2,] ,
                      dims=c(r, 1)*m)

plot(coast)
points(proj.grid.mat$lattice$loc, pch=20, cex=0.5)

### --- clean (set NA to the values outside boundary) --- ###
ov <- over(SpatialPoints(proj.grid.mat$lattice$loc, coast@proj4string),
           coast)

### --- check grid points inside the map --- ###
i.map <- is.na(ov)

### Plot the points where we will predict ###
par(mar=c(0,0,0,0))
plot(sps)
points(proj.grid.mat$lattice$loc[!i.map,], col="red", cex=0.2)
points(proj.grid.mat$lattice$loc[i.map,], col="blue", cex=0.2)

### --- consider only those inside map --- ###
proj.grid.mat$lattice$loc[i.map, ]

### --- Project the values of the mean and sd of the spatial effect --- ###
mean.g <- inla.mesh.project(proj.grid.mat, model.est$summary.random$spatial$mean)
sd.g <- inla.mesh.project(proj.grid.mat, model.est$summary.random$spatial$sd)
quantile_0.025 <- inla.mesh.project(proj.grid.mat, model.est$summary.random$spatial$`0.025quant`)
quantile_0.975 <- inla.mesh.project(proj.grid.mat, model.est$summary.random$spatial$`0.975quant`)


sd.g[i.map] <- mean.g[i.map] <- quantile_0.025[i.map] <- quantile_0.975[i.map] <- NA


#pdf("spatial_effect.pdf", width=10, height = 10)
### --- Spatial effect --- ###
par(mfrow=c(2,2))
par(mar=c(2,3,3,6))


### --- Posterior mean --- ###
plot(sps2, col="gray", main="Spatial mean")
image.plot(proj.grid.mat$x, 
           proj.grid.mat$y,
           mean.g, add=TRUE)
plot(cat_rec, add=TRUE)

### --- Posterior sd --- ###
plot(sps2, col="gray", main="Spatial sd")
image.plot(proj.grid.mat$x, 
           proj.grid.mat$y,
           sd.g, add=TRUE)
plot(cat_rec, add=TRUE)


### --- Posterior q0.025 --- ###
plot(sps2, col="gray", main="Spatial q0.025")
image.plot(proj.grid.mat$x, 
           proj.grid.mat$y,
           quantile_0.025, add=TRUE)
plot(cat_rec, add=TRUE)

### --- Posterior q0.975 --- ###
plot(sps2, col="gray", main="Spatial q0.975")
image.plot(proj.grid.mat$x, 
           proj.grid.mat$y,
           quantile_0.975, add=TRUE)
plot(cat_rec, add=TRUE)
#dev.off()

###########################################################################
###########################################################################
######################### --- Prediction --- ##############################
###########################################################################
###########################################################################
### --- Matrix which link the mesh with coordinates to predict--- ###
A.pred <- inla.spde.make.A(mesh, loc=proj.grid.mat$lattice$loc[!i.map, ])

### --- Stack to predict --- ###
stk.pred <- inla.stack(data=list(y=NA),
                       A=list(A.pred, 1), 
                       effects=list(spatial=1:spde$n.spde,
                                    data.frame(beta0 = 1, 
                                               extract(predictors2, 
                                                      proj.grid.mat$lattice$loc[!i.map, ]))),
                       tag='pred')

stk <- inla.stack(stk.est, stk.pred)


#### --- model --- ###
model.pred <- inla(formula.1, 
                   data=inla.stack.data(stk), family="binomial",
                   control.predictor=list(A=inla.stack.A(stk), compute=TRUE, link=1), #link:link is a vector of
                   #length given by the size of the response variable with values 1 if the corresponding
                   #data is missing and NA otherwise
                   control.inla=list(strategy = "simplified.laplace"), # Strategy
                   control.mode=list(theta=model.est$mode$theta, restart=TRUE), #Mode 
                   control.results=list(return.marginals.random=FALSE,
                                        return.marginals.predictor=FALSE), # Avoid some marginals
                   num.threads = 3,
                   verbose=FALSE)

#saveRDS(model.pred, "models/model_pred.rds")
model.pred<-readRDS("models/model_pred.rds")

####################################################
########### --- Plot predictions --- ###############
####################################################
### index for the prediction data
idx <- inla.stack.index(stk, 'pred')$data

summary(model.pred$summary.fitted.val$mean[idx])

### --- Organize probabilities into a matrix to visualize --- ###
prob.mean <- prob.sd <- prob.0.025<- prob.0.975 <- matrix(NA, proj.grid.mat$lattice$dims[1],
                               proj.grid.mat$lattice$dims[2])
prob.mean[!i.map] <- c(model.pred$summary.fitted.val$mean[idx])
prob.sd[!i.map] <- c(model.pred$summary.fitted.val$sd[idx])
prob.0.025[!i.map] <- c(model.pred$summary.fitted.val$`0.025quant`[idx])
prob.0.975[!i.map] <- c(model.pred$summary.fitted.val$`0.975quant`[idx])

#### --- plot --- ###

#pdf("predictive.pdf", width=10, height = 10)
### --- Spatial effect --- ###
par(mfrow=c(2,2))
par(mar=c(2,3,3,6))

### --- posterior predictive mean --- ###
plot(sps2, col="gray", main="Predictive mean")
image.plot(proj.grid.mat$x, 
           proj.grid.mat$y,
           prob.mean, add=TRUE, zlim=c(0,1))

plot(cat_rec, add=TRUE)

# points(data[,2:3][which(data$presence==1),], col="red", pch=20)
# points(data[,2:3][which(data$presence==0),], col="blue", pch=20)

### --- posterior predictive sd --- ###
plot(sps2, col="gray", main="Predictive sd")
image.plot(proj.grid.mat$x, 
           proj.grid.mat$y,
           prob.sd, add=TRUE)
# 
# points(data[,2:3][which(data$presence==1),], col="red", pch=20)
# points(data[,2:3][which(data$presence==0),], col="blue", pch=20)

plot(cat_rec, add=TRUE)


### --- posterior predictive q0.025 --- ###
plot(sps2, col="gray", main="Predictive q0.025")
image.plot(proj.grid.mat$x, 
           proj.grid.mat$y,
           prob.0.025, add=TRUE, zlim=c(0,1))
plot(cat_rec, add=TRUE)

### --- posterior predictive q0.975 --- ###
plot(sps2, col="gray", main="Predictive q0.975")
image.plot(proj.grid.mat$x, 
           proj.grid.mat$y,
           prob.0.975, add=TRUE, zlim=c(0,1))
plot(cat_rec, add=TRUE)

#dev.off()









