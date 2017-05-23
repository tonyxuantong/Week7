## Part1
## Libraries
library(raster)
## Load data
load("data/GewataB1.rda")
load("data/GewataB5.rda")
load("data/GewataB7.rda")
load("data/GewataB2.rda")
load("data/GewataB3.rda")
load("data/GewataB4.rda")
load("data/vcfGewata.rda")

par(mfrow = c(1, 1)) # reset plotting window
gewata <- brick(GewataB1, GewataB5, GewataB7, vcfGewata)
hist(gewata)

pairs(gewata)
plot()

## compare bands
fit1 <- lm(vcfGewata ~ GewataB1, data=gewata)
fit2 <- lm(vcfGewata ~ GewataB2, data=gewata)
fit3 <- lm(vcfGewata ~ GewataB3, data=gewata)

## stack bands
s1 <- stack(vcfGewata, GewataB1)
s2 <- stack(vcfGewata, GewataB5)
s3 <- stack(vcfGewata, GewataB7)

## fit and summary
fit1 <- lm(s1[1:30] ~ s1[31:60])
summary(fit1)
fit2 <- lm(s2[1:30] ~ s2[31:60])
summary(fit2)
fit3 <- lm(s3[1:30] ~ s3[31:60])
summary(fit3)

## Part 2
## Build a brick containing all data
alldata <- brick(GewataB1, GewataB2, GewataB3, GewataB4, GewataB5, GewataB7, vcfGewata)
names(alldata) <- c("band1", "band2", "band3", "band4", "band5", "band7", "VCF")

## Extract all data to a data.frame
df <- as.data.frame(getValues(alldata))

plot(alldata)

par(mfrow = c(1, 1)) # reset plotting window

ndvi <- overlay(GewataB4, GewataB3, fun=function(x,y){(x-y)/(x+y)})
plot(ndvi)

covs <- addLayer(gewata, ndvi, vcfGewata)
plot(covs)

load("data/trainingPoly.rda")
plot(trainingPoly, add = TRUE)
trainingPoly@data
trainingPoly@data$Class
trainingPoly@data$Code <- as.numeric(trainingPoly@data$Class)
trainingPoly@data
classes <- rasterize(trainingPoly, ndvi, field='Code')
## Plotting
# Define a colour scale for the classes (as above)
# corresponding to: cropland, forest, wetland
cols <- c("orange", "dark green", "light blue")
## Plot without a legend
plot(classes, col=cols, legend=FALSE)
## Add a customized legend
legend("topright", legend=c("cropland", "forest", "wetland"), fill=cols, bg="white")
covmasked <- mask(covs, classes)
plot(covmasked)
## Combine this new brick with the classes layer to make our input training dataset
names(classes) <- "class"
trainingbrick <- addLayer(covmasked, classes)
plot(trainingbrick)
## Extract all values into a matrix
valuetable <- getValues(trainingbrick)
valuetable <- na.omit(valuetable)
valuetable <- as.data.frame(valuetable)
head(valuetable, n = 10)
tail(valuetable, n = 10)
valuetable$class <- factor(valuetable$class, levels = c(1:3))
val_crop <- subset(valuetable, class == 1)
val_forest <- subset(valuetable, class == 2)
val_wetland <- subset(valuetable, class == 3)

## 1. NDVI
par(mfrow = c(3, 1))
hist(val_crop$NDVI, main = "cropland", xlab = "NDVI", xlim = c(0, 1), ylim = c(0, 4000), col = "orange")
hist(val_forest$NDVI, main = "forest", xlab = "NDVI", xlim = c(0, 1), ylim = c(0, 4000), col = "dark green")
hist(val_wetland$NDVI, main = "wetland", xlab = "NDVI", xlim = c(0, 1), ylim = c(0, 4000), col = "light blue")
par(mfrow = c(1, 1))
## 2. VCF
par(mfrow = c(3, 1))
hist(val_crop$VCF, main = "cropland", xlab = "% tree cover", xlim = c(0, 100), ylim = c(0, 7500), col = "orange")
hist(val_forest$VCF, main = "forest", xlab = "% tree cover", xlim = c(0, 100), ylim = c(0, 7500), col = "dark green")
hist(val_wetland$VCF, main = "wetland", xlab = "% tree cover", xlim = c(0, 100), ylim = c(0, 7500), col = "light blue")
par(mfrow = c(1, 1))
## 3. Bands 3 and 4 (scatterplots)
plot(band4 ~ band3, data = val_crop, pch = ".", col = "orange", xlim = c(0, 0.2), ylim = c(0, 0.5))
points(band4 ~ band3, data = val_forest, pch = ".", col = "dark green")
points(band4 ~ band3, data = val_wetland, pch = ".", col = "light blue")
legend("topright", legend=c("cropland", "forest", "wetland"), fill=c("orange", "dark green", "light blue"), bg="white")
## Construct a random forest model
# Covariates (x) are found in columns 1 to 5 of valuetable
# Training classes (y) are found in the 'class' column of valuetable
## Caution: this step takes fairly long!
# but can be shortened by setting importance=FALSE
library(randomForest)
modelRF <- randomForest(x=valuetable[ ,c(1:5)], y=valuetable$class,
                        importance = TRUE)
## Inspect the structure and element names of the resulting model
modelRF
class(modelRF)
str(modelRF)
names(modelRF)
## Inspect the confusion matrix of the OOB error assessment
modelRF$confusion
# to make the confusion matrix more readable
colnames(modelRF$confusion) <- c("cropland", "forest", "wetland", "class.error")
rownames(modelRF$confusion) <- c("cropland", "forest", "wetland")
modelRF$confusion
varImpPlot(modelRF)
## Double-check layer and column names to make sure they match
names(covs)
names(valuetable)
## Predict land cover using the RF model
predLC <- predict(covs, model=modelRF, na.rm=TRUE)
## Plot the results
# recall: 1 = cropland, 2 = forest, 3 = wetland
cols <- c("orange", "dark green", "light blue")
plot(predLC, col=cols, legend=FALSE)
legend("bottomright", 
       legend=c("cropland", "forest", "wetland"), 
       fill=cols, bg="white")
