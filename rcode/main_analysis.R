##'# Procedimento geral de classificacao de imagens Landsat
#'# Paulo E. Cardoso 27/07/2015
#'# Classificação da Cena Landsat de 2014 LC82040522014030LGN00
#'# para correcao da are ade mangal na ilha de bolama

#'# Load/Install Packages ##########################################################
kpacks <- c('raster', 'sp', 'rgdal', 'rgeos', 'ggplot2'
            ,'dplyr', 'foreign', 'reshape2', 'lubridate', 'scales'
            ,'randomForest', 'caret', 'glcm', 'tidyr', 'e1071'
            ,'mlbench', 'pROC', 'rpart', 'caretEnsemble', 'gbm'
            ,'RStoolbox')
new.packs <- kpacks[!(kpacks %in% installed.packages()[ ,"Package"])]
if(length(new.packs)) install.packages(new.packs)
lapply(kpacks, require, character.only=T)
remove(kpacks, new.packs)
#' To test
#' # browseURL() # Load a given URL into a WWW browse
#' # https://everydropr.wordpress.com/2012/09/22/how-to-efficiently-mask-raster-layer-in-r/
#' # http://web.nateko.lu.se/timesat/timesat.asp?cat=6 # TIMESAT NDVI profile
#'##################################################################################
save.image("D:/Programacao/RLandsat/Run_pntc_mangrove_main0.RData")
load('D:/Programacao/RLandsat/R/tascap1.rda')
#'# Projections ####################################################################
p.utm28n <- CRS("+init=epsg:32628") # 32733 UTM 28N Landsat Images
p.wgs84 <- CRS("+init=epsg:4326") # WGS84 Long Lat

#' Parameters ######################################################################
#' # LC82040522015145LGN00
#' DATE_ACQUIRED = 2015-05-25
#' SUN_AZIMUTH = 65.24042599
#' SUN_ELEVATION = 64.56521845
#' CLOUD_COVER = 4.19
#' #############################################
#' LC82040522015081LGN00
#' DATE_ACQUIRED = 2015-03-22
#' SUN_AZIMUTH = 110.15370427
#' SUN_ELEVATION = 60.99211548
#' CLOUD_COVER = 0.05
#' #############################################
#' # LC82040522015049LGN00
#' DATE_ACQUIRED = 2015-02-18
#' SUN_AZIMUTH = 128.07522029
#' SUN_ELEVATION = 53.24490938
#' CLOUD_COVER = 0.02
#' #############################################
#' # LC82040522015001LGN00
#' DATE_ACQUIRED = 2015-01-01
#' SUN_AZIMUTH = 143.56592769
#' SUN_ELEVATION = 47.09703778
#' CLOUD_COVER = 0.03
#' #############################################
#' LC82040522014318LGN00
#' DATE_ACQUIRED = 2014-11-14
#' SUN_AZIMUTH = 144.52272103
#' SUN_ELEVATION = 53.58396593
#' CLOUD_COVER = 0.61
#' #############################################

scene1 <- 'LC82040522014030LGN00' # 2014-11-14 <- ok
#' Working Folders -----------------------------------------------------------------
dir_work <- 'D:/Sig/Raster/landsat'
dir_landsat <- 'dos1'
#dir_landsat <- 'tif'
#dir.data <- 'D:/Sig/Bissau/Cacheu/sig/vetor'
dir_sub <- 'sub'
#path_to_dos1 <- file.path(dir.work, dir.landsat)

#' Pick a Reference band -----------------------------------------------------------
pick_scene <- scene1
band <- raster(file.path(dir_work, pick_scene, dir_landsat,
                         grep(".tif$",
                              list.files(file.path(dir_work, pick_scene, dir_landsat),
                                         all.files = F),
                              ignore.case = TRUE, value = TRUE)[1])
               ,values = F, package = "raster")
values(band) <- 1
plot(band)
#' Study Site - Will be used for mask scenes ---------------------------------------
frames <- readOGR(dsn = 'D:/Sig/Bissau/Mangal'
              , layer = 'frames_mangais7_1')
sp::is.projected(frames)
if(is.na(sp::proj4string(frames))) proj4string(frames) <- p.utm28n
#spTransform(ae, p.utm33n)
ae <- frames[frames$gr==5, ]
ae_ext <- extent(ae)
ae_lext <- ae_ext + 50 # large extent for texture calculation
ae1 <- as(ae_lext, 'SpatialPolygons')
proj4string(ae1) <- sp::proj4string(ae)
ae1wgs84 <- spTransform(ae1, CRS=p.wgs84)

#' Centroide da area de estudo - obtencao de dados auxiliares
ae_ctd <- coordinates(rgeos::gCentroid(ae1wgs84))

#' Function pack -------------------------------------------------------------------
func <- grep("^f_", ls(), perl=TRUE, value = T)
save(list=func, file='D:/Programacao/RLandsat/functions5.RData')
load('D:/Programacao/RLandsat/functions5.RData')

#' Import your ROI for image cropping/extract --------------------------------------
#' Masks 
mask_ae <- f_createRoiMask(maskpoly = ae, maskv = 1, band = band)
mask_ae1 <- f_createRoiMask(maskpoly = ae1, maskv = 1, band = band)
plot(mask_ae, axes = T)
#'Stack DOS1  ----------------------------------------------------------------------
stk_dos1 <- f_stkDOS1(roi = ae, fpath = file.path(dir_work, scene1, dir_landsat))
plot(stk_dos1)
dem <- getData('SRTM', lon=ae_ctd[1], lat=ae_ctd[2])
dem_1 <- raster::crop(dem, ae1wgs84)
dem_ae1 <- raster::projectRaster(dem_1, mask_ae)#crs=p.utm28n)
#dem_ae1 <- raster::resample(dem_ae1, mask_ae)
compareRaster(dem_ae1, mask_ae)
#stk_mask <- f_applmask(stk = stk_dos1, mask = mask_ae)

#' PCA on Bands --------------------------------------------------------------------
pca_obj1 <- f_pca(stk = stk_dos1, corr = F, comps = 3)[[2]]
names(pca_obj1) <- c('comp1_1','comp2_1', 'comp3_1') 

#' TasseledCap ----------------------------------------------------------------
tascap <- tasseledCap(stk_dos1, sat = "Landsat8OLI")
names(tascap) <- c('brigthness', 'greenness', 'soilwetness')
plot(tascap)

writeRaster(tascap, filename='D:/Sig/Bissau/bijagos/sub/tasscap'
            ,bylayer=T
            ,format="GTiff"
            ,overwrite=T)
#' Veg Index -----------------------------------------------------------------------
#' RSToolbox Vegetation Index
vegindx <- RStoolbox::spectralIndices(stk_dos4
                                      ,blue='b2_ae', green = 'b3_ae'
                                      , red = "b4_ae"
                                      , nir = "b5_ae", swir1 = 'b6_ae'
                                      , swir2 = 'b7_ae'
                                      , indices = c('EVI', 'NDVI'))

plot(vegindx)

stk_veg <- f_vegind(stk = stk_dos1)

plot(stk_veg)
#' Stack para a construcao do Modelo -----------------------------------------------
stk_modelRF <- raster::stack(
                            pca_obj1
                           ,stk_veg[['nbr']]
                           ,stk_veg[['ndvi']]
                           ,stk_veg[['ndmi']]
                           #,tascap
                           #,tascap1
                           #, dem_ae1
                           )
names(stk_modelRF)
#' Crop stack com ae: atencao ao stk !!
stk_modelRF <- raster::crop(stk_modelRF, ae)
#write.table(as.data.frame(names(stk_model)), file='clipboard', sep = '\t')
#writeRaster(stk_model
#            ,filename='D:/Sig/Bissau/Cacheu/usosolo/sub/vars/stk_model2000'
#            ,format="GTiff"
#            ,overwrite=T
#            ,bylayer=T)
#' Signature Develop ---------------------------------------------------------------
#' fora: "floresta_palmar" ,"lala_palmar"
train <- rgdal::readOGR(dsn = file.path('D:/Sig/Bissau/Mangal')
                        ,layer = 'areas_treino', stringsAsFactors = F)

as.data.frame(table(train$classe_n2))
#' Subset training classes ---------------------------------------------------------
#classes <- c("pampam","lala","floresta_aberta", "savana_arborea"
#             ,"palmar","caju","aberto", "mangal"
#             ,"bolanha", "agua", 'floresta_palmar')
#train <- train[train@data$classe %in% classes, ]
datasets <- f_sign(stk = stk_modelRF, train = train, class = 'classe_n2', subsets = F)
head(datasets)
set.seed(998)
inTraining <- createDataPartition(datasets$classe, p = .80, list = FALSE)
training <- datasets[ inTraining, ]
testing  <- datasets[-inTraining, ]
as.data.frame(table(training$classe))
head(training)

#' Random Forest model with caret --------------------------------------------------
fitControl <- caret::trainControl(
  ## 5-fold CV
  method = "repeatedcv",
  number = 5,  #5
  ## repeated n times
  repeats = 2) #5
set.seed(825)
rfFit1 <- caret::train(factor(classe) ~ ., data = training
                       ,method = "rf"
                       ,tuneLength = 5
                       #,tuneGrid=grid_c
                       ,metric = "Accuracy"
                       ,trControl = fitControl
                       ## This last option is actually one
                       ## for gbm() that passes through
                       ,verbose = FALSE
                       ,allowParallel=TRUE)
rfFit1
#' Predict Random Forest------------------------------------------------------------
rf_mod <- raster::predict(stk_modelRF, rfFit1,  type='raw', progress='window')
rf_mod_main <- raster::mask(rf_mod, ae)

#' RF Model Outputs ----------------------------------------------------------------
print(rfFit1$finalModel)
varImp(rfFit1)
plot(varImp(rfFit1))
plot(rfFit1)
plot(rfFit1$finalModel)
#library('rattle')
#fancyRpartPlot(rfFit1$finalModel)

#' Confusion matrix ----------------------------------------------------------------
rf_test <- raster::predict(rfFit1, newdata = testing)
matrizc<- caret::confusionMatrix(data = rf_test, testing$classe)
matrizc

#' Export Random Forest surface
writeRaster(rf_mod_main
            , filename=file.path('D:/Sig/Raster/landsat/mangal'
                                 ,'rforest_4_LC82040522014030LGN00_5_2.tif')
            , format="GTiff", datatype="INT1U", overwrite=TRUE)

#' Write KML -----------------------------------------------------------------------
#rf_mod_pntc_wgs84 <- projectRaster(rf_mod_pntc, crs=p.wgs84,method='ngb')
#KML(rf_mod_pntc_wgs84, file = file.path('D:/Sig/Bissau/Cacheu/usosolo/sub'
#                                        , 'rforest10cl_4.kml')
#    , col = c('#ffffe6', '#a6cee3', '#dbfd3f', '#c13755'
#              ,'#4bb455', '#50e230', '#5115ba', '#3b6b38'
#              ,'#a36e48', '#ffae2c', '#d1cfc3'), blur=40
#    , overwrite=T, zip = 'C:/Program Files/7-Zip/7zFM.exe')

#' SVM model -----------------------------------------------------------------------
ctrlSvm <- trainControl(method = "cv", savePred=T, classProb=T)
svm.model <- train(factor(classe) ~ ., data = training
                   , method = "svmLinear", trControl = ctrlSvm)
print(svm.model$finalModel)

#' Predict SVM ---------------------------------------------------------------------
svm_mod <- predict(stk_model, svm.model, progress="text")    
writeRaster(svm_mod
            , filename=file.path('D:/Sig/Bissau/Cacheu/usosolo/sub'
                                 , 'svm10cl_3.tif')
            , format="GTiff", datatype="INT1U", overwrite=TRUE)

#' gbm fitting --------------------------------------------------------------------
set.seed(123)
fitControl <- trainControl(method = 'cv', number = 5, summaryFunction=defaultSummary)
gridctrl <- expand.grid( n.trees = seq(200,500,100)
                         , interaction.depth = c(30)
                         , shrinkage = c(0.1)
                         , n.minobsinnode = 20)

gbmFit1 <- caret::train(factor(classe) ~ ., data = training
                        , method = 'gbm'
                        ,distribution='multinomial'
                        #, trControl=fitControl
                        ,tuneGrid=gridctrl
                        #,metric='RMSE'
                        ,maximize=FALSE)
gbmFit1
trellis.par.set(caretTheme())
plot(gbmFit1)
plot(varImp(gbmFit1))

#' Confusion matrix ----------------------------------------------------------------
gbm_test <- predict(gbmFit1, newdata = testing)
matrizc<- confusionMatrix(data = gbm_test, testing$classe)
matrizc
#' Predict Random Forest------------------------------------------------------------
gbm_mod <- predict(stk_model, gbmFit1,  type='raw', progress='window')
gbm_mod_pntc <- raster::mask(gbm_mod, pntc_b)

#' Export GBM
writeRaster(gbm_mod_pntc
            , filename=file.path('D:/Sig/Bissau/Cacheu/usosolo/sub'
                                 , 'gbm10cl_3.tif')
            , format="GTiff", datatype="INT1U", overwrite=TRUE)

#' Boosted tree
ctr <- trainControl(method = "cv", number = 10)
boostFit1 <- train(factor(class) ~ ., data = training,
                   method='bstTree',
                   preProc=c('center','scale'),
                   trControl=ctr)

#' MSE 
mean((rf_mod - testing$classe)^2)
