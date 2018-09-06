library("zoo")
library("reshape")
library("plotrix")
library("relaimpo")
library("MASS")
library("car")
library("ggplot2")
library("multcomp")
library("class")
library("Hmisc")
library("scatterplot3d")
library("reshape")
library("fitdistrplus")
library("compare")
library ("truncreg")
library("truncdist")
library("msm") ## for dtnorm
library("actuar")
library("moments")
library("vegan")
library("stats")


##Get the data##
###TREES FOR PSP AND TSP
#######################################################################################################################################
####Do separately the tree and sapling PSP
TreePSP <- read.csv("C:\\Users\\psladmin\\Documents\\ModelOutput\\Data\\Tree_list.csv",colClasses = "character" )
SaplingPSP <- read.csv("C:\\Users\\psladmin\\Documents\\ModelOutput\\Data\\Yos_result2.csv",colClasses = "character", header=FALSE )
Tr<-subset(TreePSP, select=c("ID_PEP_MES", "ESSENCE","DBH"))
length(1:dim(Tr)[1])

names(SaplingPSP)[names(SaplingPSP)=="V2"] = "ID_PEP_MES"
names(SaplingPSP)[names(SaplingPSP)=="V3"] = "ESSENCE"
names(SaplingPSP)[names(SaplingPSP)=="V4"] = "DBH"
SaplingPSP$V5 <- NULL
SaplingPSP$V1 <- NULL
head(SaplingPSP)
length(1:dim(SaplingPSP)[1])
PS<-rbind(SaplingPSP,Tr)
head(PS)
length(1:dim(PS)[1])

###############Temporal###
TreeT<- read.table("C:\\Users\\psladmin\\Documents\\ModelOutput\\Data\\TemporalTrees.txt",
                   colClasses = "character",sep = ",", quote="\"")
length(1:dim(TreeT)[1])
head(TreeT)
TreeT$DBH<-as.numeric(TreeT$DBH)
TreeT<-subset(TreeT, select=c("ID_PET_MES", "ESSENCE","DBH"))
names(TreeT)[names(TreeT)=="ID_PET_MES"] = "ID_PEP_MES"
############merge PSP trees with Temporal trees
Tree<-rbind(PS,TreeT)
head(TreeT)
length(1:dim(Tree)[1])
Tree$DBH<-as.numeric(Tree$DBH)
Tree$ESSENCE<-as.factor(Tree$ESSENCE)
Tree<-Tree[-which(Tree$DBH > 90),]
Tree<-na.omit(Tree)
#write.csv(Tree,file ="C:\\Users\\yomiq\\Documents\\YOSDATA\\LAVAL\\ModelOutput2\\Stratified_trees.txt")


#### INVENTORY PLOTS PATHWAYS###############################################

PoolPlots<-read.csv("C:\\Users\\psladmin\\Documents\\ModelOutput3ClimateChange\\InventoryPlotsPathway.csv",colClasses = "character", header=TRUE,sep = ",")
str(PoolPlots)
PoolPlots[,24:34] <- sapply(PoolPlots[,24:34],as.numeric)
length(PoolPlots[,1])  ##4882...3714
PoolPlots2<-PoolPlots[which(PoolPlots$stands== "EnEn"),]
length(1:dim(PoolPlots2)[1])##3714



###GET INTENSITY FILE FOR THE C-2 FUEL TYPE#####################################
Intensities<-read.csv("C:\\Users\\psladmin\\Documents\\ModelOutput3ClimateChange\\Data\\Sopfeu.csv",colClasses = "character", header=TRUE,sep = ",")
range(Intensities$INT)
Intensities<- Intensities[which(Intensities$JA== "1"),]
length(Intensities[,1])##1563
Intensities<- Intensities[which(Intensities$Comb== "C2"),]
Intensities<- Intensities[which(Intensities$Domaine== "6"),]
Intensities<-Intensities[which(Intensities$SupFin> 0),]
length(Intensities[,1])##1112
###Apply Catchpole and size weighted
knownpoints<-data.frame(x<-c(0.04,0.1,0.13,0.16,0.21,0.28,0.37,0.48,0.63,0.7,0.8,0.9,1),
                        y<-c(1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.18,0.15,0.1,0.05))
aim<-runif(24000,0,1)
yos<-approx(knownpoints$x,knownpoints$y, xout=aim, rule=1)
ProportionalIntensities<-yos$y
ProportionalIntensities<-na.omit(ProportionalIntensities)

PI<-sample(ProportionalIntensities, 240000, replace = T, prob = NULL)
Intensities$SupFin<-as.numeric(Intensities$SupFin)
TotalFireSize<-sum(Intensities$SupFin)
Intensities$Weight<-Intensities$SupFin/TotalFireSize
IntensitiesWeighted<-sample(Intensities$INT,size=240000, replace=T, prob=Intensities$Weight)
IntensitiesWeighted <- as.numeric(IntensitiesWeighted)
length(IntensitiesWeighted)##
mean(IntensitiesWeighted)

##BioSimdata
load_data <- function(path) {
  files1 <- dir("C:\\Users\\psladmin\\Documents\\ModelOutput3ClimateChange\\WeatherBioSim\\Historical", pattern = '\\.csv', full.names = TRUE)
  tables <- lapply(files1, read.csv)
  do.call(cbind, tables)
}

Historical <- load_data(files1)
Years<-"Historical"
Historical$Year<-Years
names(Historical)
str(Historical)
Historical[c(4,5,7,8,10,11,13,14,16,17)] <- list(NULL)
colnames(Historical) <- c("Latitude","Longitude","TMean","PreciSept","PreciMayAug","TMaxFeb","TMinAug","TMinJuAu","Year")
head(Historical)
length(Historical[,1]) ##146460
Historical$x <-rep(1:30,each=4882)
Historical$y<-rep(1:30, 4882)
Historical$Index<-paste(Historical$x,Historical$y, sep = "_")
Historical$Concatanate <- paste(Historical$Latitude,Historical$Longitude, sep = "_")
length(Historical[,1])

#
#load_data <- function(path) {
#  files2 <- dir("C:\\Users\\psladmin\\Documents\\ModelOutput3ClimateChange\\WeatherBioSim\\First", pattern = '\\.csv', full.names = TRUE)
#  tables <- lapply(files2, read.csv)
#  do.call(cbind, tables)
#}
#
#First <- load_data(files2)
#Years<-"First"
#First$Year<-Years
#names(First)
#str(First)
#First[c(4,5,7,8,10,11,13,14,16,17)] <- list(NULL)
#colnames(First) <- c("Latitude","Longitude","TMean","PreciSept","PreciMayAug","TMaxFeb","TMinAug","TMinJuAu","Year")
#head(First)
#length(First[,1])#146460
#First$x <-rep(1:30,each=4882)
#First$y<-rep(1:30, 4882)
#First$Index<-paste(First$x,First$y, sep = "_")
#
#load_data <- function(path) {
#  files3 <- dir("C:\\Users\\psladmin\\Documents\\ModelOutput3ClimateChange\\WeatherBioSim\\Second", pattern = '\\.csv', full.names = TRUE)
#  tables <- lapply(files3, read.csv)
#  do.call(cbind, tables)
#}
#
#Second <- load_data(files3)
#Years<-"Second"
#Second$Year<-Years
#names(Second)
#str(Second)
#Second[c(4,5,7,8,10,11,13,14,16,17)] <- list(NULL)
#colnames(Second) <- c("Latitude","Longitude","TMean","PreciSept","PreciMayAug","TMaxFeb","TMinAug","TMinJuAu","Year")
#head(Second)
#Second$x <-rep(1:30,each=4882)
#Second$y<-rep(1:30, 4882)
#Second$Index<-paste(Second$x,Second$y, sep = "_")
#
#load_data <- function(path) {
#  files <- dir("C:\\Users\\psladmin\\Documents\\ModelOutput3ClimateChange\\WeatherBioSim\\Third", pattern = '\\.csv', full.names = TRUE)
#  tables <- lapply(files, read.csv)
#  do.call(cbind, tables)
#}
#
#Third <- load_data(files)
#Years<-"Third"
#Third$Year<-Years
#names(Third)
#str(Third)
#Third[c(4,5,7,8,10,11,13,14,16,17)] <- list(NULL)
#colnames(Third) <- c("Latitude","Longitude","TMean","PreciSept","PreciMayAug","TMaxFeb","TMinAug","TMinJuAu","Year")
#head(Third)
#Third$x <-rep(1:30,each=4882)
#Third$y<-rep(1:30, 4882)
#Third$Index<-paste(Third$x,Third$y, sep = "_")
#
#Weather<-rbind(Historical, First, Second, Third)
#Weather$Concatanate <- paste(Weather$Latitude,Weather$Longitude, sep = "_")
#length(Weather[,1])



subpop <- function(path, d= Plots){
  res <- subset(d, Pathway %in% path ,select=c("Concatanate","ID_PEP_MES","Pathway","Latitude","Longitude"))
  return(res)
}


GetTrees <- function(T,Plots){
  res <- NULL
  for (i in Plots){
    res<-rbind(res,subset(T,ID_PEP_MES %in% i))
  }
  return(res)
}

GetClimate<-function(x,y){
  res<-NULL
  for (i in y){
    res<-rbind(res,subset(x,Concatanate %in% i)) 
  }
  return(res)
}


subintensity <- function(reg, Int, season){
  intense <- subset(Int, FIRE_CODE %in% reg&SprSummer%in%season,select=c("INT","IH","ISec","Season","Dates"))
  return(intense)
}


DiameterGrowthRate<-function(stand,dbhq,baq,growthindex){ #increase in cm per year
  gi <- growthindex
  m2ft2 <- 10.7639
  cminch <- 2.54
  mfoot <- 3.28
  dbhinches <- dbhq/cminch
  b1 <- 0.0008236 ##parameters taken from Teck and Hilt
  b2 <- 0.0549439
  st <- b1*(1-exp(-1*b2*dbhinches))
  pba<- (gi*mfoot)*st #potential basal area growth (foot2/year)
  ###convert to diameter growth rate
  b3 <- 0.012 ##modifier
  bad<-(stand*baq)* m2ft2
  BAL<-(rev(cumsum(rev(bad)))-bad)
  bagr <- as.vector(pba*(exp(-b3 *BAL)))   ##basal area growth rate
  diametergrowthrate <- ((0.00545415*dbhinches^2 + bagr)/0.00545415)^.5 - dbhinches ##diameter growth rate inches
  ddgr <- diametergrowthrate * cminch    ##transform to cm
  return (ddgr)##returns diameter growth rate in cm
}

}

##ARTEMIS
mortality <- function(g,dbh,ba){
  tmp1 <- g*ba
  cba <- rev(cumsum(rev(tmp1)))-tmp1 #cumulative basal area
  p <- 1-exp(-exp(-1.624+(-2.6181)+(0.1229*dbh)+(-0.8280*log(dbh))+
                    (0.0098*0*log(10))+(0.0208*(cba))+log(10)))
  pa <- p^1/10
  survan <- 1-pa
  return(survan)
}

###Recruitment after fire as a function of seed abscission schedule,
#granivory, seed mortality by fire and optimal seedbeds. Predicts number of 3 yr old seedlings (1.5 and 7 cm), unless seedbed#and drainage are optimal, then it may exceptionally reach 10 to 15 cm) (VanBoagaert et al. 2015)

SeedProd<-function(stand,baq){
  m2Ha <-1e4
  Bd <- sum(baq[5:15]*stand[5:15])/10000 #pre-fire basal area/area
  Qd <- 163400*Bd^0.95 #germinable seeds/m2 in the aerial seed bank
  Pq <- 1   #proportion of seed abscised
  Sas <- 0.58     # fraction of seed surviving pass fire
  m <- 0.0012 ##black spruce seed mass (g)
  w <- 0.14 #proportion of optimal seedbeds(Boifin and Munson 2013)
  Sj <- 0.43*(w*(1-exp(-1.83*m^0.43))+(1-w)*(1-exp(-0.33*m^0.76)))
  Fd <- Qd*Sj*Sas*Pq #the number of expected 3 year recuits/m2
  return(round(Fd*m2Ha)) ##per m2
}


#Peng
Height <- function(dbh){
  1.3+  1.065*(dbh^0.8868)
}


## Holdaway
crownratio <- function(dbh,BA){
  b1 <- 5.54
  b2 <- 0.0072
  b3 <- 4.2
  b4 <- 0.053
  y <- b1/(1+(b2*BA)) + (b3*(1-exp(-b4*dbh)))
  round((y-.45)/10, digits=4)
}


updateCR <- function(pcr,ccr,dbh,TopHeight,di,BA,DBA,stand){
  b1 <- 5.54
  b2 <- 0.0072
  b3 <- 4.2
  b4 <- 0.053
  if(sum(is.na(stand))>0)
    stop(message="updateCR stand NA")
  if(sum(is.na(di))>0)
    stop(message="updateCR di NA")
  if(sum(is.na(DBA))>0)
    stop(message="updateCR DBA NA")
  d1 <- (b3 * b4 * exp(-b4 * dbh))* di
  d2 <-(-b1 * b2/((1 + b2 * BA)^2))* DBA
  xH <- TopHeight+di #Heightgrowth
  dH <- xH-TopHeight  #deltaheight
  MaximumCR <-(ccr*TopHeight+dH)/xH   #maximumCR
  ccr <- ccr+ d1 + ifelse(ccr < pcr, 0, d2)
  DBA <- rep(DBA,15)
  newcr <- pmin(MaximumCR,ccr)
  ccr <- ifelse(DBA<0,newcr,ccr)
  ccr <- ifelse(stand>0,ccr,pcr) #deals with all/ANY zeros
  ccr
}



####CARBON DYNAMICS##
Stemwood <- function(DBH){
  stembiomass <- 0.0477*DBH^2.5147
  return(stembiomass)
}

Bark <- function(DBH){
  barkbiomass <- 0.0153*DBH^2.2429
  return(barkbiomass)
}

Branches <- function(DBH){
  branchbiomass <- 0.0278*DBH^2.0839
  return(branchbiomass)
}

Needles <- function(DBH){
  needlesbiomass <-0.1648*DBH^1.4143
  return(needlesbiomass)
}


####Coarse root equation is from Ouimet et al. 2008.
Coarse <- function(DBH){ ####only roots >5 mm
  rootbiomass <-0.0085*1.036*(DBH^2.87) ###all classes
  fineroot <- 0.00153*2.40*(DBH^1.123) #### 2mm-5mm
  coarsebiomass <-rootbiomass-fineroot
  return(coarsebiomass)
}
Fineroots <- function(dbh){  #### <5mm. kg/tree Chen et al. 2004
  finerootbiomass <-0.011*(dbh^1.9748)
  return(finerootbiomass)
}
#

SnagsCarbon<-function(mortality,BioMassCarbon){
  snags <- mortality
  Stemwoodbig <- BioMassCarbon[1,5:15]%*%(snags[5:15])
  Barkmerchantable <- BioMassCarbon[2,5:15]%*%(snags[5:15])
  Stemwoodsmall <- BioMassCarbon[1,1:4]%*%(snags[1:4])
  Barksmall <- BioMassCarbon[2,1:4]%*%(snags[1:4])
  Branches <- BioMassCarbon[3,1:15]%*%(snags[1:15])
  Foliage <- BioMassCarbon[4,1:15]%*%(snags[1:15])
  CRoot <- BioMassCarbon[5,1:15]%*%(snags[1:15])
  FRoot <- BioMassCarbon[6,1:15]%*%(snags[1:15])
  SnagC <- Stemwoodbig+Barkmerchantable
  SnagbranchC <- Branches+Stemwoodsmall+Barksmall
  SnagFoliage <- Foliage
  SnagCoarse <- CRoot
  SnagFine <- FRoot
  return(list(SnagC=SnagC,SnagbranchC=SnagbranchC,SnagFoliage=SnagFoliage,SnagCoarse=SnagCoarse,SnagFine=SnagFine))
}



UpdateCarbonPool3 <-function(CPool,ICB,GrowthCBiomass,biomass,Snags,Snagbranches,SnagFoliage,SnagCoarse,SnagFine,CarbonPoolTransferMatrix){
  
  NPP <- (GrowthCBiomass-ICB)##NPP
  tmp <- as.vector(t(CPool)%*%CarbonPoolTransferMatrix)
  ctmp <- as.vector(t(biomass)%*%Input_Matrix2)
  BiomassLostTurnover <- sum(ctmp[8:14])##only litterfall and root turnover
  Biomass <- sum(biomass)-BiomassLostTurnover
  ctmp[6] <-ctmp[6]+Snags ##adding C from snags
  ctmp[7] <-ctmp[7]+Snagbranches
  ctmp[9] <-ctmp[9]+SnagCoarse*(0.5)  ##adding C from coarse roots to AG fast
  ctmp[10] <-ctmp[10]+SnagFine*(0.5)+SnagFoliage
  ctmp[12] <-ctmp[12]+SnagFine*(0.5)
  ctmp[13] <-ctmp[13]+SnagCoarse*(0.5)##BG fast
  Inputs <-ctmp[6:14]  ##how much C is incorporated into the DOMCpools including snags (mortality)
  CPool <- tmp[2:10]+Inputs
  SoilCAtmFlux <- tmp[1]
  NEP <-(NPP-SoilCAtmFlux)
  return(list(NPP=NPP,NEP=NEP,Biomass=Biomass,BiomassLostTurnover=BiomassLostTurnover,CPool=CPool,
              SoilCAtmFlux=SoilCAtmFlux,Inputs=Inputs))
}




#MAT<-sample(MeanAnnualTemperature$température.moyenne.annuelle, 1, replace = FALSE, prob = NULL)
Decayrates <- function(MAT){
  REFT <-10 ##reference temperature
  Reduction <- MAT-REFT ##reduction
  Q10 <- c(2,2,2,2,2.65,2.65,2,2,1)
  BDR <- c(0.0187,0.072,0.034,0.1435,0.355,0.015,0.5,0.1435,0.0033) #base decay rates at 10C
  TEMPMOD <- exp(1)^((Reduction)*log(Q10)*0.1)
  STANDMOD <- 1
  ADR <- BDR*TEMPMOD*STANDMOD
  ADR
}

Input_Matrix2 <- matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                          0,(1-0.04),0,0,0,0,(0.04*0.25),0,(0.04*0.75),0,0,0,0,0,
                          0,0,(1-0.16),0,0,0,0,0,0,(0.16*1),0,0,0,0,
                          0,0,0,(1-0.02),0,0,0,0,(0.02*0.5),0,0,0,(0.02*0.5),0,
                          0,0,0,0,(1-0.64),0,0,0,0,(0.64*0.5),0,(0.64*0.5),0,0
),nrow=5, ncol=14, byrow=TRUE)
colnames(Input_Matrix2)<-c("Merchantable","Otherwood","Needles","Coarse","Fine","Snags","Snagbranch","Medium","AGfast","AGveryfast","AGslow","BGveryfast","BGfast","BGslow")
rownames(Input_Matrix2)<-c("Merchantable","Otherwood","Needles","Coarse","Fine")


###Fire module###

CanopyFuelStocksbs <- function(dbh){
  w <- 0.6329+0.02*dbh^2.2
  return(w)
}


FuelClass <- function(stand,dbh){  ## fuel per DBH class kg/ha  using Stocks allometry
  w <- CanopyFuelStocksbs(dbh)*stand
  return(w)
}

VerticalFuelProfile <- function(c,Top,Base){ #Base<Top must be quaranteed
  MaxTop <- max(Top)#maximum height of dbh class
  f <- numeric(MaxTop)#number of 1m height sections up to the maximum height
  for (i in seq(1:length(Top))){ #for each  dbh class we do the following:
    v <- seq(Base[i]+1,Top[i]) # 1 m sections per diameter class
    x <- 1 / length(v)  #here we approximate the shape as cylinder. each section gets the same proportion of biomass fuel
    f[v]<- f[v]+(c[i]*x) ## the vector f accumulates  for each 1m section the fuel of every dbh classes per hectare
  }
  f/10000 # since the units of f are  kg/ha*m we have to convert to kg/m3  that means convert hectares to square meters...  (fuel in a section of 1 ha in area and 1m in height)
}

## this gives me the critical crown base height for crowning
#based on the initial fire intensity (Van Wagner)
Zc <- function(I){
  C <- 0.010
  m <- 100  #moisture
  h <- 460 + 26*m
  criticalbase <- I^(2/3)/(C*h)
  return(criticalbase)
}

Crowning <- function(I,Top,Base,Bulk,b,DenCrit){
  zc <- ceiling(Zc(I))   #Critical base height at intensity (I)
  cb <- Base[Base<=zc]   # which dbh classes are burning
  ct <- Top[Base<=zc]     # flames extend to top of ct...
  CrownLayer <- 0        # initialise
  k <-length(cb)         # k= number of dbh classes that affected by fire
  if (k > 0){
    fc <- ct[k]              # Height of the last dbh class affected by surface fire. Returns the  flame length value that corresponds to the highest strata with fire. Maxflamelength
    ct <- Top[Base<=fc]      #  We evaluate again the CBH of dbh classes. Which ones have a CBH lower than the flame height. return the dbh classes with crown base height smaller than  max flame
    kk <- length(ct)
    if (kk>0){    ###We evaluate then if crowning can be sustained
      for (i in seq(kk,1,by=-1)){   #we are checking from top to down
        t <- max(ct[i],b)    ##maximum height of the last dbh class with crowning
        den <- sum(Bulk[b:t])/(t-b+1)   #running mean "Available CBD for combustion"
        if (den > DenCrit){
          CrownLayer <- i
          break
        }
      }
    }
  }
  CrownLayer
}


UpdateCrownLayer<-function(cl,Top,Base,Bulk,b,DenCrit){
  newcl<-0
  if (cl > 0) {
    for (i in seq(cl,1,by=-1)){
      t <- max(Top[i],b)
      den <- sum(Bulk[b:t])/(t-b+1)
      if (den>DenCrit){
        newcl<-i
        break
      }
    }
  }
  newcl
}


FlameLength<-function(I){
  0.0775*((I)^0.46)
}

UpdateIntensity<-function(I,Ht){
  h<-0
  h<-h+(Ht*0.5)	#add 1/2 mean canopy height for crown fires (Byram)
  259.83*(h^2.174) ##flame length intensity relationship
}

ScorchHeight<-function(I){
  0.1483*(I^0.667) #Van Wagner (1973)
}

CrownKill<-function(I,Top,CCR){
  z<-ScorchHeight(I)   ##sh
  ht<-Top
  cbh<-ht*(1-CCR)     #CBH
  cl<-ht-cbh
  tmp<-z-(ht-cl)
  tcls<-ifelse(z<cbh,0,tmp)
  cls<-ifelse(z>ht,cl,tcls)
  cvscy<-100*(cls/cl)
  return(cvscy)
}

ScorchMortality<-function(Bark,CK){
  bc<-6.316*(1-exp(-Bark))
  cc<- -0.000535*(CK^2) #
  1/(1+exp(-1.941+bc+ cc ))
}

#####Total Basal area lost
Basalost<-function (stand, baq, newstand){
  a<-sum(stand*baq)    ### Total basal area before
  b<-sum(newstand*baq)
  Percentage_T_Basalost<- 100*(1-(b/a))
  return(Percentage_T_Basalost)
}


###Historical climate and FRI

exe<-function(stand,Y,FRI){
  
  
  GI<-function(ClimateData.List){ ##climate sensitive growth index
    b0 <- 6.2
    b1 <- 0.2357
    b2 <- 0.01003
    b3 <- 0.3663
    b4 <- 1.426
    b5 <- -0.129
    b6 <- 0.00942
    b7 <- 0.00483
    growthIndex <- rep(0,length(ClimateData.List[,3])-1)
    for (i in seq(from=2, to=length(ClimateData.List[,3]))){
      gi <- b0 + b1*ClimateData.List[i,6] + b2*ClimateData.List[i,6]^2 + b3*ClimateData.List[i,7] + b4*ClimateData.List[i,8]+b5*ClimateData.List[i,8]^2+b6*ClimateData.List[i-1,4]+b7*ClimateData.List[i,5]
      growthIndex[i]<- gi
    }
    growthIndex
  }
  
  CorrectDupli<-function(x){
    newdata<-NULL
    newdata <- subset(x, Year=="Historical",select=Latitude:Concatanate)
    ClimateData<-subset(newdata[1:30,]) ##first 30 of the time series
    return(ClimateData)
  }
  
  
  Sampled <- subpop (path=c("Pathway17","Pathway19","Pathway20","Pathway21"), d=PoolPlots2)
  Sampled <- Sampled[sample(1:dim(Sampled)[1], size=1, replace=T),]
  PlotID<-Sampled[1,2]
  WeatherPlot<-Sampled[1,1]
  Latitude <- Sampled[1,4]
  Longitude <- Sampled[1,5]
  ClimateData<-GetClimate(Historical,WeatherPlot)
  ClimateData.List <- CorrectDupli(ClimateData) ##to fix duplicates 
  GrowthIndex <- GI(ClimateData.List)  ##growthindex 
  GrowthIndex[1]<- 14  
  
  Period <- ClimateData.List[,9] 
  Tree.List <- GetTrees (Tree,Plots=Sampled)
  as.character(Tree.List$DBH)
  as.factor(Tree.List$ESSENCE)
  Tree.List <- as.data.frame(lapply(Tree.List[,],function(x)rep(x,25)))
  range.DBH<-c(seq(1,30, by=2), 100)
  #Resume the results by class
  Tree.List$DBH <- as.numeric(Tree.List$DBH)
  stand <- table(cut(Tree.List$DBH, breaks=range.DBH, labels=seq(1,15)))
  stand[1:4] <- stand[1:4]*10
  ###Partition the basal area of big trees >31 cm and add number of trees that the surplus of basal area represents
  basal_big_class <- 0.0707905544
  BAB <- rep(0,100)
  TBA <- 3.142*(Tree.List[Tree.List[,3]>31,3]/200)^2
  BAB <- round(TBA/basal_big_class,digits=0)
  y <- sum(BAB)
  stand[15] <- stand[15]+y
  n <- length(stand)
  N1s <- rep(1,n)
  N0s <- rep(0,n)
  BurnMortP <- numeric(n)
  BurnMortPpar <- numeric(n)
  shannon <- 0
  ScH <- 0
  cl <- 0
  ###################################
  ############################################################
  dbhl <-c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29) ###dbh lower limit
  dbhu <- c(3,5,7,9,11,13,15,17,19,21,23,25,27,29,31) ##dbh upper limit
  dbhq <-sqrt((dbhu^3-dbhl^3)/((dbhu-dbhl)*3)) ###assuming a uniform dbh distribution over interval#and the basal area at the stage mean dbh is...
  baq <-(dbhq/2)^2*pi/1e4
  ######INITIALISE VARIABLES (TOP,BASE,CROWN RATIOS,INTENSITY,BARK)
  Top <-ceiling(Height(dbhq))##ceiling gives the upper integer of the height at diameter class i
  Top <-ifelse(Top<2,rep(2,length(Top)),Top)	#Top>=2m
  HeightUpdated <-Top
  PCR <-crownratio(dbhq,sum(baq*stand))
  CCR <-PCR
  Base<-floor(HeightUpdated*(1-CCR)) #####gives crown base height
  Bark <-0.032*2.54*dbhq ##inches to cm!! Black spruce Equation taken from Behave plus fire modelling system  (thickness in cm)
  hatom2 <-1e4
  b <-3
  DenCrit <-0.11
  iw <-dbhu-dbhl
  BarkThickness <-0.032*2.54*dbhq
  n <-length(stand)
  N1s <-rep(1,n)
  N0s <-rep(0,n)
  n <- length(stand)
  m2Ha <- 1e4
  biocar_factor <- 0.5
  MgKg <- 1000
  SaplingSurvival <- 0.98 ###sapling survival of different initial sizes (Matthias et al. 2003)
  pBurn <-1/FRI
  NF <-0
  fire.year <- NULL
  RegenLagTime <- 27   ## Assume it takes a natural seedling 30 years to grow to 214 cm (DBH 1.0 cm), assume a 27 year delay for a 3-yr old natural black spruce seedling( 20 years was rather optimistic; VanBoagaert et al. 2015)
  RegenRegular <- 60 ##calibration exercise 11/11/2014
  RegenIrregular <- 55 ##calibration exercise 11/11/2014
  shannon <- diversity(stand[5:15], index="shannon", MARGIN=1, base=exp(1)) ##updated shannon
  RegenCohorts <- rpois(RegenLagTime, ifelse(shannon<1.7,RegenRegular,RegenIrregular))
  carbon <-c(5.8,2.1,24,9.89,9.0,36,1.69,3.75,74) * MgKg
  SoilCAtmFlux <- 3 * MgKg
  Snags <- 0
  Snagbranches <- 0
  SnagFoliage <- 0
  SnagCoarse <- 0
  SnagFine <- 0
  ICB <- 3.11 * MgKg
  NEP <- 0
  ICB <- 52* MgKg
  
  #This is to get the simulated plots stand structure characteristics
  ###Initialize object variables to save simulation results####
  
  Size<-matrix(0,Y,n,byrow=TRUE)
  PreFireStand<-matrix(0,Y,n,byrow=TRUE)
  Recruits<-numeric(Y)
  BA<-numeric(Y)
  FireDeaths<-matrix(0,Y,n,byrow=TRUE)
  Senescence<-matrix(0,Y,n,byrow=TRUE)
  Transition<-matrix(0,Y,n,byrow=TRUE) ##transition probability
  Crecimiento<-matrix(0,Y,n,byrow=TRUE)
  Parcela<-matrix(0,Y,n,byrow=TRUE)
  Muertos<-matrix(0,Y,n,byrow=TRUE)
  Mortality<-matrix(0,Y,n,byrow=TRUE)##probability of mortality
  Structure<-numeric(Y)
  InitialIntensity<-numeric(Y)
  FireSeason<-numeric(Y)
  BALost<-numeric(Y)
  Delta_BA<-numeric(Y)
  DeltaN<-matrix(0,Y,n,byrow=TRUE)
  CR<-matrix(0,Y,n,byrow=TRUE)
  Heights<-matrix(0,Y,n,byrow=TRUE)
  Heights<-matrix(0,Y,n,byrow=TRUE)
  ShiftCrownratio<-matrix(0,Y,n,byrow=TRUE)
  ShiftHeights<-matrix(0,Y,n,byrow=TRUE)
  DiameterGrowth<-matrix(0,Y,n,byrow=TRUE)
  SnagCProduction<-numeric(Y)
  Turnover<-numeric(Y) ##turnover
  DOMC_Pool<-matrix(0,Y,length(carbon),byrow=TRUE)
  DOM_Flux<-numeric(Y)
  DOM_Inputs<-matrix(0,Y,9,byrow=TRUE)
  BioMass<-matrix(c(Stemwood(dbhq),Bark(dbhq),Branches(dbhq),Needles(dbhq),Coarse(dbhq),Fineroots(dbhq)),nrow=6,
                  ncol=length(dbhq),byrow=TRUE)
  BioMassCarbon<-BioMass*biocar_factor ##Biomass C per diameter class
  Merchantable<-c(0,0,0,0,Stemwood(dbhq[5:15])+Bark(dbhq[5:15]))
  OtherWood<-c(Stemwood(dbhq[1:4])+Bark(dbhq[1:4]),rep(0,11))+ Branches(dbhq)
  InitialCBiomass<-sum(BioMassCarbon%*%as.matrix(stand)) ##bio1 for the site
  ICB<- InitialCBiomass
  CarbonBiomass1<-numeric(Y)#after mortality and turnover
  NetPrimaryProductivity<-numeric(Y)
  CarbonBiomass2<-numeric(Y)##after growth and recruitment
  TotalLiveBiomass<-numeric(Y)
  ACB<- matrix(0,Y,3,byrow=TRUE)
  BCB<- matrix(0,Y,2,byrow=TRUE)
  AppliedDecayRates<-matrix(0,Y,9,byrow=TRUE)
  EmpiricalTemperature <- matrix(0,Y,120, byrow=TRUE)
  NetEcosystemProduction<-numeric(Y)
  Rh<-numeric(Y)
  EmpiricalTemperature <- numeric(Y)
  PlotIDs <-numeric(Y)
  LPeriod <-as.numeric(Y)
  Periods <-matrix(0,120,Y, byrow=TRUE)  
  CCGrowthIndex <- numeric(Y)
  MATs<-as.numeric(ClimateData.List[,3])  ##mean annual temperatures 1981-2100 
  for (y in 1:Y){ #Things that have to be reinitialized and updated
    MAT<-sample(MATs, size=1, replace=T) 
    GI<-sample(GrowthIndex, size=1, replace=T)
    LPeriod[y]<-length(ClimateData.List[,9])
    Periods[,y] <- Period
    EmpiricalTemperature[y] <- MAT
    PlotIDs[y]<- WeatherPlot      
    HeadIntensity <- sample(IntensitiesWeighted, size=1, replace=F)
    I<-HeadIntensity
    Istart <-I
    shannon<-diversity(stand[5:15], index="shannon", MARGIN=1, base=exp(1)) ##updated shannon
    Structure[y]<-shannon
    Firedeaths<-rep(0,15)
    BALost[y]<-0
    bay<-sum(stand*baq)
    AppDecayRates<-Decayrates(MAT)
    AppliedDecayRates[y,]<-AppDecayRates
    Decayrate<-rep(0,9)
    Decayrate[1]<-AppDecayRates[1]
    Decayrate[2]<-AppDecayRates[2]
    Decayrate[3]<-AppDecayRates[3]
    Decayrate[4]<-AppDecayRates[4]
    Decayrate[5]<-AppDecayRates[5]
    Decayrate[6]<-AppDecayRates[6]
    Decayrate[7]<-AppDecayRates[7]
    Decayrate[8]<-AppDecayRates[8]
    Decayrate[9]<-AppDecayRates[9]
    
    CarbonPoolTransferMatrix<-matrix(c(
      Decayrate[1]*0.83,(1-Decayrate[1]-0.08),0,0.08,0,0,Decayrate[1]*(1-0.83),0,0,0,
      Decayrate[2]*0.83,0,(1-Decayrate[2]-0.10),0,0.10,0,Decayrate[2]*(1-0.83),0,0,0,
      Decayrate[3]*0.83,0,0,1-Decayrate[3],0,0,Decayrate[3]*(1-0.83),0,0,0,
      Decayrate[4]*0.83,0,0,0,(1-Decayrate[4]),0,Decayrate[4]*(1-0.83),0,0,0,
      Decayrate[5]*0.815,0,0,0,0,(1-Decayrate[5]),Decayrate[5]*(1-0.815),0,0,0,
      Decayrate[6]*1,0,0,0,0,0,(1-Decayrate[6]-0.006),0,0,0.006,
      Decayrate[7]*0.83,0,0,0,0,0,0,(1-Decayrate[7]),0,Decayrate[7]*(1-0.83),
      Decayrate[8]*0.83,0,0,0,0,0,0,0,(1-Decayrate[8]),Decayrate[8]*(1-0.83),
      Decayrate[9]*1,0,0,0,0,0,0,0,0,1-Decayrate[9]
    ),nrow=9,ncol=10,byrow=TRUE)
    colnames(CarbonPoolTransferMatrix)<-c("Atm","Snags","Snagbranch","Medium","AGfast","AGveryfast","AGslow","BGveryfast","BGfast","BGslow")
    rownames(CarbonPoolTransferMatrix)<-c("Snags","Snagbranch","Medium","AGfast","AGveryfast","AGslow","BGveryfast",     "BGfast","BGslow")
    
    ##FIRE MODULE
    #Evaluates if a fire arrives its intensity, crowning, and post-fire tree mortality and regeneration
    
    if(runif(1) < pBurn){ ## Determine if a fire happens
      fire.year <- c(fire.year,y)
      NF <- NF+1 ## Update number of fires that occurred during the simulation
      RegenCohorts<-rep(0,RegenLagTime) ## KILL ALL regenerating trees
      prefirestand<-stand
      PreFireStand[y,]<-prefirestand
      NewRegen<-SeedProd(prefirestand,baq)
      Fuel<-as.numeric(FuelClass(stand,dbhq)) #Kg/ha per dbh class
      VF<-VerticalFuelProfile(Fuel,HeightUpdated,Base)
      cl<-Crowning(I,HeightUpdated,Base,VF,b,DenCrit)
      #if there is crowning, then update the crown layer affected by fire
      cl<-ifelse (cl>0,UpdateCrownLayer(cl,HeightUpdated,Base,VF,b,DenCrit),
                  Crowning(I,HeightUpdated,Base,VF,b,DenCrit))
      u<-UpdateIntensity(I,HeightUpdated[cl])
      I<-ifelse(cl>0,max(UpdateIntensity(I,HeightUpdated[cl])),Istart)
      ScH<-ScorchHeight(I)
      CK<-CrownKill(I,HeightUpdated,CCR)
      BurnMortP<-ScorchMortality(BarkThickness,CrownKill(I,HeightUpdated,CCR))
      Firedeaths<-rbinom(N1s,stand,BurnMortP)
      Firedeaths<-ifelse(is.na(Firedeaths),N0s,Firedeaths)
      FireDeaths[y,]<-Firedeaths
      newstand<-stand-Firedeaths
      severity<-Basalost(stand,baq,newstand)
      stand<-newstand ###update stand after a fire
      InitialIntensity[y]<-I #adjusted intensity using Catchpole et.al 1992
      #FireSeason[y]<-Season
      BALost[y]<-severity
      tmp3<-carbon[3]*0.392140  ##carbon consumed in the medium carbon pool
      tmp4<-carbon[4]*0.6415 ##carbon consumed in the Ag fast carbon pool
      tmp5<-carbon[5]*0.968533 ##carbon consumed in the Ag very fast pool
      tmp6<-carbon[6]*0.09001 ##carbon consumed in the Ag slowpool
      carbon[3]<-carbon[3]-tmp3
      carbon[4]<-carbon[4]-tmp4
      carbon[5]<-carbon[5]-tmp5
      carbon[6]<-carbon[6]-tmp6
    }
    else {
      NewRegen <- rpois(1,ifelse(shannon<1.7,RegenRegular,RegenIrregular))
      
    }
    
    ##extract this year
    Recruits[y]<-RegenCohorts[RegenLagTime]
    RegenCohorts<-c(NewRegen,RegenCohorts[1:RegenLagTime-1])
    stand[1] <- stand[1] + Recruits[y]
    CCRRecruits<-Recruits[y]*PCR[1]
    HeightRecruits<-Recruits[y]*Top[1]
    Parcela[y,]<-stand ##stand after regeneration..should capture pulses
    
    #GROWTH
    CCGrowthIndex[y] <- GI
    annual_diam_incre <- DiameterGrowthRate(stand,dbhq,baq,GI)
    adi<-annual_diam_incre
    graduating<- 1/(iw/adi) #Transition probabilities
    growth<- rbinom(N1s, stand, graduating) #stohastically
    stand<-stand-growth
    stand<-stand+c(0,growth[1:n-1]) #after growth
    BAIncrement<-sum(growth*baq)    #patch BA increment due to growth
    CCRgrowth<-growth*CCR           #the trees that grow bring their crown ratio
    Heightgrowth<-growth*HeightUpdated
    GrowthCBiomass<- sum(BioMassCarbon%*%as.matrix(stand)) ##calculate biomass due to growth..Biomass that has not been lost to turnover or mortality
    
    
    ##MORTALITY
    surviving<-mortality(stand,dbhq,baq)
    surviving[1:4]<-SaplingSurvival
    Senescencedeaths<-rbinom(N1s,stand,1-surviving)
    deaths<-Senescencedeaths+Firedeaths
    Muertos[y,]<-deaths
    Senescence[y,]<-Senescencedeaths
    SnagCpools<-SnagsCarbon(deaths,BioMassCarbon)
    Snags<-SnagCpools$SnagC
    Snagbranches<-SnagCpools$SnagbranchC
    SnagFoliage<-SnagCpools$SnagFoliage
    SnagCoarse<-SnagCpools$SnagCoarse
    SnagFine<-SnagCpools$SnagFine
    stand<-stand-Senescencedeaths  ###fire deaths are taken care in the fire module
    DeltaBA<- BAIncrement-sum(deaths*baq) ##net basal area increment after mortality (fire and senescence) and growth
    DeltaStand<-sum(stand)  ##net change in density after mortality (fire and senescence) and growth
    if(sum(stand<0)>0)
      stop(message="neg count 3")
    Delta_BA[y]<-DeltaBA
    DeltaN[y,]<-deaths+growth
    
    
    ##Dynamically updating crown ratios and heights of natural (all,regenerated and fire derived recruits)
    
    xH <- Top+adi ##Heightgrowth
    dH <- xH-Top  ##deltaheight
    MaximumCR <- (CCR*Top+dH)/xH
    TotalN <- stand+c(0,growth[1:n-1])
    TotalN[1]<-TotalN[1]+Recruits[y]
    CCR <- updateCR(PCR,CCR,dbhq,Top,adi,bay,DeltaBA,stand)##CCR after recruitment,growth,mortality
    CCRnow <- stand*CCR
    ShiftCR <- CCRnow+c(0,CCRgrowth[1:n-1])
    ShiftCR[1] <- ShiftCR[1]+CCRRecruits
    CCR <- ifelse(stand>0,pmin(MaximumCR,ShiftCR/TotalN),PCR)
    ##Updating heights
    Top <- Height(dbhq)
    HeightUpdated <- Top+adi ##Heightgrowth
    Heightnow<-stand*HeightUpdated
    ShiftHeight<-Heightnow+c(0,Heightgrowth[1:n-1])
    ShiftHeight[1]<-ShiftHeight[1]+HeightRecruits
    HeightUpdated<-ifelse(stand>0,ShiftHeight/TotalN,Top)
    Base <- HeightUpdated*(1-CCR) #####gives crown base height
    ##Carbon Dynamics
    Stemwoodsmall<-BioMassCarbon[1,1:4]%*%(stand[1:4])
    Barkmerchantable<-BioMassCarbon[2,5:15]%*%(stand[5:15])
    LiveBiomassCPools<-BioMassCarbon%*%(stand) #biomass C kg/ha
    ##Correction made October 2014
    LiveBiomassCPoolsCorrected<-matrix(0,nrow=5,ncol=1)
    LiveBiomassCPoolsCorrected[1,1]<- LiveBiomassCPools[1,1]-Stemwoodsmall+Barkmerchantable# Merchantable+Bark
    LiveBiomassCPoolsCorrected[2,1]<- LiveBiomassCPools[2,1]-Barkmerchantable+LiveBiomassCPools[3,1]+
      Stemwoodsmall   #Otherwood+Bark
    LiveBiomassCPoolsCorrected[3,1]<-LiveBiomassCPools[4,1]#Foliage
    LiveBiomassCPoolsCorrected[4,1]<-LiveBiomassCPools[5,1]#Coarse
    LiveBiomassCPoolsCorrected[5,1]<-LiveBiomassCPools[6,1]#Fine
    
    #Corrected October to match IPCC Good Practice Guidance
    ACB[y,]<-LiveBiomassCPoolsCorrected[1:3]
    BCB[y,]<-LiveBiomassCPoolsCorrected[4:5]
    TotalLiveBiomass[y]<-sum(LiveBiomassCPoolsCorrected)
    DOMCPools<-UpdateCarbonPool3(carbon,ICB,GrowthCBiomass,LiveBiomassCPoolsCorrected,Snags,Snagbranches,SnagFoliage,SnagCoarse,SnagFine,CarbonPoolTransferMatrix)
    ICB<-DOMCPools$Biomass
    CarbonBiomass1[y]<- GrowthCBiomass  #Net biomass gain
    NetPrimaryProductivity[y]<-DOMCPools$NPP
    CarbonBiomass2[y]<-ICB
    carbon<-DOMCPools$CPool
    Inputs<-DOMCPools$Inputs #includes carbon turnover+carbon from snags
    turnover<-DOMCPools$BiomassLostTurnover #Different from Inputs. only contains C from litterfall and root turnover. use this when calculating NPP
    DOMC_Pool[y,]<-carbon #Nine DOM carbon pools (snag=snagbracnhes etc.)
    DOM_Flux[y]<-DOMCPools$SoilCAtmFlux
    DOM_Inputs[y,]<-Inputs #includes carbon turnover+carbon from snags
    Turnover[y]<-turnover #how much biomassC is lost due only to litterfall
    NetEcosystemProduction[y]<-DOMCPools$NEP
    Rh[y]<-DOMCPools$SoilCAtmFlux
    BA[y]<-bay
    CR[y,]<-CCR
    Heights[y,]<-HeightUpdated
    ShiftCrownratio[y,]<-ShiftCR
    ShiftHeights[y,]<-ShiftHeight
    DiameterGrowth[y,]<-adi ##annual diamater increment
    Transition[y,]<-graduating ##transition probabilities
    Mortality[y,]<-1-surviving ###probability of mortality per DBH class
    Crecimiento[y,]<-growth
    Muertos[y,]<-deaths
    Size[y,]<-stand ###  after all processes..FINAL size..or initial one
    
  }
  
  res<-list(Parcela=Parcela,Size=Size,DeltaN=DeltaN,BA=BA,PreFireStand=PreFireStand,Senescence=Senescence,AppliedDecayRates=AppliedDecayRates,
            EmpiricalTemperature=EmpiricalTemperature,FireDeaths=FireDeaths,Muertos=Muertos,Structure=Structure,Crecimiento=Crecimiento
            ,Mortality=Mortality,Transition=Transition,Structure=Structure,DiameterGrowth=DiameterGrowth,CR=CR,
            ShiftCrownratio=ShiftCrownratio,ShiftHeights=ShiftHeights,Heights=Heights,BALost=BALost, NF=NF,InitialIntensity=InitialIntensity,fire.year=fire.year,Recruits=Recruits,ACB=ACB,BCB=BCB,TotalLiveBiomass=TotalLiveBiomass,CarbonBiomass1=CarbonBiomass1,CarbonBiomass2=CarbonBiomass2,Turnover=Turnover,SnagCProduction=SnagCProduction,DOMC_Pool=DOMC_Pool,DOM_Flux=DOM_Flux,DOM_Inputs=DOM_Inputs,NetPrimaryProductivity=NetPrimaryProductivity, Rh=Rh,NetEcosystemProduction=NetEcosystemProduction,LPeriod=LPeriod,PlotIDs=PlotIDs,Periods=Periods,PlotIDs=PlotIDs,CCGrowthIndex =CCGrowthIndex)
  
  return(res)
}


abc<-exe(stand,100,916)


n.iter<-1000#plots to check
Y<-2400
FRI<-916
#stand dynamics
Ba_s<- matrix(0,n.iter,Y,byrow=T)
Size_list <- vector("list", n.iter) # create list
PreFireStand_s<-vector("list", n.iter)
Transition_list <- vector("list", n.iter)
Structure_s<-matrix(0,n.iter,Y,byrow=T)
Recruits_s<-matrix(0,n.iter,Y,byrow=T)
Heights_list<-vector("list", n.iter)
CR_list<-vector("list", n.iter)
GrowthIndex_list <- vector ("list", n.iter)
DiameterGrowth_list<- vector("list", n.iter)
DecayRates_list<-vector("list", n.iter)
#carbon dynamics
CarbonBiomass1_s<- matrix(0,n.iter,Y,byrow=T)
CarbonBiomass2_s<- matrix(0,n.iter,Y,byrow=T)
DOM_Pool_list<- vector("list", n.iter)
DOM_Inputs_list<- vector("list", n.iter)
DOM_Flux_s<-matrix(0,n.iter,Y,byrow=T)
PrimaryProductivity_s<-matrix(0,n.iter,Y,byrow=T)
NetEcosystemProduction_s<-matrix(0,n.iter,Y,byrow=T)
Rh_s<-matrix(0,n.iter,Y,byrow=T)
Turnover_s<-matrix(0,n.iter,Y,byrow=T)
TotalLiveBiomass_s<-matrix(0,n.iter,Y,byrow=T)
AbovegroundLiveBiomass_list<-vector("list",n.iter)
BelowgroundLiveBiomass_list<-vector("list",n.iter)
###Fire disturbance
Severity_s<-matrix(0,n.iter,Y,byrow=T)
InitialIntensity_s<-matrix(0,n.iter,Y,byrow=T)
NF_s<-numeric(n.iter)

for(i in 1:n.iter){
  Size_list[[i]] <- matrix(0,Y,15,byrow=T)
  DiameterGrowth_list[[i]]<- matrix(0,Y,15,byrow=T)
  PreFireStand_s[[i]]<-matrix(0,Y,15,byrow=T)
  Transition_list[[i]]<-matrix(0,Y,15,byrow=T)
  Heights_list[[i]]<-matrix(0,Y,15,byrow=T)
  CR_list[[i]]<-matrix(0,Y,15,byrow=T)
  AbovegroundLiveBiomass_list[[i]]<-matrix(0,Y,3,byrow=T)
  BelowgroundLiveBiomass_list[[i]]<-matrix(0,Y,2,byrow=T)
  DOM_Pool_list[[i]]<- matrix(0,Y,9,byrow=T)
  DOM_Inputs_list[[i]]<- matrix(0,Y,9,byrow=T)
  DecayRates_list[[i]]<- matrix(0,Y,9,byrow=T)
}

do.call(rbind, Size_list)
do.call(rbind, Transition_list)
do.call(rbind, DiameterGrowth_list)
do.call(rbind, Heights_list)
do.call(rbind, CR_list)
do.call(rbind, AbovegroundLiveBiomass_list)
do.call(rbind, BelowgroundLiveBiomass_list)
do.call(rbind, DOM_Pool_list)
do.call(rbind, DOM_Inputs_list)
do.call(rbind, DecayRates_list)

for (i in 1:n.iter){
  CarbonModel<-exe(stand,Y,FRI)
  Ba_s[i,]<-CarbonModel$BA
  Structure_s[i,]<-CarbonModel$Structure
  Recruits_s[i,]<-CarbonModel$Recruits
  Size_list[[i]]<-CarbonModel$Size
  DiameterGrowth_list[[i]] <-CarbonModel$DiameterGrowth
  PreFireStand_s[[i]] <-CarbonModel$PreFireStand
  Transition_list[[i]]<-CarbonModel$Transition
  Heights_list[[i]]<-CarbonModel$Heights
  CR_list[[i]] <-CarbonModel$CR
  TotalLiveBiomass_s[i,]<-CarbonModel$TotalLiveBiomass
  AbovegroundLiveBiomass_list[[i]]<-CarbonModel$ACB
  BelowgroundLiveBiomass_list[[i]]<-CarbonModel$BCB
  EmpiricalTemperature_list[[i]] <-CarbonModel$EmpiricalTemperature
  GrowthIndex_list [[i]] <- CarbonModel$CCGrowthIndex
  CarbonBiomass1_s[i,] <-CarbonModel$CarbonBiomass1
  CarbonBiomass2_s[i,]<-CarbonModel$CarbonBiomass2
  PrimaryProductivity_s[i,]<-CarbonModel$NetPrimaryProductivity
  NetEcosystemProduction_s[i,]<-CarbonModel$NetEcosystemProduction
  Rh_s[i,]<-CarbonModel$Rh
  Turnover_s[i,]<-CarbonModel$Turnover
  DOM_Pool_list[[i]]<-CarbonModel$DOMC_Pool
  DecayRates_list[[i]]<-CarbonModel$AppliedDecayRates
  DOM_Flux_s[i,] <-CarbonModel$DOM_Flux
  DOM_Inputs_list[[i]] <-CarbonModel$DOM_Inputs
  Severity_s[i,]<-CarbonModel$BALost
  InitialIntensity_s[i,]<-CarbonModel$InitialIntensity
  NF_s[i]<-CarbonModel$NF
  print(i)
}




Year<-2400
plots<-1000
BasalArea<-Ba_s[,Year]
#Recruitment<-Recruits_s[,Year]
BiomassTurnover<-Turnover_s[1:plots,Year]
ICB<-CarbonBiomass2_s[1:plots,Year]
Inputs<-sapply(DOM_Inputs_list, rowSums)
Litterfall<-Inputs[Year,]
MerchantableStemwood<-sapply(AbovegroundLiveBiomass_list, function(m) m[Year,1])
Otherwood<-sapply(AbovegroundLiveBiomass_list, function(m) m[Year,2])
Foliage<-sapply(AbovegroundLiveBiomass_list, function(m) m[Year,3])
CoarseRoots<-sapply(BelowgroundLiveBiomass_list, function(m) m[Year,1])
FineRoots<-sapply(BelowgroundLiveBiomass_list, function(m) m[Year,2])
AGTotal<-MerchantableStemwood+Otherwood+Foliage #this should be the same as in sapplyTotallivem[500,1]
BGTotal<-CoarseRoots+FineRoots #@yr500
BiomassLiveCStock<-AGTotal+BGTotal
BiomassLiveCStock1<-TotalLiveBiomass_s[1:plots,Year] #1:1000 plots @ yr500
NPP<-PrimaryProductivity_s[1:plots,Year]
SoilRespiration<-Rh_s[1:plots,Year]
NEP<-NetEcosystemProduction_s[1:plots,Year]
StandStructure<-Structure_s[1:plots,Year]
#Intensities<-InitialIntensity_s[1:plots,1:Year]
Snags<-sapply(DOM_Pool_list, function(m) m[Year,1]) #
SnagBranch<-sapply(DOM_Pool_list, function(m) m[Year,2])
AGMedium<-sapply(DOM_Pool_list, function(m) m[Year,3])
AGfast<-sapply(DOM_Pool_list, function(m) m[Year,4])
AGveryfast<-sapply(DOM_Pool_list, function(m) m[Year,5])
AGslow<-sapply(DOM_Pool_list, function(m) m[Year,6])
BGveryfast<-sapply(DOM_Pool_list, function(m) m[Year,7])
BGfast<-sapply(DOM_Pool_list, function(m) m[Year,8])
BGslow<-sapply(DOM_Pool_list, function(m) m[Year,9])
SoilCStock<-Snags+SnagBranch+AGMedium+AGfast+AGveryfast+AGslow+BGveryfast+BGfast+BGslow
EcosystemCStock<-BiomassLiveCStock1+SoilCStock
FireReturnInterval<-rep(FRI,plots)





Initial152NoClimate <- cbind(FireReturnInterval,StandStructure,ICB,BasalArea,MerchantableStemwood,Otherwood,Foliage,CoarseRoots,FineRoots,AGTotal,BGTotal,BiomassTurnover, Litterfall,NPP,SoilRespiration,NEP,Snags,SnagBranch,AGMedium,AGfast,AGveryfast,AGslow,BGveryfast,BGfast,BGslow,BiomassLiveCStock,SoilCStock,EcosystemCStock)
save("Initial152NoClimate", file="C:\\Users\\psladmin\\Documents\\ModelOutput3ClimateChange\\InitialCstocks\\Initial152NoClimate.RData")
