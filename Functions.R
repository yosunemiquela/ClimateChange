########################
## Sampling functions ##
########################
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

GetClimate<-function(ClimateData,SamplePlotID){
  res <- NULL
  for (i in SamplePlotID){
    res <- rbind(res,subset(ClimateData,Concatanate %in% i))
  }
  return(res)
}

GetDC<-function(DroughtCodeData,SamplePlotID){
  res <- NULL
  for (i in SamplePlotID){
    res <- rbind(res,subset(DroughtCodeData,Concatanate %in% i))
  }
  return(res)
}

subintensity <- function(reg, Int, season){
  intense <- subset(Int, FIRE_CODE %in% reg&SprSummer%in%season,select=c("INT","IH","ISec","Season","Dates"))
  return(intense)
}

CorrectDupli<-function(x){
  newdata<-NULL
  newdata2<-NULL
  newdata3<-NULL
  newdata4<-NULL 
  newdata <- subset(x, Year=="Historical",select=Latitude:Concatanate)
  newdata2 <- subset(x, Year=="First",select=Latitude:Concatanate)
  newdata3 <- subset(x, Year=="Second",select=Latitude:Concatanate)
  newdata4 <- subset(x, Year=="Third",select=Latitude:Concatanate)  
  Histnewdata<-subset(newdata[1:30,]) ##first 30 of the time series
  Firstnewdata<-subset(newdata2[1:30,]) ##first 30 of the time series
  Secondnewdata<-subset(newdata3[1:30,]) ##first 30 of the time series
  Thirdnewdata<-subset(newdata4[1:30,]) ##first 30 of the time series
  ClimateData.List<-rbind(Histnewdata,Firstnewdata,Secondnewdata,Thirdnewdata)
  return(ClimateData.List)
}

###########################
## Stand dynamics module ##
###########################

GrowthIndexfunction <- function(ClimateData.List){ ##climate sensitive growth index
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


DiameterGrowthRate<-function(stand,dbhq,baq,growthindex){ #increase in basal area per year
  gi <- growthindex
  m2ft2 <- 10.7639
  cminch <- 2.54
  mfoot <- 3.28
  dbhinches <- dbhq/cminch
  b1 <- 0.0008236 ##parameters taken from Teck and Hilt
  b2 <- 0.0549439
  st <- b1*(1-exp(-1*b2*dbhinches))
  pba <- (gi*mfoot)*st #potential basal area growth (foot2/year)
  #pbasqcm <- pba * 929.03 ##potential basal area growth (cm2) ##to compare with Teck and Hilt 1991
  ###convert to diameter growth rate
  b3 <- 0.012 ##modifier
  bad<-(stand*baq)* m2ft2
  BAL<-(rev(cumsum(rev(bad)))-bad)
  bagr <- as.vector(pba*(exp(-b3 *BAL)))   ##basal area growth rate
  diametergrowthrate <- ((0.00545415*dbhinches^2 + bagr)/0.00545415)^.5 -dbhinches ##diameter growth rate inches
  ddgr <- diametergrowthrate * cminch    ##transform to cm
  return (ddgr)##returns diameter growth rate in cm
}








# From ARTEMIS
mortality <- function(g, dbh, ba){
  tmp1 <- g * ba
  cba <- rev(cumsum(rev(tmp1)))-tmp1  # cumulative basal area
  p <- 1-exp(-exp(-1.624 + (-2.6181) + (0.1229 * dbh) + (-0.8280 * log(dbh)) +
                    (0.0098 * 0 * log(10)) + (0.0208 * (cba)) + log(10)))
  pa <- p^1/10
  survan <- 1-pa
  return(survan)
}


##########################
## Postfire recruitment ##
##########################
SeedProd <- function(stand, baq){
  m2Ha <- 1e4
  Bd <- sum(baq[5:15] * stand[5:15]) / m2Ha  # pre-fire basal area/area
  Qd <- 163400 * Bd^0.95                     # germinable seeds/m2 in the aerial seed bank
  Pq <- 1                                    # proportion of seed abscised
  Sas <- 0.58                                # fraction of seed surviving pass fire
  m <- 0.0012                                # black spruce seed mass (g)
  w <- 0.14                                  # proportion of optimal seedbeds (Boiffin and Munson 2013)
  Sj <- 0.43 * (w * (1 - exp(-1.83 * m^0.43)) + (1 - w) * (1 - exp(-0.33 * m^0.76)))
  Fd <- Qd * Sj * Sas * Pq                   # the number of expected 3 year recuits/m2
  return(round(Fd * m2Ha))                   # per m2
}



# To estimate tree height, Peng, 1999 ##
Height <- function(dbh){
  1.3 + 1.065 * (dbh^0.8868)
}



# To estimate crown ratios, Holdaway ##
crownratio <- function(dbh, BA){
  b1 <- 5.54
  b2 <- 0.0072
  b3 <- 4.2
  b4 <- 0.053
  y <- b1 / (1 + (b2 * BA)) + (b3 * (1 - exp(-b4 * dbh)))
  round((y - 0.45) / 10, digits = 4)
}


# updates the crown ratios
updateCR <- function(pcr, ccr, dbh, TopHeight, di, BA, DBA, stand) {
  b1 <- 5.54
  b2 <- 0.0072
  b3 <- 4.2
  b4 <- 0.053
  if (sum(is.na(stand)) > 0)
    stop(message = "updateCR stand NA")
  if (sum(is.na(di)) > 0)
    stop(message = "updateCR di NA")
  if (sum(is.na(DBA)) > 0)
    stop(message = "updateCR DBA NA")
  d1 <- (b3 * b4 * exp(-b4 * dbh)) * di
  d2 <- (-b1 * b2 / ((1 + b2 * BA)^2)) * DBA
  xH <- TopHeight + di                       # Heightgrowth
  dH <- xH-TopHeight                         # deltaheight
  MaximumCR <- (ccr * TopHeight + dH) / xH   # maximumCR
  ccr <- ccr + d1 + ifelse(ccr < pcr, 0, d2)
  DBA <- rep(DBA, 15)
  newcr <- pmin(MaximumCR, ccr)
  ccr <- ifelse(DBA < 0, newcr, ccr)
  ccr <- ifelse(stand > 0, ccr, pcr)         # deals with all/ANY zeros
  ccr
}



#####################
## Carbon dynamics ##
#####################

Stemwood <- function(DBH) {
  stembiomass <- 0.0477 * DBH^2.5147
  return(stembiomass)
}


Bark <- function(DBH) {
  barkbiomass <- 0.0153 * DBH^2.2429
  return(barkbiomass)
}


Branches <- function(DBH) {
  branchbiomass <- 0.0278 * DBH^2.0839
  return(branchbiomass)
}


Needles <- function(DBH) {
  needlesbiomass  <- 0.1648 * DBH^1.4143
  return(needlesbiomass)
}

# Coarse root equation. Ouimet et al. 2008.
Coarse <- function(DBH){  # only roots >5 mm
  rootbiomass <- 0.0085 * 1.036 * (DBH^2.87)  # all classes
  fineroot <- 0.00153 * 2.40 * (DBH^1.123)     # 2mm-5mm
  coarsebiomass <- rootbiomass-fineroot
  return(coarsebiomass)
}


Fineroots <- function(dbh){  # <5mm. kg/tree Chen et al. 2004
  finerootbiomass <- 0.011 * (dbh^1.9748)
  return(finerootbiomass)
}


SnagsCarbon<-function(mortalitys,mortalityf,BioMassCarbon){
  kpob <- 0.25
  kpofr<- 0.045
  kpofo <- 1
  snags <- mortalitys
  snagsf <- mortalityf
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
  ###fire
  StemwoodbigF <- BioMassCarbon[1,5:15]%*%(snagsf[5:15])
  BarkmerchantableF <- BioMassCarbon[2,5:15]%*%(snagsf[5:15])
  StemwoodsmallF <- BioMassCarbon[1,1:4]%*%(snagsf[1:4])
  BarksmallF <- BioMassCarbon[2,1:4]%*%(snagsf[1:4])
  BranchesF <- BioMassCarbon[3,1:15]%*%(snagsf[1:15])
  FoliageF <- BioMassCarbon[4,1:15]%*%(snagsf[1:15])
  CRootF <- BioMassCarbon[5,1:15]%*%(snagsf[1:15])
  FRootF <- BioMassCarbon[6,1:15]%*%(snagsf[1:15])
  SnagCF <- StemwoodbigF+BarkmerchantableF
  SnagbranchCF <- BranchesF+StemwoodsmallF+BarksmallF
  SnagFoliageF <- FoliageF
  SnagCoarseF <- CRootF
  SnagFineF <- FRootF
  tmp1 <- FoliageF * kpofo
  tmp2 <-  FRootF * kpofr
  tmp3 <-   SnagbranchCF  *kpob
  SnagFoliageF <- SnagFoliageF -tmp1
  SnagFineF <- SnagFineF - tmp2
  SnagbranchCF <-  SnagbranchCF -  tmp3
  CE <- tmp1+tmp2+tmp3 
  
  return(list(SnagC=SnagC,SnagbranchC=SnagbranchC,SnagFoliage=SnagFoliage,SnagCoarse=SnagCoarse,SnagFine=SnagFine,
              SnagCF=SnagCF, SnagbranchCF=SnagbranchCF,SnagFoliageF=SnagFoliageF,
              SnagCoarseF=SnagCoarseF,SnagFineF=SnagFineF,CE=CE ))
}





#  Kurz et al 2013
Decayrates <- function(MAT) {
  REFT <- 10                                    # reference temperature
  Reduction <- MAT-REFT                         # reduction
  Q10 <- c(2, 2, 2, 2, 2.65, 2.65, 2, 2, 1)
  BDR <- c(0.0187, 0.072, 0.034, 0.1435, 0.355, 0.015, 0.5, 0.1435, 0.0033) #base decay rates at 10C
  TEMPMOD <- exp((Reduction) * log(Q10) * 0.1)
  STANDMOD <- 1
  ADR <- BDR * TEMPMOD * STANDMOD
  ADR
}



#################
## Fire module ##
#################
CanopyFuelStocksbs <- function(dbh) {
  w <- 0.6329 + 0.02 * dbh^2.2
  return(w)
}


# fuel per DBH class kg/ha  using Stocks allometry
FuelClass <- function(stand, dbh) {
  w <- CanopyFuelStocksbs(dbh) * stand
  return(w)
}


VerticalFuelProfile <- function(c, Top, Base) { # Base < Top must be guaranteed
  MaxTop <- max(Top)                         # maximum height of dbh class
  f <- numeric(MaxTop)                       # number of 1m height sections up to the maximum height
  for (i in seq(1:length(Top))) {             # for each  dbh class we do the following:
    v <- seq(Base[i] + 1, Top[i])          # 1 m sections per diameter class
    x <- 1 / length(v)                     # here we approximate the shape as cylinder.
    # each section gets the same proportion of biomass fuel
    f[v] <-  f[v] + (c[i] * x)             # the vector f accumulates  for each 1m section the fuel
    # of every dbh class per hectare
  }
  f / 10000                                  # convert to kg/m3
}


# critical crown base height for crowning based on the initial fire intensity (Van Wagner, 1973)
Zc <- function(I) {
  C <- 0.010
  m <- 100  #moisture
  h <- 460  +  26 * m
  criticalbase <- I^(2/3)/(C * h)
  return(criticalbase)
}

Crowning <- function(I, Top, Base, Bulk, b, DenCrit) {
  zc <- ceiling(Zc(I))                       #Critical base height at intensity (I)
  cb <- Base[Base <= zc]                     # which dbh classes are burning
  ct <- Top[Base <= zc]                      # flames extend to top of ct...
  CrownLayer <- 0                            # initialise
  k <- length(cb)                            # k =  number of dbh classes that affected by fire
  if (k > 0) {
    fc <- ct[k]                            # Height of the last dbh class affected by surface fire.
    # Returns the  flame length value that corresponds to
    # the highest strata with fire. Maxflamelength
    ct <- Top[Base <= fc]                  #  We evaluate again the CBH of dbh classes. Which ones
    # have a CBH lower than the flame height. return the dbh
    # classes with crown base height smaller than  max flame
    kk <- length(ct)
    if (kk > 0) {                          # We evaluate then if crowning can be sustained
      for (i in seq(kk, 1, by = -1)) {   # we are checking from top to down
        t <- max(ct[i], b)             # maximum height of the last dbh class with crowning
        den <- sum(Bulk[b:t])/(t-b + 1)   #running mean "Available CBD for combustion"
        if (den > DenCrit) {
          CrownLayer <- i
          break
        }
      }
    }
  }
  CrownLayer
}


UpdateCrownLayer <- function(cl, Top, Base, Bulk, b, DenCrit) {
  newcl <- 0
  if (cl > 0) {
    for (i in seq(cl, 1, by = -1)) {
      t <- max(Top[i], b)
      den <- sum(Bulk[b:t])/(t-b + 1)
      if (den>DenCrit) {
        newcl <- i
        break
      }
    }
  }
  newcl
}


FlameLength <- function(I) {
  0.0775 * ((I)^0.46)
}

UpdateIntensity <- function(I, Ht) {
  h <- 0
  h <- h + (Ht * 0.5)                        # add 1/2 mean canopy height for crown fires (Byram)
  259.83 * (h^2.174)                         # flame length intensity relationship
}

ScorchHeight <- function(I) {
  0.1483 * (I^0.667) # Van Wagner (1973)
}

CrownKill <- function(I, Top, CCR) {
  z <- ScorchHeight(I)
  ht <- Top
  cbh <- ht * (1 - CCR)
  cl <- ht-cbh
  tmp <- z - (ht - cl)
  tcls <- ifelse(z < cbh, 0, tmp)
  cls <- ifelse(z > ht, cl, tcls)
  cvscy <- 100 * (cls / cl)
  return(cvscy)
}

ScorchMortality <- function(Bark, CK) {
  bc <- 6.316 * (1 - exp(-Bark))
  cc <-  -0.000535 * (CK^2) #
  1 / (1 + exp(-1.941 + bc +  cc ))
}

# severity of fire
Basalost <- function (stand, baq, newstand) {
  a <- sum(stand * baq)    # Total basal area before
  b <- sum(newstand * baq) # Total basal area after fire
  Percentage_T_Basalost <- 100 * (1 - (b / a))
  return(Percentage_T_Basalost)
}
