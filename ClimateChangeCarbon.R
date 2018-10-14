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



exe<-function(stand,Y,Weather){
  
  #FRIs <- c(rep(916,30), rep(716,30), rep(458,30), rep(300,30)) ##Path1
  #FRIs <- c(rep(916,30), rep(304,30), rep(241,30), rep(170,30))  ##Path2
  #FRIs <- c(rep(783,30), rep(212,30), rep(250,30), rep(248,30))  ##Path3
  #FRIs <- c(rep(152,30), rep(182,30), rep(289,30), rep(145,30))  ##Path4
  FRIs <- c(rep(783,30), rep(1112,30), rep(735,30), rep(404,30))  ##Path5 
  
  MgKg <- 1000    
  #CPool <- c(4.9,1.5,24.9,9.0,9.1,35.1,1.6,3.4,95.0)* MgKg Path1
  #CPool <- c(5.0,1.6,27.0,9.7,10.0,36.8,1.7,3.7,97.8)* MgKg ##Path2
  #CPool <- c(5.0,1.5,24.9,8.9,8.8,35.3,1.5,3.4,95.5)* MgKg ##Path3
  #CPool <- c(5.0,1.8,24.8,9.2,8.7,35.7,1.5,3.2,94.9)* MgKg ##Path4
  CPool <- c(4.8,1.6,22.3,8.2,8.1,33.1,1.4,3.1,91.4)* MgKg ##Path5
  

  Sampled <- subpop (path=c("Pathway1"), d=PoolPlots2)
  Sampled <- Sampled[sample(1:dim(Sampled)[1], size=1, replace=T),]
  PlotID <- Sampled[1,2]
  WeatherPlot <- Sampled[1,1]
  Latitude <- Sampled[1,4]
  Longitude <- Sampled[1,5]
  ClimateData <- GetClimate (Weather,WeatherPlot)
  DroughtCodes <- GetDC(DroughtCode, WeatherPlot) 
  DroughtCode.List <- CorrectDupli(DroughtCodes)
  ClimateData.List <- ClimateData 
  ClimateData.List <- CorrectDupli(ClimateData) ##to fix duplicates
  MAT <- as.numeric(ClimateData.List[,3]) ## Mean annual temperatures 1981-2100
  DC <- as.numeric(DroughtCode.List[,3]) ## Projected Drought Codes 1981-2100
  GrowthIndex <- GI(ClimateData.List)  ##Growth Index
  GrowthIndex[1] <- 13.7
  Tree.List <- GetTrees (Tree,Plots=PlotID)
  as.character(Tree.List$DBH)
  as.factor(Tree.List$ESSENCE)
  Tree.List <- as.data.frame(lapply(Tree.List[,],function(x)rep(x,25)))
  range.DBH <- c(seq(1,30, by=2), 100)
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
  N1s <- rep(1, n)
  N0s <- rep(0, n)
  BurnMortP <- numeric(n)
  BurnMortPpar <- numeric(n)
  shannon <- 0
  ScH <- 0
  cl <- 0
  dbhl  <- c(1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29)  # dbh lower limit
  dbhu <- c(3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31)  # dbh upper limit
  dbhq  <- sqrt((dbhu^3-dbhl^3)/((dbhu-dbhl) * 3))  # assuming a uniform dbh distribution over interval
  baq  <- (dbhq/2)^2 * pi/1e4
  
  # Initialise variables
  Top <- ceiling(Height(dbhq))
  Top <- ifelse(Top<2, rep(2, length(Top)), Top)
  HeightUpdated <- Top
  PCR <- crownratio(dbhq, sum(baq * stand))
  CCR <- PCR
  Base <- floor(HeightUpdated * (1-CCR))
  Bark <- 0.032 * 2.54 * dbhq ##inches to cm!! Black spruce Equation taken from Behave plus fire modelling system  (thickness in cm)
  hatom2 <- 1e4
  b <- 3
  DenCrit <- 0.11
  iw <- dbhu-dbhl
  BarkThickness <- 0.032 * 2.54 * dbhq
  m2Ha  <- 1e4
  biocar_factor <- 0.5
  MgKg <- 1000
  SaplingSurvival <- 0.98  # sapling survival of different initial sizes (Matthias et al. 2003)
  NF <- 0
  fire.year <- NULL
  RegenLagTime <- 27   # assume it takes a 3-yr old natural seedling 27 yearss to grow to
  # 214 cm (DBH 1.0 cm), VanBoagaert et al. 2015)
  RegenRegular <- 60   # natural regeneration
  RegenIrregular <- 55 # natural regeneration
  shannon <- diversity(stand[5:15], index = "shannon", MARGIN = 1, base = exp(1))  # shannon index
  RegenCohorts <- rpois(RegenLagTime, ifelse(shannon<1.7, RegenRegular, RegenIrregular))
  Snags <- 0
  Snagbranches <- 0
  SnagFoliage <- 0
  SnagCoarse <- 0
  SnagFine <- 0
  CC <- 0
  
  ###Initialise object variables to save simulation results####
  
  Size <- matrix(0, Y, n, byrow = TRUE)
  PreFireStand <- matrix(0, Y, n, byrow = TRUE)
  Recruits <- numeric(Y)
  BA <- numeric(Y)
  FireDeaths <- matrix(0, Y, n, byrow = TRUE)
  Senescence <- matrix(0, Y, n, byrow = TRUE)
  Transition <- matrix(0, Y, n, byrow = TRUE)  # transition probability
  Crecimiento <- matrix(0, Y, n, byrow = TRUE)
  Parcela <- matrix(0, Y, n, byrow = TRUE)
  Muertos <- matrix(0, Y, n, byrow = TRUE)
  Mortality <- matrix(0, Y, n, byrow = TRUE)  # probability of mortality
  Structure <- numeric(Y)
  InitialIntensity <- numeric(Y)
  Latitudes <- numeric(Y)
  Longitudes <- numeric(Y)
  DroughtCode <- numeric(Y)
  EmpiricalTemperature <- numeric(Y)
  CCGrowthIndex <- numeric(Y)
  FireSeason <- numeric(Y)
  BALost <- numeric(Y)
  Delta_BA <- numeric(Y)
  DeltaN <- matrix(0, Y, n, byrow = TRUE)
  CR <- matrix(0, Y, n, byrow = TRUE)
  Heights <- matrix(0, Y, n, byrow = TRUE)
  ShiftCrownratio <- matrix(0, Y, n, byrow = TRUE)
  ShiftHeights <- matrix(0, Y, n, byrow = TRUE)
  DiameterGrowth <- matrix(0, Y, n, byrow = TRUE)
  SnagCProduction <- numeric(Y)
  SnagCProductionS <- numeric(Y)
  SnagCProductionF <- numeric(Y)
  Turnover <- numeric(Y)
  DOMC_Pool <- matrix(0, Y, length(CPool), byrow = TRUE)
  DOM_Flux <- numeric(Y)
  DOM_Inputs <- matrix(0, Y, 9, byrow = TRUE)
  
  
  BioMass <- matrix(c(Stemwood(dbhq), Bark(dbhq), Branches(dbhq), Needles(dbhq), Coarse(dbhq), Fineroots(dbhq)), nrow = 6,
                    ncol = length(dbhq), byrow = TRUE)
  BioMassCarbon <- BioMass * biocar_factor  # Biomass C per diameter class
  
  
  Merchantable <- c(0,0,0,0,Stemwood(dbhq[5:15])+Bark(dbhq[5:15]))
  OtherWood <- c(Stemwood(dbhq[1:4])+Bark(dbhq[1:4]),rep(0,11))+ Branches(dbhq)
  
  InitialCBiomass <- sum(BioMassCarbon%*%as.matrix(stand)) ##bio1 for the site
  ICB <- InitialCBiomass
  NetPrimaryProductivity <- numeric(Y)
  TotalLiveBiomass <- numeric(Y)
  AppliedDecayRates <- matrix(0, Y, 9, byrow = TRUE)
  NetEcosystemProduction <- numeric(Y)
  Rh <- numeric(Y)
  MFRI <- numeric(Y)
  CarbonCombusted <- numeric(Y)
  AnnualBiomassRecruits <- numeric(Y)
  CarbonEmissions1 <- numeric(Y) 
  CarbonEmissions2 <- numeric(Y) 
  CarbonEmissions3 <- numeric(Y) 
  NetBiomeProduction <- numeric(Y) 
  FuelConsumed <- numeric(Y) 
  PlotIDs <- numeric(Y)
  LPeriod <- as.numeric(Y)
  Periods <- matrix(0,120,Y, byrow=TRUE)  
  
  
  for (y in 1:Y){ #Things that have to be reinitialized and updated
    
    HeadIntensity <- sample(IntensitiesWeighted, size=1, replace=F)
    I <- HeadIntensity
    Istart <- I
    shannon <- diversity(stand[5:15], index="shannon", MARGIN=1, base=exp(1)) ##updated shannon
    Structure[y] <- shannon
    Firedeaths <- rep(0,15)
    BALost[y] <- 0
    bay <- sum(stand*baq)
    SnagsS <- 0
    SnagbranchesS <- 0
    SnagFoliageS <- 0
    SnagCoarseS <- 0
    SnagFineS <- 0
    SnagsFire <- 0
    SnagbranchesFire <- 0
    SnagFoliageFire <- 0
    SnagCoarseFire <- 0
    SnagFineFire <- 0
    carbon.consumed.medium <- 0
    carbon.consumed.agfast <- 0
    carbon.consumed.agvfast <- 0
    carbon.consumed.agslow <- 0
    carbon.emissions <- 0
    CE <- 0
    C <- 0
    CC2<- 0
    BALost[y] <- 0
    bay <- sum(stand * baq)
    GI <- GrowthIndex[y]
    AppDecayRates <- Decayrates(MAT[y])
    EmpiricalTemperature[y] <- MAT[y]
    FRI <- FRIs[y]
    MFRI[y] <- FRI
    pBurn <- 1/FRI
    Latitudes [y] <- Latitude
    Longitudes [y] <- Longitude
    Decayrate <- rep(0, 9)
    Decayrate[1] <- AppDecayRates[1]
    Decayrate[2] <- AppDecayRates[2]
    Decayrate[3] <- AppDecayRates[3]
    Decayrate[4] <- AppDecayRates[4]
    Decayrate[5] <- AppDecayRates[5]
    Decayrate[6] <- AppDecayRates[6]
    Decayrate[7] <- AppDecayRates[7]
    Decayrate[8] <- AppDecayRates[8]
    Decayrate[9] <- AppDecayRates[9]
    
    # Apply decay rates to CPools
    CarbonPoolTransferMatrix <- matrix(c(
      Decayrate[1] * 0.83, (1-Decayrate[1]-0.032), 0, 0.032, 0, 0, Decayrate[1] * (1-0.83), 0, 0, 0,
      Decayrate[2] * 0.83, 0, (1-Decayrate[2]-0.10), 0, 0.10, 0, Decayrate[2] * (1-0.83), 0, 0, 0,
      Decayrate[3] * 0.83, 0, 0, 1-Decayrate[3], 0, 0, Decayrate[3] * (1-0.83), 0, 0, 0,
      Decayrate[4] * 0.83, 0, 0, 0, (1-Decayrate[4]), 0, Decayrate[4] * (1-0.83), 0, 0, 0,
      Decayrate[5] * 0.815, 0, 0, 0, 0, (1-Decayrate[5]), Decayrate[5] * (1-0.815), 0, 0, 0,
      Decayrate[6] * 1, 0, 0, 0, 0, 0, (1-Decayrate[6]-0.006), 0, 0, 0.006,
      Decayrate[7] * 0.83, 0, 0, 0, 0, 0, 0, (1-Decayrate[7]), 0, Decayrate[7] * (1-0.83),
      Decayrate[8] * 0.83, 0, 0, 0, 0, 0, 0, 0, (1-Decayrate[8]), Decayrate[8] * (1-0.83),
      Decayrate[9] * 1, 0, 0, 0, 0, 0, 0, 0, 0, 1-Decayrate[9]
    ), nrow = 9, ncol = 10, byrow = TRUE)
    colnames(CarbonPoolTransferMatrix) <- c("Atm", "Snags", "Snagbranch", "Medium",
                                            "AGfast", "AGveryfast", "AGslow", "BGveryfast", "BGfast", "BGslow")
    rownames(CarbonPoolTransferMatrix) <- c("Snags", "Snagbranch", "Medium", "AGfast", "AGveryfast", "AGslow",
                                            "BGveryfast", "BGfast", "BGslow")
    
    tmp <- as.vector(t(CPool)%*%CarbonPoolTransferMatrix)
    SoilCAtmFlux <- tmp[1]
    
    
    
    # Stand dynamics
    CCGrowthIndex[y] <- GI
    annual_diam_incre <- DiameterGrowthRate(stand,dbhq,baq,GrowthIndex[y])
    adi<-annual_diam_incre
    graduating <- 1 / (iw / adi) #Transition probabilities
    growth <- rbinom(N1s, stand, graduating) #stochastically
    stand <- stand - growth
    stand <- stand + c(0, growth[1:n - 1]) #after growth
    LastDCT <- growth[15]  # prevents loosing trees in last diameter class
    stand [15] <- stand[15] + LastDCT
    BAIncrement <- sum(growth* baq)    # patch BA increment due to growth
    CCRgrowth <- growth * CCR           # the trees that grow bring their crown ratio
    Heightgrowth <- growth * HeightUpdated
    GrowthCBiomass<- sum(BioMassCarbon%*%as.matrix(stand)) ##calculate biomass due to growth..Biomass that has not been lost to turnover or mortality
    
    # Calculate Biomass due to growth
    GrowthCBiomass <- sum(BioMassCarbon %*% as.matrix(stand))  # calculate biomass due to growth
    delta <- (GrowthCBiomass-ICB) # Biomass that has not been lost to turnover or mortality
    
    # Apply turnover
    # match IPCC Good Practice Guidance
    Stemwoodsmall <- BioMassCarbon[1, 1:4] %*% (stand[1:4])
    Barkmerchantable <- BioMassCarbon[2, 5:15] %*% (stand[5:15])
    LiveBiomassCPools <- BioMassCarbon %*% (stand) # biomass C kg/ha
    LiveBiomassCPoolsCorrected  <- matrix(0, nrow = 5, ncol = 1)
    LiveBiomassCPoolsCorrected[1, 1] <- LiveBiomassCPools[1, 1] - Stemwoodsmall + Barkmerchantable
    LiveBiomassCPoolsCorrected[2, 1] <- LiveBiomassCPools[2, 1] - Barkmerchantable + LiveBiomassCPools[3, 1] +
      Stemwoodsmall   #Otherwood + Bark
    LiveBiomassCPoolsCorrected[3, 1] <- LiveBiomassCPools[4, 1]  # Foliage
    LiveBiomassCPoolsCorrected[4, 1] <- LiveBiomassCPools[5, 1]  # Coarse
    LiveBiomassCPoolsCorrected[5, 1] <- LiveBiomassCPools[6, 1]  # Fine
    
    
    # Apply transfer rates
    Input_Matrix2 <- matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                              0, (1 - 0.04), 0, 0, 0, 0, (0.04 * 0.25), 0, (0.04 * 0.75), 0, 0, 0, 0, 0,
                              0, 0, (1 - 0.16), 0, 0, 0, 0, 0, 0, (0.16 * 1), 0, 0, 0, 0,
                              0, 0, 0, (1 - 0.02), 0, 0, 0, 0, (0.02 * 0.5), 0, 0, 0, (0.02 * 0.5), 0,
                              0, 0, 0, 0, (1 - 0.64), 0, 0, 0, 0, (0.64 * 0.5), 0, (0.64 * 0.5), 0, 0), nrow = 5, ncol = 14, byrow = TRUE)
    
    colnames(Input_Matrix2) <- c("Merchantable", "Otherwood", "Needles",
                                 "Coarse", "Fine", "Snags", "Snagbranch",
                                 "Medium", "AGfast", "AGveryfast", "AGslow",
                                 "BGveryfast", "BGfast", "BGslow")
    rownames(Input_Matrix2) <- c("Merchantable", "Otherwood",
                                 "Needles", "Coarse", "Fine")
    
    
    ctmp <- as.vector(t(LiveBiomassCPoolsCorrected)%*%Input_Matrix2)  # turnover
    BiomassLostTurnover <- sum(ctmp[6:14])
    
    
    # Carbon fluxes
    NPP <- delta  +  BiomassLostTurnover
    NEP <- (NPP-SoilCAtmFlux)
    
    
    # Mortality (senescence)
    surviving <- mortality(stand, dbhq, baq)
    surviving[1:4] <- SaplingSurvival
    Senescencedeaths <- rbinom(N1s, stand, 1 - surviving)
    Senescence[y, ] <- Senescencedeaths
    stand <- stand-Senescencedeaths # fire deaths are taken care in the fire module
    
    
    
    ##FIRE MODULE
    # Fire module
    # Evaluates if a fire arrives, its intensity, crowning, and post-fire tree mortality and regeneration.
    if(runif(1) < pBurn) {  # Determine if a fire happens
      FuelLoad <- sum(CPool[5],CPool[6])  # only AGslow and AGvf
      fire.year <- c(fire.year, y)
      NF <- NF + 1  # Update number of fires that occurred during the simulation
      RegenCohorts <- rep(0, RegenLagTime)  # KILL ALL regenerating trees
      prefirestand <- stand
      PreFireStand[y, ] <- prefirestand
      NewRegen <- SeedProd(prefirestand, baq)
      Fuel <- as.numeric(FuelClass(stand, dbhq))  # Kg/ha per dbh class
      VF <- VerticalFuelProfile(Fuel, HeightUpdated, Base)
      cl <- Crowning(I, HeightUpdated, Base, VF, b, DenCrit) #if there is crowning, update the crown layer affected by fire
      cl <- ifelse (cl>0, UpdateCrownLayer(cl, HeightUpdated, Base, VF, b, DenCrit),
                    Crowning(I, HeightUpdated, Base, VF, b, DenCrit))
      u <- UpdateIntensity(I, HeightUpdated[cl])
      I <- ifelse(cl > 0, max(UpdateIntensity(I, HeightUpdated[cl])), Istart)
      ScH <- ScorchHeight(I)
      CK <- CrownKill(I, HeightUpdated, CCR)
      BurnMortP <- ScorchMortality(BarkThickness, CrownKill(I, HeightUpdated, CCR))
      Firedeaths <- rbinom(N1s, stand, BurnMortP)
      Firedeaths <- ifelse(is.na(Firedeaths), N0s, Firedeaths)
      FireDeaths[y, ] <- Firedeaths
      newstand <- stand-Firedeaths
      severity <- Basalost(stand, baq, newstand)
      stand <- newstand # update stand after a fire
      InitialIntensity[y] <- I # adjusted intensity using Catchpole et.al 1992
      DroughtCode[y] <- DC[y]
      BALost[y] <- severity
      C <- 1.185*exp(-4.252)*exp(0.671*log(FuelLoad))*exp(0.71*log(DC[y]))
      pFF <- min(C/(CPool[5] + CPool[6]),1)
      carbon.consumed.medium <- CPool[3]*0.392140  #carbon consumed in the medium carbon pool
      carbon.consumed.agfast <- CPool[4]*0.6415    #carbon consumed in the Ag fast carbon pool
      carbon.consumed.agvfast <- pFF * CPool[5]    #carbon consumed in the Ag very fast carbon pool
      carbon.consumed.agslow <- pFF * CPool[6]     #carbon consumed in the Ag slow carbon pool
      CPool[3] <- CPool[3]- carbon.consumed.medium
      CPool[4] <- CPool[4]- carbon.consumed.agfast
      CPool[5] <- CPool[5]- carbon.consumed.agvfast #carbon left in the Ag very fast pool
      CPool[6] <- CPool[6]- carbon.consumed.agslow  #carbon left in the Ag slowpool
      FuelConsumed [y] <- C
    } else {
      NewRegen <- rpois(1, ifelse(shannon < 1.7, RegenRegular, RegenIrregular))
    }
    
    # Calculate carbon from fire-derived and senescencent trees
    
    deaths <- Senescencedeaths + Firedeaths
    Muertos[y,] <- deaths
    DeltaBA <- BAIncrement - sum(deaths*baq) ##net basal area increment after mortality (fire and senescence) and growth
    DeltaStand <- sum(stand)  ##net change in density after mortality (fire and senescence) and growth
    if(sum(stand<0)>0)
      stop(message="neg count 3")
    Delta_BA[y] <- DeltaBA
    DeltaN[y,] <- deaths+growth
    SnagCpools <- SnagsCarbon(Senescencedeaths,Firedeaths,BioMassCarbon)
    SnagsS <- SnagCpools$SnagC
    SnagbranchesS <- SnagCpools$SnagbranchC
    SnagFoliageS <- SnagCpools$SnagFoliage
    SnagCoarseS <- SnagCpools$SnagCoarse
    SnagFineS <- SnagCpools$SnagFine
    CMortalityS <- SnagsS+SnagbranchesS+SnagFoliageS+SnagCoarseS+ SnagFineS 
    SnagsFire <- SnagCpools$SnagCF
    SnagbranchesFire <- SnagCpools$SnagbranchCF
    SnagFoliageFire <- SnagCpools$SnagFoliageF
    SnagCoarseFire <- SnagCpools$SnagCoarseF
    SnagFineFire <- SnagCpools$SnagFineF
    CMortalityF <-  SnagsFire+SnagbranchesFire+SnagFoliageFire+SnagCoarseFire+SnagFineFire
    Snags <- SnagsS+SnagsFire
    Snagbranches <- SnagbranchesS+SnagbranchesFire
    SnagFoliage <- SnagFoliageS+SnagFoliageFire
    SnagCoarse <- SnagCoarseS+SnagCoarseFire
    SnagFine <-  SnagFineS+SnagFineFire
    CMortality <- Snags + Snagbranches + SnagFoliage + SnagCoarse + SnagFine
    CE <-  SnagCpools$CE
    carbon.emissions <- carbon.consumed.medium+carbon.consumed.agfast+ carbon.consumed.agvfast+carbon.consumed.agslow+carbon.emissions ##fire carbon emissions total
    CC2 <- carbon.consumed.medium+carbon.consumed.agfast+ carbon.consumed.agvfast+carbon.consumed.agslow 
    NBP <- NEP - carbon.emissions
    
    # Add recruitment
    
    Recruits[y] <- RegenCohorts[RegenLagTime]
    RegenCohorts <- c(NewRegen, RegenCohorts[1:RegenLagTime - 1])
    stand[1] <- stand[1] + Recruits[y]
    CCRRecruits  <- Recruits[y] * PCR[1]
    HeightRecruits  <- Recruits[y] * Top[1]
    Parcela[y, ] <- stand  # stand after regeneration, captures regeneration pulses
    
    
    
    # Update biomass
    Biomass <- sum(BioMassCarbon %*% as.matrix(stand))
    ICB <- Biomass
    
    # Dynamically updating crown ratios and heights of recruits (natural regenerated and fire derived)
    xH <- Top + adi  # Heightgrowth
    dH <- xH - Top     # deltaheight
    MaximumCR <- (CCR * Top + dH) / xH
    TotalN <- stand + c(0, growth[1:n - 1])
    TotalN[1] <- TotalN[1] + Recruits[y]
    CCR <- updateCR(PCR, CCR, dbhq, Top, adi, bay, DeltaBA, stand)  # CCR after recruitment, growth, mortality
    CCRnow <- stand * CCR
    ShiftCR <- CCRnow + c(0, CCRgrowth[1:n - 1])
    ShiftCR[1] <- ShiftCR[1] + CCRRecruits
    CCR <- ifelse(stand > 0, pmin(MaximumCR, ShiftCR / TotalN), PCR)
    
    # Updating heights
    Top <- Height(dbhq)
    HeightUpdated <- Top + adi  # Heightgrowth
    Heightnow <- stand * HeightUpdated
    ShiftHeight <- Heightnow + c(0, Heightgrowth[1:n-1])
    ShiftHeight[1] <- ShiftHeight[1] + HeightRecruits
    HeightUpdated <- ifelse(stand>0, ShiftHeight/TotalN, Top)
    Base <- HeightUpdated * (1-CCR)
    
    
    # Distribute turnover to carbon pools
    ctmp[6] <- ctmp[6] + Snags             # adding C from snags
    ctmp[7] <- ctmp[7] + Snagbranches      # adding C from small trees
    ctmp[9] <- ctmp[9] + SnagCoarse * (0.5) # adding C from coarse roots to AG fast
    ctmp[10] <- ctmp[10] + SnagFine * (0.5) + SnagFoliage 
    ctmp[12] <- ctmp[12] + SnagFine * (0.5)
    ctmp[13] <- ctmp[13] + SnagCoarse * (0.5)
    Inputs <- ctmp[6:14]  # how much C is incorporated into the DOMCpools including snags (mortality)
    CPool <- tmp[2:10] + Inputs # update the pools
    
    
    ## Save results into objects
    Delta_BA[y] <- DeltaBA
    DeltaN[y, ] <- deaths + growth
    Muertos[y, ] <- deaths
    SnagCProduction [y] <- CMortality
    AppliedDecayRates[y, ] <- AppDecayRates
    SnagCProductionS [y] <- CMortalityS
    SnagCProductionF [y] <- CMortalityF
    CarbonEmissions3[y] <- carbon.emissions
    TotalLiveBiomass[y] <-   Biomass
    DOMC_Pool[y, ] <- CPool # Nine DOM carbon pools
    DOM_Flux[y] <- SoilCAtmFlux
    DOM_Inputs[y, ] <- Inputs
    Turnover[y] <- BiomassLostTurnover
    NetPrimaryProductivity[y] <- NPP
    NetEcosystemProduction[y] <- NEP
    NetBiomeProduction[y] <- NBP
    Rh[y] <- SoilCAtmFlux
    BA[y] <- bay
    CR[y, ] <- CCR
    Heights[y, ] <- HeightUpdated
    ShiftCrownratio[y, ] <- ShiftCR
    ShiftHeights[y, ] <- ShiftHeight
    DiameterGrowth[y, ] <- adi
    Transition[y, ] <- graduating
    Mortality[y, ] <- 1-surviving
    Crecimiento[y, ] <- growth
    Muertos[y, ] <- deaths
    Size[y, ] <- stand # final size after all processes
    
  }
  
  res<-list(MFRI=MFRI, Latitudes= Latitudes, Longitudes=Longitudes, Parcela=Parcela,Size=Size,
            DeltaN=DeltaN,BA=BA,PreFireStand=PreFireStand,Senescence=Senescence,
            AppliedDecayRates=AppliedDecayRates,EmpiricalTemperature=EmpiricalTemperature, 
            CCGrowthIndex =CCGrowthIndex,FireDeaths=FireDeaths,Muertos=Muertos,Structure=Structure,
            Crecimiento=Crecimiento,Mortality=Mortality,Transition=Transition,DiameterGrowth=DiameterGrowth,
            CR=CR,ShiftCrownratio=ShiftCrownratio,ShiftHeights=ShiftHeights,Heights=Heights,
            BALost=BALost, NF=NF,InitialIntensity=InitialIntensity,fire.year=fire.year,
            Recruits=Recruits,TotalLiveBiomass=TotalLiveBiomass,Turnover=Turnover,
            SnagCProduction=SnagCProduction,DOMC_Pool=DOMC_Pool,DOM_Flux=DOM_Flux,
            DOM_Inputs=DOM_Inputs,NetPrimaryProductivity=NetPrimaryProductivity, 
            Rh=Rh,NetEcosystemProduction=NetEcosystemProduction, CarbonEmissions1=CarbonEmissions1,
            CarbonEmissions2= CarbonEmissions2,CarbonEmissions3=CarbonEmissions3,NetBiomeProduction=NetBiomeProduction,
            FuelConsumed=FuelConsumed,SnagCProductionS=SnagCProductionS,SnagCProductionF=SnagCProductionF,
            SnagCProduction=SnagCProduction, DroughtCode=DroughtCode)
  
  return(res)
}









####rate of change
#
#
#ratechange<-function (NetEcosystemProduction_s){
#    rc<-matrix(0,1000,29,byrow=TRUE)
#    for (i in  1:nrow(NetEcosystemProduction_s)){
#    change<-diff(NetEcosystemProduction_s[i,])
#    rc[i,]<- change
#    }
#    rc
#    }
#
#RatechangeHis<-ratechange(NetEcosystemProduction_s)
#save("RatechangeHis", file="F:\\Universite_Laval\\YOSDATA\\LAVAL\\ModelOutput3\\Simulations\\RatechangeHis.RData")
#RatechangeSecond<-ratechange(NetEcosystemProduction_s)
#save("RatechangeSecond", file="F:\\Universite_Laval\\YOSDATA\\LAVAL\\ModelOutput3\\Simulations\\Group2\\Second\\RatechangeSecond.RData")
#RatechangeThird<-ratechange(NetEcosystemProduction_s)
#save("RatechangeThird", file="F:\\Universite_Laval\\YOSDATA\\LAVAL\\ModelOutput3\\Simulations\\Group1\\Third\\RatechangeThird.RData")
#RatechangeFourth<-ratechange(NetEcosystemProduction_s)
#save("RatechangeThird", file="F:\\Universite_Laval\\YOSDATA\\LAVAL\\ModelOutput3\\Simulations\\Group1\\Fourth\\RatechangeFourth.RData")
#
#
#
#cc<-colMeans(RatechangeThird)
#mean(cc)
#
#yos2<-exe(stand,Y=10000,FRI=120000)
#KgMg<-(0.001)
#NPP<-(yos2$NetPrimaryProductivity)*KgMg # after mortality and turnover
#HeterotrophicRespiration<-(yos2$Rh)*KgMg
#NEP<-(yos2$NetEcosystemProduction)*KgMg
#
#plot(NPP, ylim=c(-5,8),xlab="Simulation period",ylab="Fluxes(MgC/ha*yr)", type="l", col="black",lwd=2)
#lines(HeterotrophicRespiration,,lty=1,col="red",lwd=2)
#lines(NEP,,lty=1,col="blue",lwd=2)
#abline(h = 0,  col = "gray60")
#legend("topleft",legend=c("NPP","Rh","NEP")
#,col=c("black","red","blue"),lty=1,lwd=2)
#
#turn<-yos2$Turnover
#icb<-yos2$CarbonBiomass2
#GrowthC<-yos2$CarbonBiomass1
#netga<-icb-turn
#netga+turn=icb
#
#

#Temperature_8.5_2011_2040<- read.csv("C:\\Users\\yomiq\\Documents\\YOSDATA\\Laval\\ModelOutput2\\TemperatureFiles\\2011-2040_HadGEM2-ES-RCP-8.5-2011-2040_LONLAT_1_mod.csv",header=TRUE,colClasses = "character",sep = ",")
#MeanAnnualTemperaturePrec8.5_2011_2040<-subset(Temperature_8.5_2011_2040, select=c("longitude", "latitude","temp?rature.moyenne.annuelle","Pr?cipitations.annuelles"))
##MeanAnnualPrecipitation8.5_2011_2040<-subset(Temperature_8.5_2011_2040, select=c("longitude", "latitude","Pr?cipitations.annuelles"))
#MeanAnnualTemperaturePrec8.5_2011_2040$TScenario<-"RCP8.5"
#MeanAnnualTemperaturePrec8.5_2011_2040$temp?rature.moyenne.annuelle<-as.numeric (MeanAnnualTemperature8.5_2011_2040$temp?rature.moyenne.annuelle)
#range(MeanAnnualTemperaturePrec8.5_2011_2040$temp?rature.moyenne.annuelle)
#MeanAnnualTemperaturePrec8.5_2011_2040$Pr?cipitations.annuelles<-as.numeric(MeanAnnualPrecipitation8.5_2011_2040$Pr?cipitations.annuelles)
#range(MeanAnnualTemperaturePrec8.5_2011_2040$Pr?cipitations.annuelles)
#
#Temperature_4.5_2011_2040<-read.csv("C:\\Users\\yomiq\\Documents\\YOSDATA\\Laval\\ModelOutput2\\TemperatureFiles\\2011-2040_HadGEM2-ES-RCP-4.5-2011-2040_LONLAT_1_mod.csv",header=TRUE,colClasses = "character",sep = ",")
#MeanAnnualTemperaturePrec4.5_2011_2040<-subset(Temperature_4.5_2011_2040, select=c("longitude", "latitude","temp?rature.moyenne.annuelle","Pr?cipitations.annuelles"))
#MeanAnnualTemperaturePrec4.5_2011_2040$TScenario<-"RCP4.5"
#MeanAnnualTemperaturePrec4.5_2011_2040$temp?rature.moyenne.annuelle<-as.numeric (MeanAnnualTemperature4.5_2011_2040$temp?rature.moyenne.annuelle)
#range(MeanAnnualTemperature4.5_2011_2040$temp?rature.moyenne.annuelle)
##MeanAnnualPrecipitation4.5_2011_2040<-subset(Temperature_4.5_2011_2040, select=c("longitude", "latitude","Pr?cipitations.annuelles"))
#MeanAnnualTemperaturePrec4.5_2011_2040$Pr?cipitations.annuelles<-as.numeric(MeanAnnualPrecipitation4.5_2011_2040$Pr?cipitations.annuelles)
#range(MeanAnnualTemperaturePrec4.5_2011_2040$Pr?cipitations.annuelles)
#
#
##################################
#
#Temperature_2.6_2011_2040<-read.csv("C:\\Users\\yomiq\\Documents\\YOSDATA\\Laval\\ModelOutput2\\TemperatureFiles\\2011-2040_HadGEM2-ES-RCP-2.6-2011-2040_LONLAT_1_mod.csv",header=TRUE,colClasses = "character",sep = ",")
#MeanAnnualTemperaturePrec2.6_2011_2040<-subset(Temperature_2.6_2011_2040, select=c("longitude", "latitude","temp?rature.moyenne.annuelle","Pr?cipitations.annuelles"))
#MeanAnnualTemperaturePrec2.6_2011_2040$TScenario<-"RCP2.6"
#MeanAnnualTemperaturePrec2.6_2011_2040$temp?rature.moyenne.annuelle<-as.numeric (MeanAnnualTemperature2.6_2011_2040$temp?rature.moyenne.annuelle)
#range(MeanAnnualTemperaturePrec2.6_2011_2040$temp?rature.moyenne.annuelle)
##MeanAnnualPrecipitation2.6_2011_2040<-subset(Temperature_2.6_2011_2040, select=c("longitude", "latitude","Pr?cipitations.annuelles"))
#MeanAnnualTemperaturePrec2.6_2011_2040$Pr?cipitations.annuelles<-as.numeric(MeanAnnualPrecipitation2.6_2011_2040$Pr?cipitations.annuelles)
#range(MeanAnnualTemperaturePrec2.6_2011_2040$Pr?cipitations.annuelles)
#
#
#
###################################
#Temperature_Historical<-read.csv("C:\\Users\\yomiq\\Documents\\YOSDATA\\Laval\\ModelOutput2\\TemperatureFiles\\1971-2000_mean_LONLAT.csv",header=TRUE,colClasses = "character",sep = ",")
#MeanAnnualTemperatureHistorical<-subset(Temperature_Historical, select=c("longitude", "latitude","temp?rature.moyenne.annuelle","Pr?cipitations.annuelles"))
#MeanAnnualTemperatureHistorical$TScenario<-"Historical"
#Temperature_1981_2010<-read.csv("C:\\Users\\yomiq\\Documents\\YOSDATA\\Laval\\ModelOutput2\\TemperatureFiles\\1981-2010_mean_LONLAT_1_mod.csv",header=TRUE,colClasses = "character",sep = ",")
#MeanAnnualTemperature2<-subset(Temperature_1981_2010, select=c("longitude", "latitude","temp?rature.moyenne.annuelle","Pr?cipitations.annuelles"))
#Histor<-merge(MeanAnnualTemperatureHistorical,MeanAnnualTemperature2)
#Histor$temp?rature.moyenne.annuelle<-as.numeric(Histor$temp?rature.moyenne.annuelle)
#mean(Histor$temp?rature.moyenne.annuelle )
#Histor$Pr?cipitations.annuelles<-as.numeric(Histor$Pr?cipitations.annuelles)
#
#
##############################################################
#Tem<-rbind(MeanAnnualTemperatureHistorical,MeanAnnualTemperaturePrec2.6_2011_2040,MeanAnnualTemperaturePrec4.5_2011_2040,MeanAnnualTemperaturePrec8.5_2011_2040)
#length(1:dim(Tem)[1])
#Tem$TScenario<-as.factor(Tem$TScenario)
#Tem$temp?rature.moyenne.annuelle<-as.numeric(Tem$temp?rature.moyenne.annuelle)
#Tem$Pr?cipitations.annuelles<-as.numeric(Tem$Pr?cipitations.annuelles)
#
#fit1<-lm(temp?rature.moyenne.annuelle~TScenario, data=Tem)
#anova(fit1)
#summary(fit1)
#aov1<-aov(temp?rature.moyenne.annuelle~TScenario,data=Tem)
#summary(aov1)
#plot(aov1)
#posthoc <- TukeyHSD(x=aov1, 'TScenario', conf.level=0.95)
#plot(posthoc)
#
#
#fit2<-lm(Pr?cipitations.annuelles~TScenario, data=Tem)
#anova(fit2)
#summary(fit2)
#
#par(mfrow = c(1, 2))
#boxplot(temp?rature.moyenne.annuelle~TScenario, data=Tem, xlab="Climatic scenario", ylab="Mean annual temperature (?C)")
#boxplot(Pr?cipitations.annuelles~TScenario, data=Tem, xlab="Climatic scenario", ylab="Mean annual precipitation (mm)")
#
#m<-tapply(Tem$temp?rature.moyenne.annuelle,Tem$TScenario,mean)
#p<-tapply(Tem$Pr?cipitations.annuelles,Tem$TScenario,mean)
#
#result <- cbind(TemperatureMaxFebruary2011_2040, TemperatureMinAugust2011_2040)
#TemperatureMinJuneAugust2011_2040,SumPrecipSeptember2011_2040,SumPrecipMayAugust2011_2040,
#head(result)
#length(1:dim(Climatic2041_2070)[1])
#FebTMax<-as.vector(tapply(Climatic2041_2070$TMaxFeb...C.,Climatic2041_2070$Year,mean))
#AugusTMin<-as.vector(tapply(Climatic2041_2070$TAugustMin...C.,Climatic2041_2070$Year,mean))
#JuneAugusTMin<-as.vector(tapply(Climatic2041_2070$TMinJuneAugust...C.,Climatic2041_2070$Year,mean))
#PrecipSept<-as.vector(tapply(Climatic2041_2070$TotalPrecipitationSeptember..mm.,Climatic2041_2070$Year,mean))
#PreciMayAugus<-as.vector(tapply(Climatic2041_2070$TotalPrecipitationMayAugust..mm.,Climatic2041_2070$Year,mean))
#


#Climatic2041_2070<-read.csv("C:\\Users\\yomiq\\Documents\\YOSDATA\\Laval\\ModelOutput3\\ClimateData2041_2070RCP8.5.csv",header=TRUE,colClasses = "character",sep = ",")
#Climatic2041_2070[,3:7] <- sapply(Climatic2041_2070[,3:7],as.numeric)
#str(Climatic2041_2070)
#length(1:dim(Climatic2041_2070)[1])
#Years<-rep(c(2041:2070),30)
#Climatic2041_2070$Year<-Years
#names(Climatic2041_2070)
#FebTMax<-as.vector(tapply(Climatic2041_2070$TMaxFeb...C.,Climatic2041_2070$Year,mean))
#AugusTMin<-as.vector(tapply(Climatic2041_2070$TAugustMin...C.,Climatic2041_2070$Year,mean))
#JuneAugusTMin<-as.vector(tapply(Climatic2041_2070$TMinJuneAugust...C.,Climatic2041_2070$Year,mean))
#PrecipSept<-as.vector(tapply(Climatic2041_2070$TotalPrecipitationSeptember..mm.,Climatic2041_2070$Year,mean))
#PreciMayAugus<-as.vector(tapply(Climatic2041_2070$TotalPrecipitationMayAugust..mm.,Climatic2041_2070$Year,mean))
#
#
###########################
###BioSIM weather data




Weather<-ClimateHistorical









#ExportID<- subset(Plots [,1:3])
#length(ExportID[,1])##4882
#ExportID$Concatanate<-paste(ExportID$LATITUDE,ExportID$LONGITUDE, sep = "_")
#ExportID$LATITUDE <- as.numeric(ExportID$LATITUDE)
#ExportID$LONGITUDE <- as.numeric(ExportID$LONGITUDE)
#
#
#write.csv(ExportID, file="D:\\ModelOutput3\\ExportID.csv")
#str(ExportID)










Climate2011_2040<-cbind(MeanTemperature2011_2040,TemperatureMaxFebruary2011_2040,TemperatureMinAugust2011_2040,TemperatureMinJuneAugust2011_2040,SumPrecipSeptember2011_2040,SumPrecipMayAugust2011_2040)
str( Climate2011_2040)
Climate2011_2040[c(5,6,7,8,10,11,12,13,15,16,17,18,20,21,22,23,25,26,27,28,30)] <- list(NULL)
Years<-rep(c(2011:2040),30)
Climate2011_2040$Year<-Years
names(Climate2011_2040)
colnames(Climate2011_2040) <- c("Latitude","Longitude","Year","TMean", "TMaxFeb","TMinAug","TMinJuAu","PreciSept", "PreciMaySep")
head(Climate2011_2040)
Climate2011_2040$Concatanate<-paste(Climate2011_2040$Latitude,Climate2011_2040$Longitude,sep = "_")


ExportID <- subset(Plots [,1:3])
length(ExportID[,1])##4882
ExportID$Concatanate <- paste(ExportID$LATITUDE,ExportID$LONGITUDE, sep = "_")
ExportID$LATITUDE <- as.numeric(ExportID$LATITUDE)
ExportID$LONGITUDE <- as.numeric(ExportID$LONGITUDE)
ed <- ExportID[order(ExportID$ID_PEP_MES, ExportID$Concatanate, decreasing=TRUE),]
ed <- ed[!duplicated(ExportID$Concatanate),]
length(ed[,1])##3330
#ed$Status <- ifelse(duplicated(ed$Concatanate)=="FALSE",0,1)
write.csv(ed, file="D:\\ModelOutput3\\ExportID.csv")
str(ed)

ClimatePlots<- read.table("C:\\Users\\yomiq\\Documents\\YOSDATA\\LAVAL\\ModelOutput3\\InventoryPlotsCC.txt",colClasses = "character", header=TRUE,
                          sep = ",", quote="\"") ##contains info regarding historical and projected MFRI for rach plot
length(ClimatePlots[,1])

PoolPlots<-cbind(Plots,ClimatePlots)
length(PoolPlots[,1])#4882
str(PoolPlots)
PoolPlots[c(4,5,6,8,9,10,11,12,13,14,15,16,18,20,21,22,23,24,25,26,27,28,30,32,34)] <- list(NULL)
PoolPlots <- PoolPlots[which(PoolPlots$stands== "EnEn"),]
length(PoolPlots[,1])#3714
#PoolPlots$Pathway<-ifelse(PoolPlots$Historic_1=="783"&PoolPlots$MFRI1=="1112"&PoolPlots$MFRI2=="735"&PoolPlots$MFRI3=="404","Pathway1",
#ifelse(PoolPlots$Historic_1=="783"&PoolPlots$MFRI1=="1112"&PoolPlots$MFRI2=="735"&PoolPlots$MFRI3=="835","Pathway2",
#ifelse(PoolPlots$Historic_1=="783"&PoolPlots$MFRI1=="1112"&PoolPlots$MFRI2=="458"&PoolPlots$MFRI3=="835","Pathway3",
##ifelse(PoolPlots$Historic_1=="783"&PoolPlots$MFRI1=="1112"&PoolPlots$MFRI2=="735"&PoolPlots$MFRI3=="835","Others",
#ifelse(PoolPlots$Historic_1=="783"&PoolPlots$MFRI1=="776"&PoolPlots$MFRI2=="458"&PoolPlots$MFRI3=="300","Pathway4",
#ifelse(PoolPlots$Historic_1=="783"&PoolPlots$MFRI1=="776"&PoolPlots$MFRI2=="458"&PoolPlots$MFRI3=="135","Pathway5",
#ifelse(PoolPlots$Historic_1=="783"&PoolPlots$MFRI1=="776"&PoolPlots$MFRI2=="241"&PoolPlots$MFRI3=="135","Pathway6",
#ifelse(PoolPlots$Historic_1=="783"&PoolPlots$MFRI1=="776"&PoolPlots$MFRI2=="1251"&PoolPlots$MFRI3=="300","Pathway7",
#ifelse(PoolPlots$Historic_1=="783"&PoolPlots$MFRI1=="776"&PoolPlots$MFRI2=="1251"&PoolPlots$MFRI3=="835","Pathway8",
#ifelse(PoolPlots$Historic_1=="783"&PoolPlots$MFRI1=="831"&PoolPlots$MFRI2=="735"&PoolPlots$MFRI3=="404","Pathway9",
#ifelse(PoolPlots$Historic_1=="783"&PoolPlots$MFRI1=="831"&PoolPlots$MFRI2=="428"&PoolPlots$MFRI3=="145","Pathway10",
#ifelse(PoolPlots$Historic_1=="783"&PoolPlots$MFRI1=="304"&PoolPlots$MFRI2=="241"&PoolPlots$MFRI3=="146","Pathway11",
#ifelse (PoolPlots$Historic_1=="783"&PoolPlots$MFRI1=="581"&PoolPlots$MFRI2=="735"&PoolPlots$MFRI3=="404","Pathway12",
#ifelse(PoolPlots$Historic_1=="783"&PoolPlots$MFRI1=="581"&PoolPlots$MFRI2=="735"&PoolPlots$MFRI3=="145","Pathway13",
#ifelse(PoolPlots$Historic_1=="783"&PoolPlots$MFRI1=="182"&PoolPlots$MFRI2=="151"&PoolPlots$MFRI3=="300","Pathway14",
#ifelse(PoolPlots$Historic_1=="783"&PoolPlots$MFRI1=="182"&PoolPlots$MFRI2=="151"&PoolPlots$MFRI3=="145","Pathway15",
#ifelse(PoolPlots$Historic_1=="783"&PoolPlots$MFRI1=="182"&PoolPlots$MFRI2=="458"&PoolPlots$MFRI3=="404","Pathway16",
#ifelse(PoolPlots$Historic_1=="916"&PoolPlots$MFRI1=="776"&PoolPlots$MFRI2=="458"&PoolPlots$MFRI3=="300","Pathway17",
#ifelse(PoolPlots$Historic_1=="916"&PoolPlots$MFRI1=="776"&PoolPlots$MFRI2=="458"&PoolPlots$MFRI3=="835","Pathway18",
#ifelse(PoolPlots$Historic_1=="916"&PoolPlots$MFRI1=="776"&PoolPlots$MFRI2=="458"&PoolPlots$MFRI3=="390","Pathway19",
#ifelse(PoolPlots$Historic_1=="916"&PoolPlots$MFRI1=="776"&PoolPlots$MFRI2=="458"&PoolPlots$MFRI3=="229","Pathway20",
#ifelse(PoolPlots$Historic_1=="916"&PoolPlots$MFRI1=="776"&PoolPlots$MFRI2=="241"&PoolPlots$MFRI3=="300","Pathway21",
#ifelse(PoolPlots$Historic_1=="916"&PoolPlots$MFRI1=="304"&PoolPlots$MFRI2=="241"&PoolPlots$MFRI3=="146","Pathway22",
#ifelse(PoolPlots$Historic_1=="916"&PoolPlots$MFRI1=="304"&PoolPlots$MFRI2=="241"&PoolPlots$MFRI3=="135","Pathway23",
#ifelse(PoolPlots$Historic_1=="916"&PoolPlots$MFRI1=="304"&PoolPlots$MFRI2=="241"&PoolPlots$MFRI3=="229","Pathway24",
##ifelse(PoolPlots$Historic_1=="916"&PoolPlots$MFRI1=="304"&PoolPlots$MFRI2=="241"&PoolPlots$MFRI3=="146","Pathway25",
#ifelse(PoolPlots$Historic_1=="152"&PoolPlots$MFRI1=="182"&PoolPlots$MFRI2=="428"&PoolPlots$MFRI3=="145","Pathway25",
#ifelse(PoolPlots$Historic_1=="152"&PoolPlots$MFRI1=="182"&PoolPlots$MFRI2=="428"&PoolPlots$MFRI3=="404","Pathway26",
#ifelse(PoolPlots$Historic_1=="152"&PoolPlots$MFRI1=="182"&PoolPlots$MFRI2=="151"&PoolPlots$MFRI3=="145","Pathway27",
#ifelse(PoolPlots$Historic_1=="2865"&PoolPlots$MFRI1=="1112"&PoolPlots$MFRI2=="1251"&PoolPlots$MFRI3=="835","Pathway28",
#ifelse(PoolPlots$Historic_1=="3589"&PoolPlots$MFRI1=="938"&PoolPlots$MFRI2=="313"&PoolPlots$MFRI3=="229","Pathway29",
#ifelse(PoolPlots$Historic_1=="3589"&PoolPlots$MFRI1=="776"&PoolPlots$MFRI2=="458"&PoolPlots$MFRI3=="229","Pathway30",
#ifelse(PoolPlots$Historic_1=="3589"&PoolPlots$MFRI1=="304"&PoolPlots$MFRI2=="313"&PoolPlots$MFRI3=="229","Pathway31",
#"Others")))))))))))))))))))))))))))))))
#
#
#Grouping <- subset (PoolPlots [c(100,3713,3714,3592,750,2548,406,2550,3256,3253,2362,84,77,3679,2328,1,2251,419,2221,2228,2375,2266,2096,2260,19,31,74,2189,2240,2236,2241),7:11])
#str(Grouping)
#H<-as.numeric(Grouping[,1])
#m1<-as.numeric(Grouping[,2])
#m2 <- as.numeric(Grouping[,3])
#m3 <- as.numeric(Grouping[,4])
#Trend <- data.frame(Pathway= rep(1:31,4),MFRI=c(H,m1,m2,m3),Period= c(rep("1",31),rep("2",31),rep("3",31),rep("4",31)))
#Trend$Pathway<-as.numeric(Trend$Pathway)
#Trend$Period<-as.numeric(Trend$Period)
#str(Trend)
#
## convert factor to numeric for convenience
#npaths <- max(Trend$Pathway)
## get the range for the x and y axis
#xrange <- range(Trend$Period)
#yrange <- range(Trend$MFRI)
#
## set up the plot
#plot(xrange, yrange, type="n", xlab="Period",
#  	ylab="MFRI (years)" )
#colors <- rainbow(npaths)
#linetype <- c(1:npaths)
#
#
## add lines
#for (i in 1:npaths) {
#  paths <- subset(Trend, Pathway==i)
# lines(paths$Period, paths$MFRI, type="b", lwd=1.5,
#    lty=linetype[i], col=colors[i])
#}
#
## add a legend
#legend(xrange[1], yrange[2], 1:npaths, cex=0.8, col=colors,
#  	 lty=linetype, title="Change in MFRI")
#
#
#
#plot(MFRI~Period , data=Trend)
#
# dev.off()
#
#plot(yos[2,])
#c(Grouping[,1],Grouping[,2])
##Grouping <- cbind(PoolPlots[,7:11], PoolPlots[c(100,3713,3714,3592,750,2548,406,2550,3256,3253,2362,84,77,3679,2328,1,2251,419,2221,2228,2375,2266,2096,2260,19,31,74,2189,2240,2236,2241),])
#G<-t(Grouping)
#as.data.frame(G)
#plot(Gr[1])
#x.sub8 <- x.df[c(1, 3), 1:3]
#yos[2,]
#






dbhl<-c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29) ###dbh lower limit
dbhu<- c(3,5,7,9,11,13,15,17,19,21,23,25,27,29,31) ##dbh upper limit
dbhq<-sqrt((dbhu^3-dbhl^3)/((dbhu-dbhl)*3)) #assuming a uniform dbh distribution over interval
baq<-(dbhq/2)^2*pi/1e4
bal <- (dbhl/2)^2*pi/1e4
bau <- (dbhu/2)^2* pi/1e4


cali<-function(stand,Y,Recr){
  
  ###############################################################
  ###############################################################
  #Sampled_ca <- Region_ca (reg=c("A2","D4","B3","C3"), Calibration)
  #Sampled_ca <- Sampled_ca[sample(1:dim(Sampled_ca)[1], size=1, replace=F),]
  Sampled_ca <- PoolPlots[sample(1:dim(PoolPlots)[1], size=1, replace=F),]
  PlotID<-Sampled_ca[1,2]
  Tree.List_ca <- GetTrees_ca (Tree,Plots=Sampled_ca)
  as.numeric(Tree.List_ca$DBH)
  as.character(Tree.List_ca$ESSENCE)
  #######saplings only EnEn
  newTree.List_ca<-Tree.List_ca
  #newTree.List_ca<-subset(Tree.List_ca, DBH==2& ESSENCE== "EPN"|(Tree.List_ca$DBH>=4), select=c("ID_PEP_MES","ESSENCE","DBH"))
  length(1:dim(newTree.List_ca)[1])
  ####################
  newTree.List_ca<-as.data.frame(lapply(newTree.List_ca[,],function(x)rep(x,25)))
  range.DBH<-c(seq(1,30, by=2), 100)
  #Resume the results by class
  newTree.List_ca$DBH<-as.numeric(newTree.List_ca$DBH)
  stand<-table(cut(newTree.List_ca$DBH, breaks=range.DBH, labels=seq(1,15)))
  stand[1:4]<-stand[1:4]*10
  ###Partition the basal area of big trees >31 cm and add number of trees that the surplus of basal area represents
  
  basal_big_class<-0.0707905544
  BAB<- rep(0,100)
  TBA<- 3.142*(newTree.List_ca[newTree.List_ca[,3]>31,3]/200)^2
  BAB<-round(TBA/basal_big_class,digits=0)
  y<-sum(BAB)
  stand[15]<-stand[15]+y
  stand<-stand
  dbhl<-c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29)
  dbhu<- c(3,5,7,9,11,13,15,17,19,21,23,25,27,29,31)
  bal<-(dbhl/2)^2*pi/1e4
  bau<-(dbhu/2)^2*pi/1e4
  xu<-((dbhu/2)^4 - (dbhl/2)^4)/4
  xl<-((dbhu/2)^3 - (dbhl/2)^3)/3
  dbhq<-numeric(length(dbhu))
  dbhq<-(xu/xl)*2
  baq<-(dbhq/2)^2*pi/1e4
  
  ######Initialize variables###
  n<-length(stand)
  N1s<-rep(1,n)
  N0s<-rep(0,n)
  deaths<-numeric(n)
  growth<-numeric(n)
  BA<-numeric(Y)
  Recruits<-numeric(Y)
  Parcela<-matrix(0,Y,n,byrow=TRUE)
  Mortality<-matrix(0,Y,n,byrow=TRUE)
  Muertos<-matrix(0,Y,n,byrow=TRUE)
  Crecimiento<-matrix(0,Y,n,byrow=TRUE)
  Growth<-matrix(0,Y,n,byrow=TRUE)
  Transition<-matrix(0,Y,n,byrow=TRUE)
  DeltaN<-matrix(0,Y,n,byrow=TRUE)
  DeltaNT<-matrix(0,Y,n,byrow=TRUE)
  CR<-matrix(0,Y,n,byrow=TRUE)
  PCR<-numeric(n)
  PCR0<-numeric(n)
  CCR<-numeric(n)
  ShiftCR<-numeric(n)
  Size<-matrix(0,Y,n,byrow=TRUE)
  iw<-dbhu-dbhl ###class width
  MeanDBH1<-dbhq[1]
  m2Ha<-1e4
  shannon<-diversity(stand[5:15], index="shannon", MARGIN=1, base=exp(1))
  #Germinants_natural<-ifelse(shannon>=2.0,RegenIrregular,RegenRegular)
  Germinants_natural<-Recr
  SaplingSurvival<-0.98
  PCR<-crownratio(dbhq,sum(baq*stand))
  PCR1<-crownratio(dbhq[1],sum(baq*stand))
  PCR0<-PCR
  CCR<-PCR
  
  
  for (y in 1:Y){
    bay<-sum(stand*baq)
    parcela<-stand
    #GI <- GrowthIndex[y]
    
    Germinants_natural<-Recr
    stand[1] <- stand[1] + Germinants_natural
    Recruits[y]<-Germinants_natural
    CCRRecruits<-Recruits[y]*PCR[1]
    
    
    
    annual_diam_incre <- DiameterGrowthRate(stand,dbhq,baq,14)
    adi<-annual_diam_incre
    graduating<- 1/(iw/adi)
    growth<- rbinom(N1s, stand, graduating) #growth is the number of trees that graduate from age class n to n+1 (note that no tree can grow out of the biggest size class)
    stand<-stand-growth
    stand<-stand+c(0,growth[1:n-1])
    BAIncrement<-sum(growth*baq) ##patch BA increment
    CCRgrowth<-growth*CCR
    
    surviving<-mortality(stand,dbhq,baq)
    surviving[1:4]<-SaplingSurvival
    deaths<-rbinom(N1s,stand,1-surviving)
    stand<-stand-deaths
    DeltaBA<- BAIncrement-sum(deaths*baq)
    DeltaN[y,]<-deaths+growth
    
    
    DeltaNT[y,]<-deaths+growth+Germinants_natural
    
    
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
    
    Parcela[y,]<-parcela
    Mortality[y,]<-1-surviving
    Crecimiento[y,]<-growth
    Muertos[y,]<-deaths
    BA[y]<-bay
    Growth[y,]<-adi
    Transition[y,]<-graduating
    Size[y,]<-stand
    CR[y,]<-CCR
  }
  
  res<-list(stand=stand,Size=Size,Parcela=Parcela,Muertos=Muertos,Crecimiento=Crecimiento,Mortality=Mortality,Growth=Growth,
            Transition,BA=BA, DeltaN=DeltaN,Recruits=Recruits,CR=CR,DeltaNT=DeltaNT)
  return(res)
}

n.iter<-180 #plots to check
Y<-500
#Ba_s<-numeric(n.iter)
Ba_s<- matrix(0,n.iter,Y,byrow=T)
Size.list <- vector("list", n.iter) # create list
for(i in 1:n.iter){Size.list[[i]] <- matrix(0,Y,15,byrow=T)}
do.call(rbind, Size.list)
Basal<- numeric(n.iter)
for (i in 1:n.iter){
  recruit_calibra<-cali(stand,500,20) ## con500 recruits/ha
  Ba_s[i,]<-recruit_calibra$BA
  Basal[i]<-mean(recruit_calibra$BA)
  Size.list[[i]]<-recruit_calibra$Size
  print(i)
}


##Get initial values at year 300 for different FRI
Year<-300
plots<-1000
#Sc<-"Historical"
#Season 1 Spring #2 Summer
BasalArea<-Ba_s[,Year]
#Recruitment<-Recruits_s[,Year]
BiomassTurnover<-Turnover_s[1:plots,Year]
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

