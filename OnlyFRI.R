



exe <- function(stand,Y) {
  FRIs <- c(rep(916,30), rep(716,30), rep(458,30), rep(300,30)) ##Path1
  #FRIs <- c(rep(916,30), rep(304,30), rep(241,30), rep(170,30))  ##Path2
  #FRIs <- c(rep(783,30), rep(212,30), rep(250,30), rep(248,30))  ##Path3
  #FRIs <- c(rep(152,30), rep(182,30), rep(289,30), rep(145,30))  ##Path4
  #FRIs <- c(rep(783,30), rep(1112,30), rep(735,30), rep(404,30))  ##Path5 
  
  CPool <- c(4.9,1.5,24.9,9.0,9.1,35.1,1.6,3.4,95.0)* MgKg #Path1
  #CPool <- c(5.0,1.6,27.0,9.7,10.0,36.8,1.7,3.7,97.8)* MgKg ##Path2
  #CPool <- c(5.0,1.5,24.9,8.9,8.8,35.3,1.5,3.4,95.5)* MgKg ##Path3
  #CPool <- c(5.0,1.8,24.8,9.2,8.7,35.7,1.5,3.2,94.9)* MgKg ##Path4
  #CPool <- c(4.8,1.6,22.3,8.2,8.1,33.1,1.4,3.1,91.4)* MgKg ##Path5
  
 Sampled <- subpop (path=c("Pathway17","Pathway19","Pathway20","Pathway21"), d=PoolPlots2) ##path 1
  #Sampled <- subpop (path=c("Pathway22","Pathway23","Pathway24"), d=PoolPlots2) ##path 2
  #Sampled <- subpop (path=c("Pathway11","Pathway14","Pathway15","Pathway16"), d=PoolPlots2) ##path 3
  #Sampled <- subpop (path=c("Pathway26","Pathway28"), d=PoolPlots2) ##path 4
 #Sampled <- subpop (path=c("Pathway1"), d=PoolPlots2) ##path 5
  Sampled <- Sampled[sample(1:dim(Sampled)[1], size=1, replace=T),]
  PlotID <- Sampled[1,2]
  WeatherPlot <- Sampled[1,1]
  Latitude <- Sampled[1,4]
  Longitude <- Sampled[1,5]
  DroughtCodes <- GetDC(DroughtCode, WeatherPlot) 
  DroughtCode.List <- CorrectDupli(DroughtCodes)
  DCs <- as.numeric(DroughtCode.List[1:30,3]) ## Projected Drought Codes 1981-2100
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
  
  
  # Initialize object variables to save simulation results
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
  DroughtCode <- numeric(Y)
  BALost <- numeric(Y)
  Delta_BA <- numeric(Y)
  DeltaN <- matrix(0, Y, n, byrow = TRUE)
  CR <- matrix(0, Y, n, byrow = TRUE)
  Heights <- matrix(0, Y, n, byrow = TRUE)
  ShiftCrownratio <- matrix(0, Y, n, byrow = TRUE)
  ShiftHeights <- matrix(0, Y, n, byrow = TRUE)
  Latitudes <- numeric(Y)
  Longitudes <- numeric(Y)
  DroughtCode <- numeric(Y)
  EmpiricalTemperature <- numeric(Y)
  CCGrowthIndex <- numeric(Y)
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
  InitialCBiomass <- sum(BioMassCarbon%*%as.matrix(stand))  # initial biomass for the site
  ICB <- InitialCBiomass  # initial biomass for the site
  NetPrimaryProductivity <- numeric(Y)
  TotalLiveBiomass <- numeric(Y)
  AppliedDecayRates <- matrix(0, Y, 9, byrow = TRUE)
  NetEcosystemProduction <- numeric(Y)
  Rh <- numeric(Y)
  MFRI <- numeric(Y)
  CarbonCombusted <- numeric(Y)
  CarbonEmissions1 <- numeric(Y) 
  CarbonEmissions2 <- numeric(Y) 
  CarbonEmissions3 <- numeric(Y) 
  NetBiomeProduction <- numeric(Y) 
  FuelConsumed <- numeric(Y) 
  PlotIDs <- numeric(Y)
  

  # main loop
  for (y in 1:Y) {
    
    HeadIntensity <- sample(IntensitiesWeighted, size=1, replace=F)
    I <- as.numeric(HeadIntensity)
    Istart <- I
    DC <-  sample(DCs, size=1, replace=F)
    shannon <- diversity(stand[5:15], index = "shannon", MARGIN = 1, base = exp(1)) ##updated shannon
    Structure[y] <- shannon
    FRI <- FRIs[y]
    MFRI[y] <- FRI
    pBurn <- 1/FRI
    Firedeaths <- rep(0, 15)
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
    Latitudes [y] <- Latitude
    Longitudes [y] <- Longitude
    MAT <- 0.36
    EmpiricalTemperature[y] <- MAT
    AppDecayRates <- Decayrates(MAT)
    AppliedDecayRates[y, ] <- AppDecayRates
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
    
  
    CarbonPoolTransferMatrix <- matrix(c(
    Decayrate[1] * 0.83, (1-Decayrate[1]-0.08), 0, 0.08, 0, 0, Decayrate[1] * (1-0.83), 0, 0, 0,
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
    annual_diam_incre <- DiameterGrowthRate(stand,dbhq,baq,14)
    adi <- annual_diam_incre
    graduating <- 1 / (iw / adi) #Transition probabilities
    growth <- rbinom(N1s, stand, graduating) #stochastically
    stand <- stand - growth
    stand <- stand + c(0, growth[1:n - 1]) #after growth
    LastDCT <- growth[15]  # prevents loosing trees in last diameter class
    stand [15] <- stand[15] + LastDCT
    BAIncrement <- sum(growth * baq)    # patch BA increment due to growth
    CCRgrowth <- growth * CCR           # the trees that grow bring their crown ratio
    Heightgrowth <- growth * HeightUpdated
    
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
                              0, 0, 0, 0, (1 - 0.64), 0, 0, 0, 0, (0.64 * 0.5), 0, (0.64 * 0.5), 0, 0
    ), nrow = 5, ncol = 14, byrow = TRUE)
    
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
    
    
    # Mortality (senescence and fire)
    surviving <- mortality(stand, dbhq, baq)
    surviving[1:4] <- SaplingSurvival
    Senescencedeaths <- rbinom(N1s, stand, 1 - surviving)
    Senescence[y, ] <- Senescencedeaths
    stand <- stand-Senescencedeaths        # fire deaths are taken care in the fire module
    
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
      DroughtCode[y] <- DC
      BALost[y] <- severity
      C <- 1.185*exp(-4.252)*exp(0.671*log(FuelLoad))*exp(0.71*log(DC))
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
    DeltaBA<- BAIncrement-sum(deaths*baq) ##net basal area increment after mortality (fire and senescence) and growth
    DeltaStand <- sum(stand) 
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


n.iter <- 1000 #plots to check
Y <- 120
MFRI_s <- matrix(0,n.iter,Y,byrow=T)
#stand dynamics
Ba_s <- matrix(0,n.iter,Y,byrow=T)
Size_list <- vector("list", n.iter) # create list
PreFireStand_s <-vector("list", n.iter)
Transition_list <- vector("list", n.iter)
Structure_s <- matrix(0,n.iter,Y,byrow=T)
Recruits_s <- matrix(0,n.iter,Y,byrow=T)
Heights_list <- vector("list", n.iter)
CR_list <- vector("list", n.iter)
DiameterGrowth_list<- vector("list", n.iter)
DecayRates_list<-vector("list", n.iter)
DOM_Pool_list <- vector("list", n.iter)
DOM_Inputs_list <- vector("list", n.iter)
DOM_Flux_s <- matrix(0,n.iter,Y,byrow=T)
PrimaryProductivity_s <- matrix(0,n.iter,Y,byrow=T)
NetEcosystemProduction_s <- matrix(0,n.iter,Y,byrow=T)
Rh_s <- matrix(0,n.iter,Y,byrow=T)
NetBiomeProduction_s <- matrix(0,n.iter,Y,byrow=T)
Turnover_s <- matrix(0,n.iter,Y,byrow=T)
TotalLiveBiomass_s <- matrix(0,n.iter,Y,byrow=T)
EmpiricalTemperature_s <- matrix(0,n.iter,Y,byrow=T)
DroughtCodes_s <- matrix(0,n.iter,Y,byrow=T)
SnagCProduction_s <- matrix(0,n.iter,Y,byrow=T)
Severity_s <- matrix(0,n.iter,Y,byrow=T)
InitialIntensity_s <- matrix(0,n.iter,Y,byrow=T)
NF_s <- numeric(n.iter)
Emissions_s <- matrix(0,n.iter,Y,byrow=T)

for(i in 1:n.iter){
  Size_list[[i]] <- matrix(0,Y,15,byrow=T)
  DiameterGrowth_list[[i]]<- matrix(0,Y,15,byrow=T)
  PreFireStand_s[[i]]<-matrix(0,Y,15,byrow=T)
  Transition_list[[i]]<-matrix(0,Y,15,byrow=T)
  Heights_list[[i]]<-matrix(0,Y,15,byrow=T)
  CR_list[[i]]<-matrix(0,Y,15,byrow=T)
  DOM_Pool_list[[i]]<- matrix(0,Y,9,byrow=T)
  DOM_Inputs_list[[i]]<- matrix(0,Y,9,byrow=T)
  DecayRates_list[[i]]<- matrix(0,Y,9,byrow=T)
  
}

do.call(rbind, Size_list)
do.call(rbind, Transition_list)
do.call(rbind, DiameterGrowth_list)
do.call(rbind, Heights_list)
do.call(rbind, CR_list)
do.call(rbind, DOM_Pool_list)
do.call(rbind, DOM_Inputs_list)
do.call(rbind, DecayRates_list)

for (i in 1:n.iter){
  CarbonModel <- exe(stand,Y)
  MFRI_s [i,] <- CarbonModel$MFRI
  Ba_s[i,]<-CarbonModel$BA
  Structure_s[i,]<-CarbonModel$Structure
  Recruits_s[i,]<-CarbonModel$Recruits
  Size_list[[i]]<-CarbonModel$Size
  DiameterGrowth_list[[i]] <-CarbonModel$DiameterGrowth
  PreFireStand_s[[i]] <-CarbonModel$PreFireStand
  Transition_list[[i]]<-CarbonModel$Transition
  Heights_list[[i]]<-CarbonModel$Heights
  CR_list[[i]] <- CarbonModel$CR
  TotalLiveBiomass_s[i,] <- CarbonModel$TotalLiveBiomass
  EmpiricalTemperature_s[i,] <-CarbonModel$EmpiricalTemperature
  DroughtCodes_s[i,] <- CarbonModel$DroughtCode
  SnagCProduction_s[i,] <- CarbonModel$SnagCProduction
  PrimaryProductivity_s[i,] <- CarbonModel$NetPrimaryProductivity
  NetEcosystemProduction_s[i,] <- CarbonModel$NetEcosystemProduction
  #NetBiomeProduction_s[i,] <- CarbonModel$NetBiomeProduction
  Rh_s[i,] <- CarbonModel$Rh
  Turnover_s[i,] <- CarbonModel$Turnover
  DOM_Pool_list[[i]] <- CarbonModel$DOMC_Pool
  DecayRates_list[[i]] <- CarbonModel$AppliedDecayRates
  DOM_Flux_s[i,] <- CarbonModel$DOM_Flux
  DOM_Inputs_list[[i]] <- CarbonModel$DOM_Inputs
  Severity_s[i,] <- CarbonModel$BALost
  InitialIntensity_s[i,] <- CarbonModel$InitialIntensity
  Emissions_s[i,] <- CarbonModel$CarbonEmissions3
  NF_s[i] <- CarbonModel$NF
  print(i)
}


Year <- 1:120
plots <- n.iter


#Season 1 Spring #2 Summer
BasalArea <- Ba_s[, Year]
BiomassLiveCStock <- TotalLiveBiomass_s[1:plots, Year]
NPP <- PrimaryProductivity_s[1:plots, Year]
SoilRespiration <- Rh_s[1:plots, Year]
NEP <- NetEcosystemProduction_s[1:plots, Year]
NBP <-  NetBiomeProduction_s[1:plots, Year]
Snags <- sapply(DOM_Pool_list, function(m) m[1:120,1]) #
SnagBranch <- sapply(DOM_Pool_list, function(m) m[1:120,2])
AGMedium <- sapply(DOM_Pool_list, function(m) m[1:120,3])
AGfast <- sapply(DOM_Pool_list, function(m) m[1:120,4])
AGveryfast <- sapply(DOM_Pool_list, function(m) m[1:120,5])
AGslow <- sapply(DOM_Pool_list, function(m) m[1:120,6])
BGveryfast <- sapply(DOM_Pool_list, function(m) m[1:120,7])
BGfast <- sapply(DOM_Pool_list, function(m) m[1:120,8])
BGslow <- sapply(DOM_Pool_list, function(m) m[1:120,9])
SoilCStock1 <- Snags+SnagBranch+AGMedium+AGfast+AGveryfast+AGslow+BGveryfast+BGfast+BGslow
SoilCStock <- t(SoilCStock1)
EcosystemCStock <- BiomassLiveCStock + SoilCStock
MineralSoil <- t(BGveryfast+BGslow)
Organic <- t(AGveryfast+AGslow)
WoodyDebris <- t(Snags+SnagBranch+AGfast+AGMedium)
sn <- t(Snags)
sb <- t(SnagBranch)
am <-t(AGMedium)
af <- t(AGfast)
avf <- t(AGveryfast)
asl <- t(AGslow)
bgvf <- t(BGveryfast)
bgf <- t(BGfast)
bgs <- t(BGslow)

Path1onlyFRI <- list(NPP,NEP,NBP,SoilRespiration,BiomassLiveCStock,SoilCStock,EcosystemCStock,MineralSoil,Organic,WoodyDebris,
               sn,sb,am,af,avf,asl,bgvf,bgf,bgs)
Path2onlyFRI <- list(NPP,NEP,NBP,SoilRespiration,BiomassLiveCStock,SoilCStock,EcosystemCStock,MineralSoil,Organic,WoodyDebris,
                     sn,sb,am,af,avf,asl,bgvf,bgf,bgs)
Path3onlyFRI <- list(NPP,NEP,NBP,SoilRespiration,BiomassLiveCStock,SoilCStock,EcosystemCStock,MineralSoil,Organic,WoodyDebris,
                     sn,sb,am,af,avf,asl,bgvf,bgf,bgs)
Path4onlyFRI <- list(NPP,NEP,NBP,SoilRespiration,BiomassLiveCStock,SoilCStock,EcosystemCStock,MineralSoil,Organic,WoodyDebris,
                     sn,sb,am,af,avf,asl,bgvf,bgf,bgs)

Path5onlyFRI <- list(NPP,NEP,NBP,SoilRespiration,BiomassLiveCStock,SoilCStock,EcosystemCStock,MineralSoil,Organic,WoodyDebris,
                     sn,sb,am,af,avf,asl,bgvf,bgf,bgs)

save("Path1onlyFRI", file = "Path1onlyFRI.RData")
save("Path2onlyFRI", file = "Path2onlyFRI.RData")
save("Path3onlyFRI", file = "Path3onlyFRI.RData")
save("Path4onlyFRI", file = "Path4onlyFRI.RData")
save("Path5onlyFRI", file = "Path5onlyFRI.RData")
