n.iter <- 1000#plots to check
Y <- 120
Weather <- Weather
Latitude_s <- matrix(0,n.iter,Y,byrow=T)
Longitude_s <- matrix(0,n.iter,Y,byrow=T)
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
GrowthIndex_list <- vector ("list", n.iter)
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
  CarbonModel <- exe(stand,Y,Weather)
  Latitude_s [i,] <- CarbonModel$Latitudes
  Longitude_s [i,] <- CarbonModel$Longitudes
  MFRI_s [i,] <- CarbonModel$MFRI
  Ba_s[i,]<-CarbonModel$BA
  Structure_s[i,]<-CarbonModel$Structure
  Recruits_s[i,]<-CarbonModel$Recruits
  Size_list[[i]]<-CarbonModel$Size
  DiameterGrowth_list[[i]] <-CarbonModel$DiameterGrowth
  PreFireStand_s[[i]] <-CarbonModel$PreFireStand
  Transition_list[[i]]<-CarbonModel$Transition
  Heights_list[[i]]<-CarbonModel$Heights
  GrowthIndex_list [[i]] <- CarbonModel$CCGrowthIndex
  CR_list[[i]] <- CarbonModel$CR
  TotalLiveBiomass_s[i,] <- CarbonModel$TotalLiveBiomass
  EmpiricalTemperature_s[i,] <-CarbonModel$EmpiricalTemperature
  DroughtCodes_s[i,] <- CarbonModel$DroughtCode
  SnagCProduction_s[i,] <- CarbonModel$SnagCProduction
  PrimaryProductivity_s[i,] <- CarbonModel$NetPrimaryProductivity
  NetEcosystemProduction_s[i,] <- CarbonModel$NetEcosystemProduction
  NetBiomeProduction_s[i,] <- CarbonModel$NetBiomeProduction
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


Lat <- Latitude_s[,1]
Lon <- Longitude_s[,1]
NEP <- colMeans(NetEcosystemProduction_s) ##means over all replicates
NPP <- colMeans(PrimaryProductivity_s)
SoilResp <- colMeans(DOM_Flux_s)
NBP <- colMeans( NetBiomeProduction_s)
BiomassLiveCStock <- colMeans(TotalLiveBiomass_s)
BasalArea <- colMeans(Ba_s)
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
SoilCStock <- rowMeans(SoilCStock1)
EcosystemCStock <- BiomassLiveCStock + SoilCStock

MineralSoil <- BGveryfast+BGslow
Organic <- AGveryfast+AGslow
WoodyDebris <- Snags+SnagBranch+AGfast+AGMedium
Org <- rowMeans(Organic)
Min <- rowMeans(MineralSoil)
WD <- rowMeans(WoodyDebris)
sn <- rowMeans(Snags)
sb <-rowMeans(SnagBranch)
am <-rowMeans(AGMedium)
af <- rowMeans(AGfast)
avf <- rowMeans (AGveryfast)
asl <- rowMeans(AGslow)
bgvf <- rowMeans(BGveryfast)
bgf <- rowMeans(BGfast)
bgs <- rowMeans (BGslow)



intensities <- rep(0, 120)
for (i in 1:120){
  I <- InitialIntensity_s[,i]
  I2 <- mean(I[I!=0]) 
  I2 <- ifesle()
  intensities[i] <- I2
  intensities <- replace(intensities,is.na(intensities),0)
}


droughts <- rep(0, 120)
for (i in 1:120){
  D <- DroughtCodes_s[,i]
  D2 <- mean(D[D!=0]) 
  droughts[i] <- D2
  droughts<- replace(droughts,is.na(droughts),0)
}

seve <- rep(0, 120)
for (i in 1:120){
  S <- Severity_s[,i]
  S2 <- mean(S[S!=0]) 
  seve[i] <- S2
  seve <- replace(seve,is.na(seve),0)
}

emi <- rep(0, 120)
for (i in 1:120){
  E <- Emissions_s[,i]
  E2 <- mean(E[E!=0]) 
  emi[i] <- E2
  emi <- replace(emi,is.na(emi),0)
}

GIs <- rowMeans(sapply(GrowthIndex_list, function(m) m[]))
AnnTemp <- colMeans(EmpiricalTemperature_s)




##DecompositionratesCarbonPools
SnagsDecorate <- rowMeans(sapply(DecayRates_list, function(m) m[1:120,1]))
SnagsBranchDecorate <- rowMeans(sapply(DecayRates_list, function(m) m[1:120,2]))
MediumDecRate <- rowMeans(sapply(DecayRates_list, function(m) m[1:120,3]))
agfastDecRate <- rowMeans(sapply(DecayRates_list, function(m) m[1:120,4]))
agvfDecRate <- rowMeans(sapply(DecayRates_list, function(m) m[1:120,5]))
agsDecRate <- rowMeans(sapply(DecayRates_list, function(m) m[1:120,6]))
bgvfDecRate <- rowMeans(sapply(DecayRates_list, function(m) m[1:120,7]))
bgfDecRate <- rowMeans(sapply(DecayRates_list, function(m) m[1:120,8]))
bgsDecRate <- rowMeans(sapply(DecayRates_list, function(m) m[1:120,9]))



fireregime <- cbind(intensities,droughts,seve,emi)
goann <- cbind(GIs,AnnTemp)
ApDecRates <- cbind(SnagsDecorate,SnagsBranchDecorate,MediumDecRate,agfastDecRate,
                    agvfDecRate,bgvfDecRate,bgfDecRate,bgsDecRate)
Path5 <- cbind(NPP,NEP,NBP,SoilResp,BiomassLiveCStock,SoilCStock,EcosystemCStock,Min,Org,WD,
               sn,sb,am,af,avf,asl,bgvf,bgf,bgs)

Path5Coor <- cbind(Lat,Lon)


save("ApDecRates", file="DecayRatePath5.RData")
save("fireregime", file="FireregimePath5.RData")
save("goann", file="GrowthTemp5.RData")
save("Path5", file="Path5.RData")
save("Path5Coor", file="Path5Coor.RData")






