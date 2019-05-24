###########################
## VALIDATION PROCEDURE ##
###########################

####Function that allows getting tree data from a particular sample plot given the inventory year
GetTrees_v <- function(T, i, invyear){
  res <- subset(T, T$plotid==i & T$Mes==invyear)
  return(res)
}

####Function that generates vector of number of tres by diameter class. Evaluates basal area, drainage class and deposit
manipulate <- function(x){
  x$DBH <- x$DHPMM/10
  as.character(c(x$DBH,x$ETAT,x$NO_ARBRE))
  as.factor(x$ESSENCE)
  newdata <- x[which(x$ETAT=="10"),]
  drainage <- newdata$CL_DRAI[1]
  deposit <- newdata$DEP_SUR [1]
  newdata$ba <- (newdata$DBH/2)^2 * pi/1e4
  baEPN<- newdata[which(newdata$ESSENCE=="EPN"),]
  pEPN <- sum(baEPN$ba)/sum(newdata$ba)
  newdata_ha<- as.data.frame(lapply(newdata[,],function(x)rep(x,25)))
  range.DBH<-c(seq(1,30, by=2), 100)
  #Resume the results by class
  newdata_ha$DBH <- as.numeric(newdata_ha$DBH)
  stand <- table(cut(newdata_ha$DBH, breaks=range.DBH, labels=seq(1,15)))
  res <- list(stand=stand,pEPN=pEPN,drainage=drainage,deposit=deposit)
  return(res)
}

manipulate_vm <- function(x){
  x$DBH <- x$DHPMM/10
  as.character(c(x$DBH,x$ETAT,x$NO_ARBRE))
  as.factor(x$ESSENCE)
  newdata <- x[which(x$ETAT=="10"),]
  newdata_ha<- as.data.frame(lapply(newdata[,],function(x)rep(x,25)))
  range.DBH<-c(seq(1,30, by=2), 100)
  #Resume the results by class
  newdata_ha$DBH <- as.numeric(newdata_ha$DBH)
  stand <- table(cut(newdata_ha$DBH, breaks=range.DBH, labels=seq(1,15)))
  res <- list(stand=stand)
  return(res)
}


######################
###Validation plots###
######################
plots <- read.csv("Plots_infoinventaire.csv", header=T) ##all plots
trees <- read.csv("Arboles.csv", header=T)## all trees
plot_domain6 <- read.csv("plots_zone6.csv", header=T) ##plots domain 6
plots$ID_PEP 
plots$ID_PEP <- as.character(plots$ID_PEP) 
plots$VERSION <- as.factor(plots$VERSION)
plotsier <- plots[which(plots$VERSION=="1er inventaire 1970 à 1974"),]
plotsdeux <-plots[which(plots$VERSION=="2e inventaire" ),]
plotstroi <-plots[which(plots$VERSION=="3e inventaire"),]

plotsier$ID_PEP <- as.character(plotsier$ID_PEP)
plotsdeux$ID_PEP <- as.character(plotsdeux$ID_PEP)
plotstroi$ID_PEP <- as.character(plotstroi$ID_PEP)


trees$ID_PEP_MES.1 <- as.character(trees$ID_PEP_MES.1)
splits <- strsplit(trees$ID_PEP_MES.1,"")
splits2 <- matrix(0, nrow=length(trees[,1]), ncol=1) 
splits3 <- matrix(0, nrow=length(trees[,1]), ncol=1) 
for (i in 1:length(trees[,1])){
  tmp <- splits[[i]][11:12]
  tmp2 <- splits[[i]][1:10]
  splits2[i] <- paste(tmp[1],tmp[2], sep="")
  splits3[i] <- paste(tmp2[1],tmp2[2],tmp2[3],
                      tmp2[4],tmp2[5],tmp2[6],
                      tmp2[7],tmp2[8],tmp2[9],
                      tmp2[10],sep="")
}


trees$Mes <- splits2
trees$plotid <- splits3
str(trees)
trees$Mes<- as.factor(trees$Mes)


plot_domain6$ID_PEP_MES.1 <- as.character(plot_domain6$ID_PEP_MES) 
splits_plots <- strsplit(plot_domain6$ID_PEP_MES.1,"")

splits2_plots <- matrix(0, nrow=length(plot_domain6[,1]), ncol=1) 
splits3_plots <- matrix(0, nrow=length(plot_domain6[,1]), ncol=1) 
for (i in 1:length(plot_domain6[,1])){
  tmp <- splits_plots[[i]][11:12]
  tmp2 <- splits_plots[[i]][1:10]
  splits2_plots[i] <- paste(tmp[1],tmp[2], sep="")
  splits3_plots[i] <- paste(tmp2[1],tmp2[2],tmp2[3],
                      tmp2[4],tmp2[5],tmp2[6],
                      tmp2[7],tmp2[8],tmp2[9],
                      tmp2[10],sep="")
}

plot_domain6$ID_PEP <- splits3_plots

##Get only plots in bioclimatic domain 6, not used for calibration or evaluation
## plots in bioclimatic domain 6 and not found in evaluation plots
inventory_1_2_zone6 <-  plot_domain6[which(plot_domain6$ID_PEP%in%inventory_1_2$ID_PEP),] ##913
inventory_2_3_zone6 <-  plot_domain6[which(plot_domain6$ID_PEP%in%inventory_2_3$ID_PEP),] ##5644

plotsmodel <- subpop (path=c("Pathway17","Pathway19","Pathway20","Pathway21","Pathway22","Pathway23","Pathway24","Pathway11",
                             "Pathway14","Pathway15","Pathway16",
                             "Pathway26","Pathway28","Pathway1"), d=PoolPlots2) ##plots used for model evaluation
inventory_2_3_zone6v <- inventory_2_3_zone6[which(!inventory_2_3_zone6$ID_PEP_MES%in%plotsmodel$ID_PEP_MES),]##3310 making sure I do not use the same plots for model evaluation
inventory_1_2_zone6v <- inventory_1_2_zone6[which(!inventory_1_2_zone6$ID_PEP_MES%in%plotsmodel$ID_PEP_MES),]##863

##PLOTS TO USE
info3depurados$invplotID_s # plots in bioclimatic zone 6, EnEn, mesic till sites
test1 <- plotsier[which(plotsier$ID_PEP%in%as.character(info3depurados$invplotID_s)),] 
test2 <- plotsdeux[which(plotsdeux$ID_PEP%in%as.character(info3depurados$invplotID_s)),]
test3 <- plotstroi[which(plotstroi$ID_PEP%in%as.character(info3depurados$invplotID_s)),]

#cleaned plots found in both inventory 3 and 2
inventory_2_3 <- test2 [which(unique(test2$ID_PEP)%in%unique(test3$ID_PEP)),] ##134 plots

#cleaned plots found in both inventory 2 and 1
inventory_1_2 <- test2 [which(unique(test2$ID_PEP)%in%unique(test1$ID_PEP)),]## 20 plots

###loop to get change in biomass among inventory plots
  n <- rep(0,15)
  N1s <- rep(1, 15)
  N0s <- rep(0, 15)
  dbhl  <- c(1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29)  # dbh lower limit
  dbhu <- c(3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31)  # dbh upper limit
  dbhq  <- sqrt((dbhu^3-dbhl^3)/((dbhu-dbhl) * 3))  # assuming a uniform dbh distribution over interval
  baq  <- (dbhq/2)^2 * pi/1e4
  # Initialise variables
  hatom2 <- 1e4
  iw <- dbhu-dbhl
  m2Ha  <- 1e4
  biocar_factor <- 0.5
  MgKg <- 1000
  BioMass <- matrix(c(Stemwood(dbhq), Bark(dbhq), Branches(dbhq), Needles(dbhq), Coarse(dbhq), Fineroots(dbhq)), nrow = 6,
                    ncol = length(dbhq), byrow = TRUE)
  BioMassCarbon <- BioMass * biocar_factor  # Biomass C per diameter class
  
  ###Initialise object variables to save simulation results####
  Y=134 ##corresponds to the number of plots in bioclimatic domain 6, mesi-till sites
  IAB3_s <- numeric(Y) 
  IAB2_s <- numeric(Y)
  IAB1_s <- numeric(Y)
  invplotID_s <- numeric(Y)
  EnEn3_s <- numeric(Y)
  EnEn2_s <- numeric(Y)
  EnEn1_s <- numeric(Y)
  Dr3_s <- numeric(Y)
  Dr2_s <- numeric(Y)
  Dr1_s <- numeric(Y)
  De3_s <- numeric(Y)
  De2_s <- numeric(Y)
  De1_s <- numeric(Y)
  
  for (y in 1:Y){ #Things that have to be reinitialized and updated
    i <- sample(unique(inventory_2_3$ID_PEP),1, replace=F)
    i3 <- GetTrees_v(trees,i,"03")# 3er inventaire
    i2 <- GetTrees_v(trees,i,"02")# 2ieme inventaire
    i1 <- GetTrees_v(trees,i,"01")# 1erie inventaire
    data3 <- manipulate(i3)  
    data2 <- manipulate(i2) 
    data1 <- manipulate(i1)
    stand3 <- data3$stand
    stand2 <- data2$stand
    stand1 <- data1$stand
    stand3EnEn <- data3$pEPN
    stand2EnEn <- data2$pEPN
    stand1EnEn <- data1$pEPN
    
    stand3deposit <-  if(is.na(data3$deposit)== TRUE) {
      data3$deposit<- NA
    } else{
      data3$deposit
    }
                        
    stand2deposit<- if(is.na(data2$deposit)== TRUE) {
      data2$deposit <- NA
    } else{
      data2$deposit
    }
    
    stand1deposit<- if(is.na(data1$deposit)== TRUE) {
      data1$deposit <- NA
    } else{
      data1$deposit
    }
    
    stand3drainage <-  if(is.na(data3$drainage)== TRUE) {
      data3$drainage <- 0
    } else{
      data3$drainage
    }
    
    stand2drainage<- if(is.na(data2$drainage)== TRUE) {
      data2$drainage <- 0
    } else{
      data2$drainage
    }
    
    stand1drainage<- if(is.na(data1$drainage)== TRUE) {
      data1$drainage <- 0
    } else{
      data1$drainage
    }
    
  
    CBiomass_3 <- sum(BioMassCarbon[1:4,5:15]%*%as.matrix(stand3[5:15]))
    CBiomass_2 <- sum(BioMassCarbon[1:4,5:15]%*%as.matrix(stand2[5:15])) 
    CBiomass_1 <- sum(BioMassCarbon[1:4,5:15]%*%as.matrix(stand1[5:15]))
    
    IAB3_s[y] <-  CBiomass_3
    IAB2_s[y] <-  CBiomass_2
    IAB1_s[y] <-  CBiomass_1
    invplotID_s[y] <- i
    EnEn3_s[y] <- stand3EnEn
    EnEn3_s[y] <- stand2EnEn 
    EnEn1_s[y] <- stand1EnEn
    Dr3_s[y] <- stand3drainage
    Dr2_s[y] <- stand2drainage
    Dr1_s[y] <- stand1drainage
    De3_s[y] <- as.character(stand3deposit) 
    De2_s[y] <- as.character(stand2deposit) 
    De1_s[y] <- as.character(stand1deposit) 
  } 
 
  



########################################################################
# Calculate change in annual AG biomass bewteen inventory 1,2 and 3

info3 <- as.data.frame(cbind(invplotID_s,De3_s,Dr3_s,EnEn3_s,IAB3_s))
info2 <- as.data.frame(cbind(invplotID_s,De2_s,Dr2_s,EnEn2_s,IAB2_s))
info1 <- as.data.frame(cbind(invplotID_s,De1_s,Dr1_s,EnEn1_s,IAB1_s))

###select plots with tree information, EnEn, drainage and deposit classes as in the model
info3$EnEn <- as.numeric(as.character(info3$EnEn3_s)) 
info3$DrCl<- as.character(info3$Dr3_s)
info3$Dep <- as.character(info3$De3_s)
info3$B3 <- as.numeric(as.character(info3$IAB3_s))/MgKg
info2$B2 <- as.numeric(as.character(info2$IAB2_s))/MgKg
info1$B1 <- as.numeric(as.character(info1$IAB1_s))/MgKg
##########################################################
# these plots are the selected based on EnEn >0.75, till mesic sites and used again to subset plots above...
info3depurados  <- subset(info3, (EnEn  >=0.75)
                          & (DrCl=="20"|DrCl=="21"|DrCl=="30"|
                               DrCl=="31"|DrCl=="34"|DrCl=="40"|   
                               DrCl=="41"|DrCl=="42") 
                          &(Dep =="1A"|Dep =="1"|Dep =="1AY"|
                              Dep =="1AM"|Dep =="M1A"|Dep =="2"),
                          select = c(1:length(info3)))

info3depurados$B3 <- as.numeric(as.character(info3depurados$IAB3_s))                          
length(info3depurados[,1]) # 406 plots

######################################
##Empirical change in annual biomass
#####################################
Biomasschange1_2 <- cbind(info2,info1)
Biomasschange2_3 <- cbind(info3,info2)
#Biomass change inventory 1 and 2
Biomasschange1_2$CB <- Biomasschange1_2$B2-Biomasschange1_2$B1
Biomasschange1_2$ACB <- Biomasschange1_2$CB/10
Biomasschange1_2nonnegative <- subset(Biomasschange1_2,(Biomasschange1_2$ACB>0),select = c(1:length(Biomasschange1_2)))
empirical1_2 <- (Biomasschange1_2nonnegative$ACB)###0.34
#Biomass change inventory 2 and 3
Biomasschange2_3$CB <- Biomasschange2_3$B3-Biomasschange2_3$B2
Biomasschange2_3$ACB <- Biomasschange2_3$CB/10
Biomasschange2_3nonnegative <- subset(Biomasschange2_3,(Biomasschange2_3$ACB>0),select = c(1:length(Biomasschange2_3)))##only get plots where there was a positive change in biomass
empirical2_3 <- (Biomasschange2_3nonnegative$ACB)###0.32

## same plots used for validation
invempirical1_2 <-  as.character(Biomasschange1_2nonnegative$invplotID_s)
invempirical2_3 <-  as.character(Biomasschange2_3nonnegative$invplotID_s)

#######################################
#Modelled annual changes in biomass
#######################################

exe <- function(Y){
  Sampled <- sample(Biomasschange2_3nonnegative$invplotID_s,1,replace=F) # 110 
  i <- as.character(Sampled)
  Tree.List <- GetTrees_v(trees,i,"01")#first inventary
  data_trees <- manipulate_vm(Tree.List)
  PlotID <- i
  WeatherPlot <- as.character(i)
  GIndex <- 14 ## mean GI for historical period
  stand <- data_trees$stand
  n <- length(stand)
  N1s <- rep(1, n)
  N0s <- rep(0, n)
  dbhl  <- c(1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29)  # dbh lower limit
  dbhu <- c(3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31)  # dbh upper limit
  dbhq  <- sqrt((dbhu^3-dbhl^3)/((dbhu-dbhl) * 3))  # assuming a uniform dbh distribution over interval
  baq  <- (dbhq/2)^2 * pi/1e4
  # Initialise variables
  hatom2 <- 1e4
  iw <- dbhu-dbhl
  m2Ha  <- 1e4
  biocar_factor <- 0.5
  MgKg <- 1000
  SaplingSurvival <- 0.98 
  NF <- 0
  fire.year <- NULL
  RegenLagTime <- 27   # assume it takes a 3-yr old natural seedling 27 yearss to grow to
  # 214 cm (DBH 1.0 cm), VanBoagaert et al. 2015)
  RegenRegular <- 0  # natural regeneration
  RegenIrregular <- 0 # natural regeneration
  RegenCohorts <- rpois(RegenLagTime, ifelse(shannon<1.7, RegenRegular, RegenIrregular))
  
  ###Initialise object variables to save simulation results####
  PlotIDs <- numeric(Y)
  IAB_s <- numeric(Y) ##initial C biomass (for validation purposes)
  FAB_s <- numeric(Y) 
  
  BioMass <- matrix(c(Stemwood(dbhq), Bark(dbhq), Branches(dbhq), Needles(dbhq), Coarse(dbhq), Fineroots(dbhq)), nrow = 6,
                    ncol = length(dbhq), byrow = TRUE)
  BioMassCarbon <- BioMass * biocar_factor  # Biomass C per diameter class
  
  
  InitialCBiomass <- sum(BioMassCarbon[1:4,5:15]%*%as.matrix(stand[5:15])) ##bio1 for the site
  ICB <- InitialCBiomass
  
  
  for (y in 1:Y){ #Things that have to be reinitialized and updated
    
    # Stand dynamics
    annual_diam_incre <- DiameterGrowthRate(stand,dbhq,baq,GIndex)
    adi<-annual_diam_incre
    graduating <- 1 / (iw / adi) #Transition probabilities
    growth <- rbinom(N1s, stand, graduating) #stochastically
    stand <- stand- growth
    stand <- stand + c(0, growth[1:n - 1]) #after growth
    LastDCT <- growth[15]  # prevents loosing trees in last diameter class
    stand [15] <- stand[15] + LastDCT
    BAIncrement <- sum(growth*baq)    # patch BA increment due to growth
    GrowthCBiomass <- sum(BioMassCarbon[1:4,5:15]%*%as.matrix(stand[5:15])) ##calculate biomass due to growth..Biomass that has not been lost to turnover or mortality
    delta <- (GrowthCBiomass-ICB) # Biomass that has not been lost to turnover or mortality
    
    # Mortality (senescence)
    surviving <- mortality(stand, dbhq, baq)
    Senescencedeaths <- rbinom(N1s, stand, 1 - surviving)
    stand <- stand-Senescencedeaths # fire deaths are taken care in the fire module
    if(sum(stand<0)>0)
      stop(message="neg count 3")
    
    # Add recruitment
    
    Recruits[y] <- RegenCohorts[RegenLagTime]
    RegenCohorts <- c(NewRegen, RegenCohorts[1:RegenLagTime - 1])
    stand[1] <- stand[1] + Recruits[y]
    
    # Update biomass
    FinalAbovegroundBiomass <- sum(BioMassCarbon[1:4,5:15] %*% as.matrix(stand[5:15]))
    ICB <- FinalAbovegroundBiomass
    Latitudes[y]<-  Latitude
    Longitudes[y]<- Longitude
    PlotIDs[y]<- WeatherPlot 
    IAB_s[y] <- InitialCBiomass #initial aboveground biomass
    FAB_s[y] <- FinalAbovegroundBiomass #final aboveground biomass
  }
  
  res <- list(PlotIDs=PlotIDs,
              IAB_s=IAB_s,FAB_s=FAB_s)
  
  return(res)
}


#run simulations
n.iter <- 110 #plots to check
Y <- 10
IAB_s_s <- matrix(0,n.iter,Y,byrow=T)
FAB_s_s <- matrix(0,n.iter,Y,byrow=T)
PlotIDs_s <-  matrix(0,n.iter,Y,byrow=T) 

for (i in 1:n.iter){
  CarbonModel <- exe(Y)
  IAB_s_s[i,] <- CarbonModel$IAB_s
  FAB_s_s[i,] <- CarbonModel$FAB_s
  PlotIDs_s[i,] <-  CarbonModel$PlotIDs
  print(i)
}

Year <- 10
plots <- 110
InitialAbovegBiomass <- IAB_s_s[1:plots,Year]
FinalAbovegBiomass <- FAB_s_s[1:plots,Year]
invplotID_s <- PlotIDs_s[1:plots,Year]
################
modelleddata <- as.data.frame(cbind(invplotID_s,InitialAbovegBiomass,FinalAbovegBiomass))
as.numeric(as.character(modelleddata$FinalAbovegBiomass))
modelleddata$AGmodelled <- as.numeric(as.character(modelleddata$FinalAbovegBiomass))-as.numeric(as.character(modelleddata$InitialAbovegBiomass))
modelleddata$ACB <- (modelleddata$AGmodelled/10)/MgKg
mean(modelleddata$ACB[modelleddata$ACB>0])#0.23
###############
modelleddata1 <- as.data.frame(cbind(invplotID_s,InitialAbovegBiomass,FinalAbovegBiomass))
as.numeric(as.character(modelleddata1$FinalAbovegBiomass))
modelleddata1$AGmodelled <- as.numeric(as.character(modelleddata1$FinalAbovegBiomass))-as.numeric(as.character(modelleddata1$InitialAbovegBiomass))
modelleddata1$ACB <- (modelleddata1$AGmodelled/10)/MgKg
mean(modelleddata1$ACB[modelleddata1$ACB>0]) #0.2184878

setwd("~/Desktop/Validation")
save("modelleddata1",file="modelledachange_2_3.RData")
save("modelleddata",file="modelledachange_1_2.RData")
save("Biomasschange1_2nonnegative",file="empiricalachange1_2.RData")
save("Biomasschange2_3nonnegative",file="empiricalachange2_3.RData")

####################################
#Validation procedure
#Analysis and plots
####################################
### Cumulative density functions
# create ECDF of data
cdf_inven <- ecdf(sort(Biomasschange1_2nonnegative$ACB))
cdf_modelled <-ecdf(sort(modelleddata1$ACB))
# cdf_invent1 <- ecdf(sort(annualchangebiomass12_normalized))
plot(cdf_inven, verticals=TRUE, do.points=FALSE, col="grey50",lwd=2, main=NULL,ylab="Cumulative distribution function", xlab="Annual change in AG biomass") 
plot(cdf_modelled  , verticals=TRUE, do.points=FALSE, col="black",lwd=2, add=TRUE) 
legend(0.4, 0.4, legend=c("Empirical", "Modelled"), col=c("grey50","black"), lty=1, bty = "n")
####################
cdf_modelled(c(0.3,0.7))
cdf_inven(c(0.3,0.7))
####################
#KS test
####################
my.ecdf.i <- cdf_inven(sort(Biomasschange1_2nonnegative$ACB))
my.ecdf.m <- cdf_inven(sort(modelleddata1$ACB))
ks.test(my.ecdf.i,my.ecdf.m)
###Boxplots
par(mfrow = c(1, 2))
boxplot(Biomasschange1_2nonnegative$ACB,ylab ="Annual change in AG biomass (MgC/ha)", main="Empirical", ylim=c(0,0.9))
boxplot(modelleddata1$ACB,ylab ="Annual change in AG biomass (MgC/ha)", main="Modelled", ylim=c(0,0.9))
t.test(empirical1_2, modelleddata1$ACB, var.equal = FALSE)
###Histograms
par(mfrow = c(1, 2))
hist(Biomasschange1_2nonnegative$ACB)
hist(modelleddata1$ACB)




