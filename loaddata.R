##Get the data##
###TREES FOR PSP AND TSP
#######################################################################################################################################
####Do separately the tree and sapling PSP
TreePSP <- read.csv("Tree_list.csv",colClasses = "character" )
SaplingPSP <- read.csv("Yos_result2.csv",colClasses = "character", header=FALSE )
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
TreeT<- read.table("TemporalTrees.txt",
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

PoolPlots<-read.csv("InventoryPlotsPathway.csv",colClasses = "character", header=TRUE,sep = ",")
str(PoolPlots)
PoolPlots[,24:34] <- sapply(PoolPlots[,24:34],as.numeric)
length(PoolPlots[,1])  ##4882...3714
PoolPlots2<-PoolPlots[which(PoolPlots$stands== "EnEn"),]
length(1:dim(PoolPlots2)[1])##3714



###GET INTENSITY FILE FOR THE C-2 FUEL TYPE#####################################
Intensities<-read.csv("Sopfeu.csv",colClasses = "character", header=TRUE,sep = ",")
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
# 
# FirstPreMayAu<- read.table("FirstPreMayAug.txt",header = TRUE, sep = ",")
# FirstPreSep<- read.table("FirstPrecSept.txt",header = TRUE, sep = ",")
# FirstTmax<- read.table("FirstMAT.txt",header = TRUE, sep = ",")
# FirstTFeb<- read.csv("FirstTMaxFeb.txt",header = TRUE, sep = ",")
# FirstTAug<- read.csv("FirstTMinAug.txt",header = TRUE, sep = ",")
# FirsTJuAu<- read.csv("FirstTMinJuneAug.txt",header = TRUE, sep = ",")


PreMayAu<- read.csv("HistPreMayAug.csv",header = TRUE, sep = ",")
PreSep<- read.csv("HistoricalPrecSept.csv",header = TRUE, sep = ",")
Tmax<- read.csv("HistoricalMAT.csv",header = TRUE, sep = ",")
TFeb<- read.csv("HistTMaxFeb.csv",header = TRUE, sep = ",")
TAug<- read.csv("HistTMinAug.csv",header = TRUE, sep = ",")
TJuAu<- read.csv("HistTMinJuneAug.csv",header = TRUE, sep = ",")
Historical <- cbind(Tmax,PreSep,PreMayAu,TFeb,TAug,TJuAu)
str(Historical)
Historical[c(4,5,7,8,10,11,13,14,16,17)] <- list(NULL)
colnames(Historical) <- c("Latitude","Longitude","TMean","PreciSept","PreciMayAug","TMaxFeb","TMinAug","TMinJuAu")
head(Historical)
length(Historical[,1]) ##146460
Historical$Year <- rep("Historical",length(Historical[,1]))
#Historical$x <-rep(1:30,each=4882)
#Historical$y<-rep(1:30, 4882)
#Historical$Index<-paste(Historical$x,Historical$y, sep = "_")
Historical$Concatanate <- paste(Historical$Latitude,Historical$Longitude, sep = "_")

PreMayAuF<- read.csv("FirstPreMayAug.csv",header = TRUE, sep = ",")
PreSepF <- read.csv("FirstPrecSept.csv",header = TRUE, sep = ",")
TmaxF <- read.csv("FirstMAT.csv",header = TRUE, sep = ",")
TFebF <- read.csv("FirstTMaxFeb.csv",header = TRUE, sep = ",")
TAugF <- read.csv("FirstTMinAug.csv",header = TRUE, sep = ",")
TJuAuF <- read.csv("FirstTMinJuneAug.csv",header = TRUE, sep = ",")
First <- cbind(TmaxF,PreSepF,PreMayAuF,TFebF,TAugF,TJuAuF)
First[c(4,5,7,8,10,11,13,14,16,17)] <- list(NULL)
colnames(First) <- c("Latitude","Longitude","TMean","PreciSept","PreciMayAug","TMaxFeb","TMinAug","TMinJuAu")
head(First)
First$Year <- rep("First",length(First[,1]))
First$x <-rep(1:30,each=4882)
First$y<-rep(1:30, 4882)
First$Index<-paste(First$x,First$y, sep = "_")
First$Concatanate <- paste(First$Latitude,First$Longitude, sep = "_")

PreMayAuS <- read.csv("SecondPreMayAug.csv",header = TRUE, sep = ",")
PreSepS <- read.csv("SecondPrecSept.csv",header = TRUE, sep = ",")
TmaxS <- read.csv("SecondMAT.csv",header = TRUE, sep = ",")
TFebS <- read.csv("SecondTMaxFeb.csv",header = TRUE, sep = ",")
TAugS <- read.csv("SecondTMinAug.csv",header = TRUE, sep = ",")
TJuAuS <- read.csv("SecondTMinJuneAug.csv",header = TRUE, sep = ",")
Second <- cbind(TmaxS,PreSepS,PreMayAuS,TFebS,TAugS,TJuAuS)
Second[c(4,5,7,8,10,11,13,14,16,17)] <- list(NULL)
colnames(Second) <- c("Latitude","Longitude","TMean","PreciSept","PreciMayAug","TMaxFeb","TMinAug","TMinJuAu")
head(Second)
Second$x <-rep(1:30,each=4882)
Second$y<-rep(1:30, 4882)
Second$Year <- rep("Second",length(Second[,1]))
Second$Index<-paste(Second$x,Second$y, sep = "_")
Second$Concatanate <- paste(Second$Latitude,Second$Longitude, sep = "_")


PreMayAuT <- read.csv("ThirdPreMayAug.csv",header = TRUE, sep = ",")
PreSepT <- read.csv("ThirdPrecSept.csv",header = TRUE, sep = ",")
TmaxT <- read.csv("ThirdMAT.csv",header = TRUE, sep = ",")
TFebT <- read.csv("ThirdTMaxFeb.csv",header = TRUE, sep = ",")
TAugT <- read.csv("ThirdTMinAug.csv",header = TRUE, sep = ",")
TJuAuT <- read.csv("ThirdTMinJuneAug.csv",header = TRUE, sep = ",")


Third <- cbind(TmaxT,PreSepT,PreMayAuT,TFebT,TAugT,TJuAuT)
Third[c(4,5,7,8,10,11,13,14,16,17)] <- list(NULL)
colnames(Third) <- c("Latitude","Longitude","TMean","PreciSept","PreciMayAug","TMaxFeb","TMinAug","TMinJuAu")
head(Third)
Third$x <-rep(1:30,each=4882)
Third$y<-rep(1:30, 4882)
Third$Index<-paste(Third$x,Third$y, sep = "_")
Third$Concatanate <- paste(Third$Latitude,Third$Longitude, sep = "_")
Third$Year <- rep("Third",length(Third[,1]))

Weather <- rbind(Historical,First,Second,Third)

#DROUGHTCODE
DCH <- read.csv("DroughtCodehistoric.csv",header = TRUE, sep = ",")
DCH$Year <- rep("Historical",length(DCH[,1]))
DCF <- read.csv("DroughCode2011.csv",header = TRUE, sep = ",")
DCF$Year <- rep("First",length(DCF[,1]))
DCS <- read.csv("DroughtCode2040.csv",header = TRUE, sep = ",")
DCS$Year <- rep("Second",length(DCS[,1]))
DCT <- read.csv("droughtcode2071time.csv",header = TRUE, sep = ",")
DCT$Year <- rep("Third",length(DCT[,1]))

DroughtCode <- rbind(DCH,DCF,DCS,DCT) 
DroughtCode$Concatanate <- paste(DroughtCode$Latitude,DroughtCode$Longitude, sep = "_")
