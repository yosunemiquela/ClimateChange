
library("reshape")
library("ggplot2")
library("vegan")
library("plotrix")
library("relaimpo")
library("MASS")
library("car")
library("multcomp")
library("class")
library("Hmisc")
library("zoo")
library("scatterplot3d")
library("fitdistrplus")
library("compare")
library ("truncreg")
library("truncdist")
library("msm") ## for dtnorm
library("actuar")
library("moments")
library("stats")
library("plyr")
library("quantmod")
library("nlme")
library("grid")


KgMg<-(0.001)

load("Path1.RData")
load("Path2.RData")
load("Path3.RData")
load("Path4.RData")
load("Path5.RData")

##get at year 30, 60, 90 120
pathfunction <- function(x,i){
  third <- x[[i]][,120]
  second <- x[[i]][,90] 
  first <- x[[i]][,60] 
  hist <- x[[i]][,30] 
  init<- x[[i]][,1] 
  data <- as.data.frame(cbind(init,hist,first,second,third))
  Pathway <- melt(data)
  Pathway<-Pathway[,2]
  return(Pathway)
}


pathfunction2 <- function(x,i){
  Df <- x[[i]][,120]
  Di <- x[[i]][,91]
  Cf <- x[[i]][,90] 
  Ci <- x[[i]][,61]
  Bf <- x[[i]][,60] 
  Bi<- x[[i]][,31]
  Af<- x[[i]][,30] 
  Ai<- x[[i]][,1] 
  data <- as.data.frame(cbind(Ai,Af,Bi, Bf, Ci, Cf, Di, Df))
  Pathway <- melt(data)
  Pathway<-Pathway[,2]
  return(Pathway)
}

PTWAY1_ <- matrix(0,8000,19,byrow=T)
PTWAY2_ <- matrix(0,8000,19,byrow=T)
PTWAY3_ <- matrix(0,8000,19,byrow=T)
PTWAY4_ <- matrix(0,8000,19,byrow=T)
PTWAY5_ <- matrix(0,8000,19,byrow=T)

for (i in 1:19){
  
  
  PTWAY1_[,i] <- pathfunction2(Path1,i)
  PTWAY2_[,i] <- pathfunction2(Path2,i)
  PTWAY3_[,i] <- pathfunction2(Path3,i)
  PTWAY4_[,i] <- pathfunction2(Path4,i)
  PTWAY5_[,i] <- pathfunction2(Path5,i)
}

PTWAY1a_ <- as.data.frame(PTWAY1_)
PTWAY1a_$Year<-rep(c("Ai","Af","Bi","Bf","Ci","Cf","Di","Df"),each=1000)
PTWAY1a_$Pathway<-rep("Path1", 1000)
colnames(PTWAY1a_)<- c("NPP","NEP","NBP","SoilRespiration",
                      "BiomassLiveCStock","SoilCStock","EcosystemCStock",
                      "MineralSoil","Organic","WoodyDebris",
                      "sn","sb","am","af","avf","asl","bgvf","bgf","bgs","Period","Pathway")

PTWAY2a_ <- as.data.frame(PTWAY2_)
PTWAY2a_$Year<-rep(c("Ai","Af","Bi","Bf","Ci","Cf","Di","Df"),each=1000)
PTWAY2a_$Pathway<-rep("Path2", 1000)
colnames(PTWAY2a_)<- c("NPP","NEP","NBP","SoilRespiration",
                      "BiomassLiveCStock","SoilCStock","EcosystemCStock",
                      "MineralSoil","Organic","WoodyDebris",
                      "sn","sb","am","af","avf","asl","bgvf","bgf","bgs","Period","Pathway")

PTWAY3a_ <- as.data.frame(PTWAY3_)
PTWAY3a_$Year<-rep(c("Ai","Af","Bi","Bf","Ci","Cf","Di","Df"),each=1000)
PTWAY3a_$Pathway<-rep("Path3", 1000)
colnames(PTWAY3a_)<- c("NPP","NEP","NBP","SoilRespiration",
                      "BiomassLiveCStock","SoilCStock","EcosystemCStock",
                      "MineralSoil","Organic","WoodyDebris",
                      "sn","sb","am","af","avf","asl","bgvf","bgf","bgs","Period","Pathway")
PTWAY4a_ <- as.data.frame(PTWAY4_)
PTWAY4a_$Year<-rep(c("Ai","Af","Bi","Bf","Ci","Cf","Di","Df"),each=1000)
PTWAY4a_$Pathway<-rep("Path4", 1000)
colnames(PTWAY4a_)<- c("NPP","NEP","NBP","SoilRespiration",
                      "BiomassLiveCStock","SoilCStock","EcosystemCStock",
                      "MineralSoil","Organic","WoodyDebris",
                      "sn","sb","am","af","avf","asl","bgvf","bgf","bgs","Period","Pathway")
PTWAY5a_ <- as.data.frame(PTWAY5_)
PTWAY5a_$Year<-rep(c("Ai","Af","Bi","Bf","Ci","Cf","Di","Df"),each=1000)
PTWAY5a_$Pathway<-rep("Path5", 1000)
colnames(PTWAY5a_)<- c("NPP","NEP","NBP","SoilRespiration",
                      "BiomassLiveCStock","SoilCStock","EcosystemCStock",
                      "MineralSoil","Organic","WoodyDebris",
                      "sn","sb","am","af","avf","asl","bgvf","bgf","bgs","Period","Pathway")


All2<-rbind(PTWAY1a_,PTWAY2a_,PTWAY3a_,PTWAY4a_,PTWAY5a_)
length(All2[,1])
head(All2)
KgMg<-(0.001)
All2[,1:19] <- with(All2, All2[,1:19]*KgMg)




########################
########################
PTWAY1 <- matrix(0,5000,19,byrow=T)
PTWAY2 <- matrix(0,5000,19,byrow=T)
PTWAY3 <- matrix(0,5000,19,byrow=T)
PTWAY4 <- matrix(0,5000,19,byrow=T)
PTWAY5 <- matrix(0,5000,19,byrow=T)

for (i in 1:19){
 
  
  PTWAY1[,i] <- pathfunction(Path1,i)
  PTWAY2[,i] <- pathfunction(Path2,i)
  PTWAY3[,i] <- pathfunction(Path3,i)
  PTWAY4[,i] <- pathfunction(Path4,i)
  PTWAY5[,i] <- pathfunction(Path5,i)
}

PTWAY1a <- as.data.frame(PTWAY1)
PTWAY1a$Year<-rep(c("I","A","B","C","D"),each=1000)
PTWAY1a$Pathway<-rep("Path1", 1000)
colnames(PTWAY1a)<- c("NPP","NEP","NBP","SoilRespiration",
                      "BiomassLiveCStock","SoilCStock","EcosystemCStock",
                      "MineralSoil","Organic","WoodyDebris",
                       "sn","sb","am","af","avf","asl","bgvf","bgf","bgs","Period","Pathway")
PTWAY2a <- as.data.frame(PTWAY2)
PTWAY2a$Year<-rep(c("I","A","B","C","D"),each=1000)
PTWAY2a$Pathway<-rep("Path2", 1000)
colnames(PTWAY2a)<- c("NPP","NEP","NBP","SoilRespiration",
                      "BiomassLiveCStock","SoilCStock","EcosystemCStock",
                      "MineralSoil","Organic","WoodyDebris",
                      "sn","sb","am","af","avf","asl","bgvf","bgf","bgs","Period","Pathway")
PTWAY3a <- as.data.frame(PTWAY3)
PTWAY3a$Year<-rep(c("I","A","B","C","D"),each=1000)
PTWAY3a$Pathway<-rep("Path3", 1000)
colnames(PTWAY3a)<- c("NPP","NEP","NBP","SoilRespiration",
                      "BiomassLiveCStock","SoilCStock","EcosystemCStock",
                      "MineralSoil","Organic","WoodyDebris",
                      "sn","sb","am","af","avf","asl","bgvf","bgf","bgs","Period","Pathway")
PTWAY4a <- as.data.frame(PTWAY4)
PTWAY4a$Year<-rep(c("I","A","B","C","D"),each=1000)
PTWAY4a$Pathway<-rep("Path4", 1000)
colnames(PTWAY4a)<- c("NPP","NEP","NBP","SoilRespiration",
                      "BiomassLiveCStock","SoilCStock","EcosystemCStock",
                      "MineralSoil","Organic","WoodyDebris",
                      "sn","sb","am","af","avf","asl","bgvf","bgf","bgs","Period","Pathway")
PTWAY5a <- as.data.frame(PTWAY5)
PTWAY5a$Year<-rep(c("I","A","B","C","D"),each=1000)
PTWAY5a$Pathway<-rep("Path5", 1000)
colnames(PTWAY5a)<- c("NPP","NEP","NBP","SoilRespiration",
                      "BiomassLiveCStock","SoilCStock","EcosystemCStock",
                      "MineralSoil","Organic","WoodyDebris",
                      "sn","sb","am","af","avf","asl","bgvf","bgf","bgs","Period","Pathway")



All<-rbind(PTWAY1a,PTWAY2a,PTWAY3a,PTWAY4a,PTWAY5a)
length(All[,1])
head(All)
KgMg<-(0.001)
All[,1:19] <- with(All, All[,1:19]*KgMg)
#All$Period<- ifelse(All$Year<31,"I",ifelse(All$Year>30&All$Year<61,"A", ifelse(All$Year>60&All$Year<91,"B","C")))
#All$Period <- All$Year                                                 
All$Period <- factor(All$Period,
                     levels=c("I",'A','B','C','D'), ordered=TRUE)


All$Pathway <- as.factor(All$Pathway)

##mean annual rate of change over each 30 year climatic period as (F-I)/30. How fast C is lost
#ratechange <- function (x) {
#Hist <- (x[30,1:19] - x[1,1:19])/30
#First <- (x[60,1:19]- x[31,1:19])/30
#Second <- (x[90,1:19]- x[61,1:19])/30
#Third <- (x[120,1:19]- x[91,1:19])/30
#return (list(Hist,First,Second,Third))
#}

ratechangeE <- function (x) {
  Hist <- (x[2] - x[1])/30
  First <- (x[3] - x[2])/30
  Second <- (x[4] - x[3])/30
  Third <- (x[5] - x[4])/30
  return (list(Hist,First,Second,Third))
}

All1<- All[which(All$Pathway=="Path1"),]
rc1 <- tapply(All1$EcosystemCStock, All1$Period, mean)
All2<- All[which(All$Pathway=="Path2"),]
rc2 <-tapply(All2$EcosystemCStock, All2$Period, mean)
All3<- All[which(All$Pathway=="Path3"),]
rc3 <-tapply(All3$EcosystemCStock, All3$Period, mean)
All4<- All[which(All$Pathway=="Path4"),]
rc4 <-tapply(All4$EcosystemCStock, All4$Period, mean)
All5<- All[which(All$Pathway=="Path5"),]
rc5 <-tapply(All5$EcosystemCStock, All5$Period, mean)

All1<- All[which(All$Pathway=="Path1"),]
rc1 <- tapply(All1$MineralSoil, All1$Period, mean)
All2<- All[which(All$Pathway=="Path2"),]
rc2 <-tapply(All2$MineralSoil, All2$Period, mean)
All3<- All[which(All$Pathway=="Path3"),]
rc3 <-tapply(All3$MineralSoil, All3$Period, mean)
All4<- All[which(All$Pathway=="Path4"),]
rc4 <-tapply(All4$MineralSoil, All4$Period, mean)
All5<- All[which(All$Pathway=="Path5"),]
rc5 <-tapply(All5$MineralSoil, All5$Period, mean)

ratechangeeco1 <- ratechangeE(rc1)
ratechangeeco2 <- ratechangeE(rc2)
ratechangeeco3 <- ratechangeE(rc3)
ratechangeeco4 <- ratechangeE(rc4)
ratechangeeco5 <- ratechangeE(rc5)


##Mean percent change in C stocks. How much % C MgC/ha was lost
##Map these results

#meanpercentchangeE <- function (x) {
#Hist <- 100*(x[30,7] - x[1,7])/x[1,7]
#First <- 100*(x[60,7]- x[31,7])/x[31,7]
#Second <- 100*(x[90,7]- x[61,7])/x[61,7]
#Third <- 100*(x[120,7]- x[91,7])/x[91,7]
#return (list(Hist,First,Second,Third))
#}

meanpercentchangeE <- function (x) {
  Hist <- 100*(x[2] - x[1])/x[1]
  First <- 100*(x[3]- x[2])/x[2]
  Second <- 100*(x[4]- x[3])/x[3]
  Third <- 100*(x[5]- x[4])/x[4]
  return (list(Hist,First,Second,Third))
}




#### ANOVA to test for relative effect of climatic period and pathway 
## on each each variable 

aov.fct <- function(x){
  aov1 <- summary(aov(x~Period+Pathway+Period*Pathway,data=All))
  period<- as.character(aov1[[1]][1,4:5])
  path <- as.character(aov1[[1]][2,4:5] )
  inter <- as.character(aov1[[1]][3,4:5])
  return(list(period=period,path=path,inter=inter))
}


summary(aov(SoilRespiration~Period+Pathway+Period*Pathway,data=All))
resanova <- lapply(as.list((All[,1:19])),aov.fct)



##Function
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval:
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

##Fluxes
npp <- summarySE(All2, measurevar="NPP", groupvars=c("Pathway","Period "))
nep <- summarySE(All2, measurevar="NEP", groupvars=c("Pathway","Period"))
nbp <- summarySE(All2, measurevar="NBP", groupvars=c("Pathway","Period"))
soilr <- summarySE(All2, measurevar="SoilRespiration", groupvars=c("Pathway","Period"))
##Stocks
ecs <- summarySE(All2, measurevar="EcosystemCStock", groupvars=c("Pathway","Period"))
soil<- summarySE(All2, measurevar="SoilCStock", groupvars=c("Pathway","Period"))
bcs <- summarySE(All2, measurevar="BiomassLiveCStock", groupvars=c("Pathway","Period"))
ol <- summarySE(All2, measurevar="Organic", groupvars=c("Pathway","Period"))
mscs <- summarySE(All2, measurevar="MineralSoil", groupvars=c("Pathway","Period"))
wd <- summarySE(All2, measurevar="WoodyDebris", groupvars=c("Pathway","Period"))
snags <- summarySE(All, measurevar="sn", groupvars=c("Pathway","Period"))
snabranch <- summarySE(All, measurevar="sb", groupvars=c("Pathway","Period"))
abovemed <- summarySE(All, measurevar="am", groupvars=c("Pathway","Period"))
abovefa <- summarySE(All, measurevar="af", groupvars=c("Pathway","Period"))
abovefa <- summarySE(All, measurevar="af", groupvars=c("Pathway","Period"))
abovevf <- summarySE(All, measurevar="avf", groupvars=c("Pathway","Period"))
abovesl<- summarySE(All, measurevar="asl", groupvars=c("Pathway","Period"))
belowveryfa <- summarySE(All, measurevar="bgvf", groupvars=c("Pathway","Period"))
belowfas <- summarySE(All, measurevar="bgf", groupvars=c("Pathway","Period"))
beloslo <- summarySE(All, measurevar="bgs", groupvars=c("Pathway","Period"))

library(grid)

my_grobA = grobTree(textGrob("A)", x=0.01,  y=0.95, hjust=0,
                             gp=gpar(col="black", fontsize=10, fontface="bold")))
my_grobB = grobTree(textGrob("B)", x=0.01,  y=0.95, hjust=0,
                             gp=gpar(col="black", fontsize=10, fontface="bold")))
my_grobC = grobTree(textGrob("C)", x=0.01,  y=0.95, hjust=0,
                             gp=gpar(col="black", fontsize=10, fontface="bold")))
my_grobD = grobTree(textGrob("D)", x=0.01,  y=0.95, hjust=0,
                             gp=gpar(col="black", fontsize=10, fontface="bold")))
my_grobE = grobTree(textGrob("E)", x=0.01,  y=0.95, hjust=0,
                             gp=gpar(col="black", fontsize=10, fontface="bold")))
my_grobF = grobTree(textGrob("F)", x=0.01,  y=0.95, hjust=0,
                             gp=gpar(col="black", fontsize=10, fontface="bold")))

black.12.text <- element_text(face="plain",color = "black", size =12)
black.bold.text <- element_text(face = "plain", color = "black", size=14)



# # # # # # # #
#Graphics
# # # # # # # #

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#####Fluxes####
p1 <- ggplot(npp, aes(x=Period, y=NPP, group=Pathway)) + ylim(2,4.5) +
  labs(x="Climatic Period",y="NPP (MgC/ha*yr)")+
  geom_errorbar(aes(ymin=NPP-ci, ymax=NPP+ci), colour="black", width=.1)
p1a <- p1+geom_point(aes(shape=factor(Pathway)), size=2) 
p2 <- p1a + scale_shape_manual(values = c(19,6,8,0,4))
p2 <- p2+ theme_bw()
p2a <- p2 + theme(legend.position="none")+ annotation_custom(my_grobA)

p3 <- ggplot(soilr, aes(x=Period, y=SoilRespiration, group=Pathway)) + ylim(2,4.5) +
  labs(x="Climatic Period",y="Rh (MgC/ha*yr)")+
  geom_errorbar(aes(ymin=SoilRespiration-ci, ymax=SoilRespiration+ci), colour="black", width=.1)
p3a <- p3+geom_point(aes(shape=factor(Pathway)), size=2) 
p4 <- p3a + scale_shape_manual(values = c(19,6,8,0,4))
p5 <- p4+ theme_bw()
p5a <- p5 + theme(legend.position="none")+ annotation_custom(my_grobB)

p6 <- ggplot(nep, aes(x=Period, y=NEP, group=Pathway)) + ylim(-0.9,1.5) +
  labs(x="Climatic Period",y="NEP (MgC/ha*yr)")+
  geom_errorbar(aes(ymin=NEP-ci, ymax=NEP+ci), colour="black", width=.1)
p6a <- p6+geom_point(aes(shape=factor(Pathway)), size=2) + geom_hline(yintercept = 0)
p6ab <- p6a + scale_shape_manual(values = c(19,6,8,0,4))
p7 <- p6ab+ theme_bw()
p7a <- p7 + theme(legend.position="none")+ annotation_custom(my_grobC)


p8 <- ggplot(nbp, aes(x=Period, y=NBP, group=Pathway)) + ylim(-1.0,1.5) +
  labs(x="Climatic Period",y="NBP (MgC/ha*yr)")+
  geom_errorbar(aes(ymin=NBP-ci, ymax=NBP+ci), colour="black", width=.1)
p8a <- p8+geom_point(aes(shape=factor(Pathway)), size=2) + geom_hline(yintercept = 0)
p8ab <- p8a + scale_shape_manual(values = c(19,6,8,0,4))
p8abc <- p8ab+ theme_bw()
p9 <- p8abc + theme(legend.position="none")+ annotation_custom(my_grobD)


fluxes <- multiplot(p2a,p7a,p5a,p9 ,cols=2)
tiff('fluxes.tiff', units="in", width=7, height=7, res=300)
fluxes <- multiplot(p2a,p7a,p5a,p9 ,cols=2)
dev.off()

#####Stocks####

p9a <- ggplot(ecs, aes(x=Period, y=EcosystemCStock, group=Pathway)) + ylim(190,262) +
  labs(x="Climatic Period",y="Ecosystem C (MgC/ha)")+
  geom_errorbar(aes(ymin=EcosystemCStock-ci, ymax=EcosystemCStock+ci), colour="black", width=.1)
p9ab <- p9a+geom_point(aes(shape=factor(Pathway)), size=2)
p9abc <- p9ab+theme_bw()
p10 <- p9abc + theme(legend.position="none", title=black.12.text, axis.title = black.12.text, text=black.12.text)
p10a <- p10+annotation_custom(my_grobA)


p10ab <- ggplot(soil, aes(x=Period, y=SoilCStock , group=Pathway)) + ylim(140,205) +
labs(x="Climatic Period",y="DOM C (MgC/ha)")+
geom_errorbar(aes(ymin=SoilCStock-ci, ymax=SoilCStock +ci), colour="black", width=.1)
p10abc <- p10ab + geom_point(aes(shape=factor(Pathway)), size=2)
p11 <- p10abc+theme_bw()
p11a <- p11 + theme(legend.position="none", title=black.12.text, axis.title = black.12.text, text=black.12.text)
p11ab <- p11a + annotation_custom(my_grobB)



p11abc <- ggplot(bcs, aes(x=Period, y=BiomassLiveCStock , group=Pathway)) + ylim(30,65) +
labs(x="Climatic Period",y="Total live Biomass C (MgC/ha)")+
geom_errorbar(aes(ymin=BiomassLiveCStock-ci, ymax=BiomassLiveCStock+ci), colour="black", width=.1)
p12 <- p11abc + geom_point(aes(shape=factor(Pathway)), size=2)
p12a <- p12+theme_bw()
p12ab <- p12a + theme(legend.position="none", title=black.12.text, axis.title = black.12.text, text=black.12.text)
p12abc <- p12ab + annotation_custom(my_grobC)



p13 <- ggplot(ol, aes(x=Period, y=Organic, group=Pathway)) + ylim(30,52) +
labs(x="Climatic Period",y="Organic C (MgC/ha)")+
geom_errorbar(aes(ymin=Organic-ci, ymax=Organic+ci), colour="black", width=.1)
p13a <- p13 + geom_point(aes(shape=factor(Pathway)), size=2)
p13ab <- p13a + theme_bw()
p13abc <- p13ab + theme(legend.position="none", title=black.12.text, axis.title = black.12.text, text=black.12.text)
p14 <- p13abc + annotation_custom(my_grobD)


p14a <- ggplot(wd, aes(x=Period, y=WoodyDebris, group=Pathway)) + ylim(25,50) +
labs(x="Climatic Period",y="Woody Debris C (MgC/ha)")+
geom_errorbar(aes(ymin=WoodyDebris-ci, ymax=WoodyDebris+ci), colour="black", width=.1)
p14ab <- p14a+geom_point(aes(shape=factor(Pathway)), size=2)
p14abc <- p14ab+theme_bw()
p15 <- p14abc + theme(legend.position="none", title=black.12.text, axis.title = black.12.text, text=black.12.text)
p15a <- p15+annotation_custom(my_grobE)


p15ab <- ggplot(mscs, aes(x=Period, y=MineralSoil, group=Pathway)) + ylim(90,105) +
labs(x="Climatic Period",y="Mineral C (MgC/ha)")+
geom_errorbar(aes(ymin=MineralSoil-ci, ymax=MineralSoil+ci), colour="black", width=.1)
p15abc <- p15ab+geom_point(aes(shape=factor(Pathway)), size=2.2)
p16 <- p15abc+theme_bw()
p16a <- p16 + theme(legend.position="none", title=black.12.text, axis.title = black.12.text, text=black.12.text)
p16ab <- p16a+annotation_custom(my_grobF)


stocks <- multiplot(p10a,p14,p11ab,p15a,p12abc,p16ab, cols=3)

tiff('stocks.tiff', units="in", width=7, height=7, res=300)
stocks <- multiplot(p10a,p14,p11ab,p15a,p12abc,p16ab, cols=3)
dev.off()



####mean % change in ecosystem C

meanchangeeco1 <- meanpercentchangeE(rc1)
meanchangeeco2 <- meanpercentchangeE(rc2)
meanchangeeco3 <- meanpercentchangeE(rc3)
meanchangeeco4 <- meanpercentchangeE(rc4)
meanchangeeco5 <- meanpercentchangeE(rc5)

P1 <- c(meanchangeeco1[[1]],meanchangeeco1[[2]],meanchangeeco1[[3]],meanchangeeco1[[4]])
P2 <- c(meanchangeeco2[[1]],meanchangeeco2[[2]],meanchangeeco2[[3]],meanchangeeco2[[4]])
P3 <- c(meanchangeeco3[[1]],meanchangeeco3[[2]],meanchangeeco3[[3]],meanchangeeco3[[4]])
P4 <- c(meanchangeeco4[[1]],meanchangeeco4[[2]],meanchangeeco4[[3]],meanchangeeco4[[4]])
P5 <- c(meanchangeeco5[[1]],meanchangeeco5[[2]],meanchangeeco5[[3]],meanchangeeco5[[4]])
A <- rbind(P1[1],P2[1],P3[1],P4[1],P5[1])
B<- rbind(P1[2],P2[2],P3[2],P4[2],P5[2])
C <- rbind(P1[3],P2[3],P3[3],P4[3],P5[3])
D <- rbind(P1[4],P2[4],P3[4],P4[4],P5[4])
periodos <- rep(c("A","B","C","D"), each=5)
periodos<- factor(periodos,
                     levels=c('A','B','C','D'), ordered=TRUE)
pathways <- rep(c("Pathway1","Pathway2","Pathway3","Pathway4","Pathway5"))
PercentC <- cbind(rbind(A,B,C,D),pathways)
PercentC <- as.data.frame(PercentC)
PercentC$Climatic <- periodos
colnames(PercentC) <- c("Perchange","Pathway","Climatic")
PercentC$Perchange <- as.numeric(as.character(PercentC$Perchange))
#write.csv(PercentC,file="PercentCchange.csv")


mr <- ggplot(PercentC, aes(x=Climatic, y=Perchange, group=Pathway))+ ylim(-11,8) +
labs(x="Climatic Period",y="Mean % change in ecosystem C (MgC/ha)")+ geom_hline(yintercept = 0)
mra <- mr+geom_point(aes(shape=factor(Pathway)), size=3)
mrab <- mra+theme_bw()
mrabc <- mrab + theme(legend.position="none", title=black.12.text, axis.title = black.12.text, text=black.12.text)

tiff('mpc.tiff', units="in", width=7, height=7, res=300)
mrabc
dev.off()

##rate of change

ratechangeeco1 <- ratechangeE(All1)
ratechangeeco2 <- ratechangeE(All2)
ratechangeeco3 <- ratechangeE(All3)
ratechangeeco4 <- ratechangeE(All4)
ratechangeeco5 <- ratechangeE(All5)


P1r <- c(ratechangeeco1[[1]],ratechangeeco1[[2]],ratechangeeco1[[3]],ratechangeeco1[[4]])
P2r <- c(ratechangeeco2[[1]],ratechangeeco2[[2]],ratechangeeco2[[3]],ratechangeeco2[[4]])
P3r <- c(ratechangeeco3[[1]],ratechangeeco3[[2]],ratechangeeco3[[3]],ratechangeeco3[[4]])
P4r <- c(ratechangeeco4[[1]],ratechangeeco4[[2]],ratechangeeco4[[3]],ratechangeeco4[[4]])
P5r <- c(ratechangeeco5[[1]],ratechangeeco5[[2]],ratechangeeco5[[3]],ratechangeeco5[[4]])
Initialr <- rbind(P1r[1],P2r[1],P3r[1],P4r[1],P5r[1])
Ar <- rbind(P1r[2],P2r[2],P3r[2],P4r[2],P5r[2])
Br <- rbind(P1r[3],P2r[3],P3r[3],P4r[3],P5r[3])
Cr <- rbind(P1r[4],P2r[4],P3r[4],P4r[4],P5r[4])
periodos <- rep(c("A","B","C","D"), each=5)
periodos<- factor(periodos,
                  levels=c('A','B','C','D'), ordered=TRUE)
pathways <- rep(c("Pathway1","Pathway2","Pathway3","Pathway4","Pathway5"))
RateC <- cbind(rbind(Initialr,Ar,Br,Cr),pathways)
RateC <- as.data.frame(RateC)
RateC$Climatic <- periodos
colnames(RateC) <- c("Rate","Pathway","Climatic")
RateC$Rate <- as.numeric(as.character(RateC$Rate))

mr2 <- ggplot(RateC, aes(x=Climatic, y=Rate, group=Pathway))+ ylim(-0.8,0.5) +
  labs(x="Climatic Period",y="Rate of change in ecosystem C (MgC/ha*yr)")+ geom_hline(yintercept = 0)
mra2 <- mr2+geom_point(aes(shape=factor(Pathway)), size=3, show.legend = TRUE)
mrab2 <- mra2+theme_bw()

#mrabc2 <- mrab2 + theme(legend.position="none", title=black.12.text, axis.title = black.12.text, text=black.12.text)

tiff('mrcec.tiff', units="in", width=7, height=7, res=300)
mrabc2
dev.off()



###Mean annual temperatures and decay rates

d1 <- as.data.frame(get(load("DecayRatePath1.RData")))
d1$Years <- rep(1:120)
d1$Pathway <-rep("Path1",120)
d2 <- as.data.frame(get(load("DecayRatePath2.RData")))
d2$Years <- rep(1:120)
d2$Pathway <-rep("Path2",120)
d3 <- as.data.frame(get(load("DecayRatePath3.RData")))
d3$Years <- rep(1:120)
d3$Pathway <-rep("Path3",120)
d4 <- as.data.frame(get(load("DecayRatePath4.RData")))
d4$Years <- rep(1:120)
d4$Pathway <-rep("Path4",120)
d5 <- as.data.frame(get(load("DecayRatePath5.RData")))
d5$Years <- rep(1:120)
d5$Pathway <-rep("Path5",120)

decayall <- as.data.frame(rbind(d1,d2,d3,d4,d5))

decaymeans  <- matrix(0,120,8,byrow=T)

for (i in 1:8){
  decaymeans[,i] <-tapply(decayall [,i],decayall $Years,mean)
}

plot(decaymeans[,1] )

###
##BioSimdata

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

Sampled <- subpop (path=c("Pathway17","Pathway19","Pathway20","Pathway21","Pathway1",
                          "Pathway22","Pathway23","Pathway24",
                          "Pathway11","Pathway14","Pathway15","Pathway16",
                          "Pathway26","Pathway28","Pathway1"), d=PoolPlots2) ##path 5



Samples <- matrix(0, 120, 10, byrow = TRUE)
listsamples <- rep(list(Samples),2495)

for (i in 1:2495){
  Samp <- Sampled[i,]
  WeatherPlot <- Samp[1,1]
  ClimateData <- GetClimate (Weather,WeatherPlot)
  ClimateDataclean <- CorrectDupli(ClimateData) 
  listsamples[[i]] <- ClimateDataclean
}

Tprojected <- matrix(0, 2495, 120, byrow = TRUE)
for (i in 1:2495){
  temperaturas <- as.vector(listsamples[[i]][[3]])
  Tprojected[i,] <- temperaturas
}
time <- seq(1:120)
meanstemp <- colMeans(Tprojected)
         
   mean(meanstemp[1:30])    
   mean(meanstemp[31:60])  
   mean(meanstemp[61:90])  
   mean(meanstemp[91:120])  
lines(Tprojected[2,],time)
lines(Tprojected[3,],time)
lines(Tprojected[4,],time)

###
##GrowthIndex

gi1 <- as.data.frame(get(load("GrowthTemp1.RData")))
gi2 <- as.data.frame(get(load("GrowthTemp2.RData")))
gi3 <- as.data.frame(get(load("GrowthTemp3.RData")))
gi4 <- as.data.frame(get(load("GrowthTemp4.RData")))
gi5 <- as.data.frame(get(load("GrowthTemp5.RData")))

years <- seq(1981,2100, by=1) 
GrowthI<- as.data.frame(cbind(gi1[,1],gi2[,1],gi3[,1],gi4[,1],gi5[,1]))

GrowthI$Time <- years
plot(GrowthI$Time,GrowthI$V1, xlab="Time period (years)", ylab="Growth index (m)",pch=1, col="black",ylim=c(5,16))
points(GrowthI$Time,GrowthI$V2, pch=1, col="black")
points(GrowthI$Time,GrowthI$V3,pch=1, col="black")
points(GrowthI$Time,GrowthI$V4,pch=1, col="black")
points(GrowthI$Time,GrowthI$V5,pch=1, col="black")
abline(v=c(2010,2040,2070), col=c("black","black","black"), lty=c(1,1,1), lwd=c(1, 1,1))
text(1995,15.5,"A",cex=1.5, font=2)
text(2025,15.5,"B",cex=1.5, font=2)
text(2055,15.5,"C",cex=1.5, font=2)
text(2085,15.5,"D",cex=1.5, font=2)

growthihist<-colMeans(GrowthI[1:30,])
mean(growthihist)
growthifirst<-colMeans(GrowthI[31:60,])
mean(growthifirst)
growthisecond<-colMeans(GrowthI[61:90,])
mean(growthisecond)
growthithird<-colMeans(GrowthI[91:120,])
mean(growthithird)



PoolPlots <- read.csv("InventoryPlotsPathway.csv",colClasses = "character", header=TRUE,sep = ",")
str(PoolPlots)
PoolPlots[,24:34] <- sapply(PoolPlots[,24:34],as.numeric)
length(PoolPlots[,1])  ##4882...3714
PoolPlots2<-PoolPlots[which(PoolPlots$stands== "EnEn"),]
length(1:dim(PoolPlots2)[1])##3714


str(PoolPlots2)
plots <- PoolPlots2[,c(3,4,30,32,34,36,38)] 
plots$FRI_Pathway <- ifelse(plots$Pathway=="Pathway17"| plots$Pathway=="Pathway19"| 
                              plots$Pathway=="Pathway20"| plots$Pathway=="Pathway21","Pathway1",
                            ifelse(plots$Pathway=="Pathway22"| plots$Pathway=="Pathway23"| 
                                     plots$Pathway=="Pathway24","Pathway2",
                                   ifelse(plots$Pathway=="Pathway11"| plots$Pathway=="Pathway14"| plots$Pathway=="Pathway15"
                                          | plots$Pathway=="Pathway16","Pathway3",ifelse(plots$Pathway=="Pathway26"
                                                                                         | plots$Pathway=="Pathway28","Pathway4",ifelse(plots$Pathway=="Pathway1","Pathway5","Other")))))                                      

plots2 <- plots[which(plots$FRI_Pathway=="Pathway1"|plots$FRI_Pathway=="Pathway2"|plots$FRI_Pathway=="Pathway3"|plots$FRI_Pathway=="Pathway4"|plots$FRI_Pathway=="Pathway5"),]
length(plots2[,1])###2495 plots

p_c <- read.csv("PercentChange.csv", header=TRUE,sep = ",")
colnames(p_c) <- c("Index","Percchange","FRI_Pathway","Climatic_Period")
str(p_c)


merging <- merge(p_c, plots2, by = "FRI_Pathway")
write.csv(merging,file="PercentCchange_spatial.csv")

###supplementary material
my_grobA = grobTree(textGrob("A)", x=0.01,  y=0.95, hjust=0,
                             gp=gpar(col="black", fontsize=12, fontface="bold")))
my_grobB = grobTree(textGrob("B)", x=0.01,  y=0.95, hjust=0,
                             gp=gpar(col="black", fontsize=12, fontface="bold")))
my_grobC = grobTree(textGrob("C)", x=0.01,  y=0.95, hjust=0,
                             gp=gpar(col="black", fontsize=12, fontface="bold")))
my_grobD = grobTree(textGrob("D)", x=0.01,  y=0.95, hjust=0,
                             gp=gpar(col="black", fontsize=12, fontface="bold")))
my_grobE = grobTree(textGrob("E)", x=0.01,  y=0.95, hjust=0,
                             gp=gpar(col="black", fontsize=12, fontface="bold")))
my_grobF = grobTree(textGrob("F)", x=0.01,  y=0.95, hjust=0,
                             gp=gpar(col="black", fontsize=12, fontface="bold")))
my_grobG = grobTree(textGrob("G)", x=0.01,  y=0.95, hjust=0,
                             gp=gpar(col="black", fontsize=12, fontface="bold")))
my_grobH = grobTree(textGrob("H)", x=0.01,  y=0.95, hjust=0,
                             gp=gpar(col="black", fontsize=12, fontface="bold")))
my_grobI = grobTree(textGrob("I)", x=0.01,  y=0.95, hjust=0,
                             gp=gpar(col="black", fontsize=12, fontface="bold")))

snags <- summarySE(All, measurevar="sn", groupvars=c("Pathway","Period"))
snabranch <- summarySE(All, measurevar="sb", groupvars=c("Pathway","Period"))
abovemed <- summarySE(All, measurevar="am", groupvars=c("Pathway","Period"))
abovefa <- summarySE(All, measurevar="af", groupvars=c("Pathway","Period"))
abovefa <- summarySE(All, measurevar="af", groupvars=c("Pathway","Period"))
abovevf <- summarySE(All, measurevar="avf", groupvars=c("Pathway","Period"))
abovesl<- summarySE(All, measurevar="asl", groupvars=c("Pathway","Period"))
belowveryfa <- summarySE(All, measurevar="bgvf", groupvars=c("Pathway","Period"))
belowfas <- summarySE(All, measurevar="bgf", groupvars=c("Pathway","Period"))
beloslo <- summarySE(All, measurevar="bgs", groupvars=c("Pathway","Period"))


p9a <- ggplot(snags, aes(x=Period, y=sn, group=Pathway)) + ylim(3,7) +
  labs(x="Climatic Period",y="Snags (MgC/ha)")+
  geom_errorbar(aes(ymin=sn-ci, ymax=sn+ci), colour="black", width=.1)
p9ab <- p9a+geom_point(aes(shape=factor(Pathway)), size=2)
p9abc <- p9ab+theme_bw()
p10 <- p9abc + theme(legend.position="none", title=black.12.text, axis.title = black.12.text, text=black.12.text)
p10a <- p10+annotation_custom(my_grobA)


p10ab <- ggplot(snabranch, aes(x=Period, y=sb , group=Pathway)) + ylim(1,2) +
  labs(x="Climatic Period",y="SnagBranch (MgC/ha)")+
  geom_errorbar(aes(ymin=sb-ci, ymax=sb +ci), colour="black", width=.1)
p10abc <- p10ab + geom_point(aes(shape=factor(Pathway)), size=2)
p11 <- p10abc+theme_bw()
p11a <- p11 + theme(legend.position="none", title=black.12.text, axis.title = black.12.text, text=black.12.text)
p11ab <- p11a + annotation_custom(my_grobB)



p<-ggplot(abovemed, aes(x=Period, y=am, group=Pathway)) + ylim(15,30) +
  labs(x="Climatic Period",y="AG Medium (MgC/ha)")+
  geom_errorbar(aes(ymin=am-ci, ymax=am+ci), colour="black", width=.1)
p1<-p+geom_point(aes(shape=factor(Pathway)), size=2)
p2<-p1+theme_bw()
p7a<-p2+theme( legend.position="none", title=black.12.text, axis.title = black.12.text, text=black.12.text)
p9<-p7a+annotation_custom(my_grobC)


p<-ggplot(abovefa, aes(x=Period, y=af, group=Pathway)) + ylim(3,12) +
  labs(x="Climatic Period",y="AG fast (MgC/ha)")+
  geom_errorbar(aes(ymin=af-ci, ymax=af+ci), colour="black", width=.1)
p1<-p+geom_point(aes(shape=factor(Pathway)), size=2)
p2<-p1+theme_bw()
p7a<-p2+theme( legend.position="none", title=black.12.text, axis.title = black.12.text, text=black.12.text)
p10<-p7a+annotation_custom(my_grobD)



p<-ggplot(abovevf, aes(x=Period, y=avf, group=Pathway)) + ylim(2,14) +
  labs(x="Climatic Period",y="AG very fast (MgC/ha)")+
  geom_errorbar(aes(ymin=avf-ci, ymax=avf+ci), colour="black", width=.1)
p1<-p+geom_point(aes(shape=factor(Pathway)), size=2)
p2<-p1+theme_bw()
p7a<-p2+theme( legend.position="none", title=black.12.text, axis.title = black.12.text, text=black.12.text)
p11<-p7a+annotation_custom(my_grobE)



p<-ggplot(abovesl, aes(x=Period, y=asl, group=Pathway)) + ylim(29,40) +
  labs(x="Climatic Period",y="AG slow (MgC/ha)")+
  geom_errorbar(aes(ymin=asl-ci, ymax=asl+ci), colour="black", width=.1)
p1<-p+geom_point(aes(shape=factor(Pathway)), size=2)
p2<-p1+theme_bw()
p7a<-p2+theme( legend.position="none", title=black.12.text, axis.title = black.12.text, text=black.12.text)
p12<-p7a+annotation_custom(my_grobF)



p<-ggplot(belowveryfa, aes(x=Period, y=bgvf,group=Pathway)) + ylim(0,2.5) +
  labs(x="Climatic Period",y="BG very fast (MgC/ha)")+
  geom_errorbar(aes(ymin=bgvf-ci, ymax=bgvf+ci), colour="black", width=.1)
p1<-p+geom_point(aes(shape=factor(Pathway)), size=2)
p2<-p1+theme_bw()
p7a<-p2+theme( legend.position="none", title=black.12.text, axis.title = black.12.text, text=black.12.text)
p13<-p7a+annotation_custom(my_grobG)



p<-ggplot(belowfas, aes(x=Period, y=bgf, group=Pathway)) + ylim(1,5.5) +
  labs(x="Climatic Period",y="BG fast (MgC/ha)")+
  geom_errorbar(aes(ymin=bgf-ci, ymax=bgf+ci), colour="black", width=.1)
p1<-p+geom_point(aes(shape=factor(Pathway)), size=2)
p2<-p1+theme_bw()
p7a<-p2+theme( legend.position="none", title=black.12.text, axis.title = black.12.text, text=black.12.text)
p14<-p7a+annotation_custom(my_grobH)




p<-ggplot(beloslo, aes(x=Period, y=bgs,group=Pathway)) + ylim(90,105) +
  labs(x="Climatic Period",y="BG slow (MgC/ha)")+
  geom_errorbar(aes(ymin=bgs-ci, ymax=bgs+ci), colour="black", width=.1)
p1<-p+geom_point(aes(shape=factor(Pathway)), size=2)
p2<-p1+theme_bw()
p7a<-p2+theme( legend.position="none", title=black.12.text, axis.title = black.12.text, text=black.12.text)
p15<-p7a+annotation_custom(my_grobI)


nnn<-multiplot(p10a,p10,p13,p11ab,p11,p14,p9,p12,p15,cols=3)

tiff('figureS5.tiff', units="in", width=9, height=9, res=300)
nnn<-multiplot(p10a,p10,p13,p11ab,p11,p14,p9,p12,p15,cols=3)
dev.off()



###Only Fire
p1 <- get(load("Path1onlyFRI.RData"))
p2 <- get(load("Path2onlyFRI.RData"))
p3 <- get(load("Path3onlyFRI.RData"))
p4 <- get(load("Path4onlyFRI.RData"))
p5 <- get(load("Path5onlyFRI.RData"))

p1 <- Path1onlyFRI
p2 <- Path2onlyFRI
p3 <- Path3onlyFRI
p4 <- Path4onlyFRI
p5 <- Path5onlyFRI

PTWAY1F <- matrix(0,5000,19,byrow=T) ##only fire
PTWAY2F <- matrix(0,5000,19,byrow=T)
PTWAY3F <- matrix(0,5000,19,byrow=T)
PTWAY4F <- matrix(0,5000,19,byrow=T)
PTWAY5F <- matrix(0,5000,19,byrow=T)

for (i in 1:19){
  
  
  PTWAY1F[,i] <- pathfunction(p1,i)
  PTWAY2F[,i] <- pathfunction(p2,i)
  PTWAY3F[,i] <- pathfunction(p3,i)
  PTWAY4F[,i] <- pathfunction(p4,i)
  PTWAY5F[,i] <- pathfunction(p5,i)
}

PTWAY1Fa <- as.data.frame(PTWAY1F)
PTWAY1Fa$Year<-rep(c("I","A","B","C","D"),each=1000)
PTWAY1Fa$Pathway<-rep("Path1", 1000)
colnames(PTWAY1Fa)<- c("NPP","NEP","NBP","SoilRespiration",
                      "BiomassLiveCStock","SoilCStock","EcosystemCStock",
                      "MineralSoil","Organic","WoodyDebris",
                      "sn","sb","am","af","avf","asl","bgvf","bgf","bgs","Period","Pathway")
PTWAY2Fa <- as.data.frame(PTWAY2F)
PTWAY2Fa$Year<-rep(c("I","A","B","C","D"),each=1000)
PTWAY2Fa$Pathway<-rep("Path2", 1000)
colnames(PTWAY2Fa)<- c("NPP","NEP","NBP","SoilRespiration",
                      "BiomassLiveCStock","SoilCStock","EcosystemCStock",
                      "MineralSoil","Organic","WoodyDebris",
                      "sn","sb","am","af","avf","asl","bgvf","bgf","bgs","Period","Pathway")
PTWAY3Fa <- as.data.frame(PTWAY3F)
PTWAY3Fa$Year<-rep(c("I","A","B","C","D"),each=1000)
PTWAY3Fa$Pathway<-rep("Path3", 1000)
colnames(PTWAY3Fa)<- c("NPP","NEP","NBP","SoilRespiration",
                      "BiomassLiveCStock","SoilCStock","EcosystemCStock",
                      "MineralSoil","Organic","WoodyDebris",
                      "sn","sb","am","af","avf","asl","bgvf","bgf","bgs","Period","Pathway")
PTWAY4Fa <- as.data.frame(PTWAY4F)
PTWAY4Fa$Year<-rep(c("I","A","B","C","D"),each=1000)
PTWAY4Fa$Pathway<-rep("Path4", 1000)
colnames(PTWAY4Fa)<- c("NPP","NEP","NBP","SoilRespiration",
                      "BiomassLiveCStock","SoilCStock","EcosystemCStock",
                      "MineralSoil","Organic","WoodyDebris",
                      "sn","sb","am","af","avf","asl","bgvf","bgf","bgs","Period","Pathway")
PTWAY5Fa <- as.data.frame(PTWAY5F)
PTWAY5Fa$Year<-rep(c("I","A","B","C","D"),each=1000)
PTWAY5Fa$Pathway<-rep("Path5", 1000)
colnames(PTWAY5Fa)<- c("NPP","NEP","NBP","SoilRespiration",
                      "BiomassLiveCStock","SoilCStock","EcosystemCStock",
                      "MineralSoil","Organic","WoodyDebris",
                      "sn","sb","am","af","avf","asl","bgvf","bgf","bgs","Period","Pathway")


All_f<-rbind(PTWAY1Fa,PTWAY2Fa,PTWAY3Fa,PTWAY4Fa,PTWAY5Fa)
length(All_f[,1])
head(All_f)
KgMg<-(0.001)
All_f[,1:19] <- with(All_f, All_f[,1:19]*KgMg)
#All$Period<- ifelse(All$Year<31,"I",ifelse(All$Year>30&All$Year<61,"A", ifelse(All$Year>60&All$Year<91,"B","C")))
#All$Period <- All$Year                                                 
All_f$Period <- factor(All_f$Period,
                     levels=c("I",'A','B','C','D'), ordered=TRUE)


All_f$Pathway <- as.factor(All_f$Pathway)
##Fluxes
npp_f <- summarySE(All_f, measurevar="NPP", groupvars=c("Pathway","Period "))
nep_f <- summarySE(All_f, measurevar="NEP", groupvars=c("Pathway","Period"))
nbp_f <- summarySE(All_f, measurevar="NBP", groupvars=c("Pathway","Period"))
soilr_f <- summarySE(All_f, measurevar="SoilRespiration", groupvars=c("Pathway","Period"))
##Stocks
ecs_f <- summarySE(All_f, measurevar="EcosystemCStock", groupvars=c("Pathway","Period"))
soil_f<- summarySE(All_f, measurevar="SoilCStock", groupvars=c("Pathway","Period"))
bcs_f <- summarySE(All_f, measurevar="BiomassLiveCStock", groupvars=c("Pathway","Period"))
ol_f <- summarySE(All_f, measurevar="Organic", groupvars=c("Pathway","Period"))
mscs_f <- summarySE(All_f, measurevar="MineralSoil", groupvars=c("Pathway","Period"))
wd_f <- summarySE(All_f, measurevar="WoodyDebris", groupvars=c("Pathway","Period"))


p9a <- ggplot(ecs_f, aes(x=Period, y=EcosystemCStock, group=Pathway)) + ylim(220,262) +
  labs(x="Climatic Period",y="Total ecosystem C (MgC/ha)")+
  geom_errorbar(aes(ymin=EcosystemCStock-ci, ymax=EcosystemCStock+ci), colour="black", width=.1)
p9ab <- p9a+geom_point(aes(shape=factor(Pathway)), size=2)
p9abc <- p9ab+theme_bw()
p10 <- p9abc + theme(legend.position="none", title=black.12.text, axis.title = black.12.text, text=black.12.text)
p10a <- p10+annotation_custom(my_grobA)


p10ab <- ggplot(soil_f, aes(x=Period, y=SoilCStock , group=Pathway)) + ylim(160,205) +
  labs(x="Climatic Period",y="DOM C (MgC/ha)")+
  geom_errorbar(aes(ymin=SoilCStock-ci, ymax=SoilCStock +ci), colour="black", width=.1)
p10abc <- p10ab + geom_point(aes(shape=factor(Pathway)), size=2)
p11 <- p10abc+theme_bw()
p11a <- p11 + theme(legend.position="none", title=black.12.text, axis.title = black.12.text, text=black.12.text)
p11ab <- p11a + annotation_custom(my_grobB)



p11abc <- ggplot(bcs_f, aes(x=Period, y=BiomassLiveCStock , group=Pathway)) + ylim(35,65) +
  labs(x="Climatic Period",y="Total live Biomass C (MgC/ha)")+
  geom_errorbar(aes(ymin=BiomassLiveCStock-ci, ymax=BiomassLiveCStock+ci), colour="black", width=.1)
p12 <- p11abc + geom_point(aes(shape=factor(Pathway)), size=2)
p12a <- p12+theme_bw()
p12ab <- p12a + theme(legend.position="none", title=black.12.text, axis.title = black.12.text, text=black.12.text)
p12abc <- p12ab + annotation_custom(my_grobC)



p13 <- ggplot(ol_f, aes(x=Period, y=Organic, group=Pathway)) + ylim(35,50) +
  labs(x="Climatic Period",y="Organic C (MgC/ha)")+
  geom_errorbar(aes(ymin=Organic-ci, ymax=Organic+ci), colour="black", width=.1)
p13a <- p13 + geom_point(aes(shape=factor(Pathway)), size=2)
p13ab <- p13a + theme_bw()
p13abc <- p13ab + theme(legend.position="none", title=black.12.text, axis.title = black.12.text, text=black.12.text)
p14 <- p13abc + annotation_custom(my_grobD)


p14a <- ggplot(wd_f, aes(x=Period, y=WoodyDebris, group=Pathway)) + ylim(30,50) +
  labs(x="Climatic Period",y="Woody Debris C (MgC/ha)")+
  geom_errorbar(aes(ymin=WoodyDebris-ci, ymax=WoodyDebris+ci), colour="black", width=.1)
p14ab <- p14a+geom_point(aes(shape=factor(Pathway)), size=2)
p14abc <- p14ab+theme_bw()
p15 <- p14abc + theme(legend.position="none", title=black.12.text, axis.title = black.12.text, text=black.12.text)
p15a <- p15+annotation_custom(my_grobE)


p15ab <- ggplot(mscs_f, aes(x=Period, y=MineralSoil, group=Pathway)) + ylim(90,105) +
  labs(x="Climatic Period",y="Mineral C (MgC/ha)")+
  geom_errorbar(aes(ymin=MineralSoil-ci, ymax=MineralSoil+ci), colour="black", width=.1)
p15abc <- p15ab+geom_point(aes(shape=factor(Pathway)), size=2.2)
p16 <- p15abc+theme_bw()
p16a <- p16 + theme(legend.position="none", title=black.12.text, axis.title = black.12.text, text=black.12.text)
p16ab <- p16a+annotation_custom(my_grobF)


stocks <- multiplot(p10a,p14,p11ab,p15a,p12abc,p16ab, cols=3)

tiff('stocks.tiff', units="in", width=7, height=7, res=300)
stocks <- multiplot(p10a,p14,p11ab,p15a,p12abc,p16ab, cols=3)
dev.off()
