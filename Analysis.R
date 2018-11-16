
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

aov.fct(All$NPP)
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
npp <- summarySE(All, measurevar="NPP", groupvars=c("Pathway","Period "))
nep <- summarySE(All, measurevar="NEP", groupvars=c("Pathway","Period"))
nbp <- summarySE(All, measurevar="NBP", groupvars=c("Pathway","Period"))
soilr <- summarySE(All, measurevar="SoilRespiration", groupvars=c("Pathway","Period"))
##Stocks
ecs <- summarySE(All, measurevar="EcosystemCStock", groupvars=c("Pathway","Period"))
soil<- summarySE(All, measurevar="SoilCStock", groupvars=c("Pathway","Period"))
bcs <- summarySE(All, measurevar="BiomassLiveCStock", groupvars=c("Pathway","Period"))
ol <- summarySE(All, measurevar="Organic", groupvars=c("Pathway","Period"))
mscs <- summarySE(All, measurevar="MineralSoil", groupvars=c("Pathway","Period"))
wd <- summarySE(All, measurevar="WoodyDebris", groupvars=c("Pathway","Period"))
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
p8abc <- p8a+ theme_bw()
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
Initial <- rbind(P1[1],P2[1],P3[1],P4[1],P5[1])
A <- rbind(P1[2],P2[2],P3[2],P4[2],P5[2])
B <- rbind(P1[3],P2[3],P3[3],P4[3],P5[3])
C <- rbind(P1[4],P2[4],P3[4],P4[4],P5[4])
periodos <- rep(c("A","B","C","D"), each=5)
periodos<- factor(periodos,
                     levels=c('A','B','C','D'), ordered=TRUE)
pathways <- rep(c("Pathway1","Pathway2","Pathway3","Pathway4","Pathway5"))
PercentC <- cbind(rbind(Initial,A,B,C),pathways)
PercentC <- as.data.frame(PercentC)
PercentC$Climatic <- periodos
colnames(PercentC) <- c("Perchange","Pathway","Climatic")
PercentC$Perchange <- as.numeric(as.character(PercentC$Perchange))

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
mra2 <- mr2+geom_point(aes(shape=factor(Pathway)), size=3)
mrab2 <- mra2+theme_bw()
mrabc2 <- mrab2 + theme(legend.position="none", title=black.12.text, axis.title = black.12.text, text=black.12.text)

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
