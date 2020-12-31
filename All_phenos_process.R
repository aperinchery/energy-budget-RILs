#read in table
setwd("C:/Users/Anna/Dropbox/AMP_EGK/Projects/Ch1_RILs/Scripts")

dat<-read.table("../RawData/Combined/completeplatesNAchanged.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE)
dat$MASS <- NULL
dat$X <- NULL
dat$X.1 <- NULL
dat$X.2 <- NULL
dat$X.3 <- NULL
dat$X.4 <- NULL

c_std<-read.table("../RawData/StdCurves/carbstdcurvev4.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE)
g_std<-read.table("../RawData/StdCurves/glystdcurvev4.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE)
p_std<-read.table("../RawData/StdCurves/protienstdcurvev5.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE)
l525_std<-read.table("../RawData/StdCurves/lipid525stdcurvev6.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE)

#avgstdc<-read.table("../RawData/StdCurves/c_standardcurve",sep="\t",header=TRUE,stringsAsFactors=FALSE)
#avgstdg<-read.table("../RawData/StdCurves/g_standardcurve",sep="\t",header=TRUE,stringsAsFactors=FALSE)
avgstdl525<-read.table("../RawData/StdCurves/l525_standardcurve",sep="\t",header=TRUE,stringsAsFactors=FALSE)
#avgstdp<-read.table("../RawData/StdCurves/p_standardcurve",sep="\t",header=TRUE,stringsAsFactors=FALSE)


#add treatmet and RIL column
dat$RIL<-as.numeric(gsub("[^[:digit:]]", "", dat$Line))#Make RIL column
dat$Treat<-gsub("[[:digit:]]","",dat$Line)#Make Treatment column

which(dat$RIL>2000)
#loop for carbohydrates
carbohydrates<-data.frame("amount"=numeric(length=0), "line"=numeric(length=0), "diet"=character(length=0), "plate"=character(length=0), "flies"=character(length=0), stringsAsFactors=FALSE)

for (plateid in unique(dat$PLATE))
{
  std.c<-subset(c_std, PLATE== plateid) #standard curve, subsetted by plate
  star<-lm(abs~amount, data=c_std) #slope model (amt vs absorbance), according to stds
  plate.s<-subset(dat,PLATE== plateid) #subsetting data by plate
  workwork<-(plate.s$CARB625-star$coefficients[1])/star$coefficients[2] #converts absobance to macromolecule amount
  carbohydrates.s<-data.frame("amount"=workwork, "line"=plate.s$RIL, "diet"=plate.s$Treat, "plate"=plate.s$PLATE, "flies"=plate.s$FLIES, stringsAsFactors=FALSE) # make dataset to fill in info
  carbohydrates<-rbind(carbohydrates, carbohydrates.s) #make new dataset from new info
}
#calculate amount per fly.
dat$c_perfly<-(carbohydrates$amount*(150/1400))/(carbohydrates$flies)

#GLYCOGEN
glycogen<-data.frame("amount"=numeric(length=0), "line"=numeric(length=0), "diet"=character(length=0), "plate"=character(length=0), "flies"=character(length=0), stringsAsFactors=FALSE)

for (plateid in unique(dat$PLATE))
{
  std.g<-subset(g_std, PLATE== plateid)
  star<-lm(abs~amount, data=c_std)
  plate.s<-subset(dat,PLATE== plateid)
  workwork<-(plate.s$GLY625-star$coefficients[1])/star$coefficients[2]
  glycogen.s<-data.frame("amount"=workwork, "line"=plate.s$RIL, "diet"=plate.s$Treat, "plate"=plate.s$PLATE, "flies"=plate.s$FLIES, stringsAsFactors=FALSE)
  glycogen<-rbind(glycogen, glycogen.s)
}

#calculate amount per fly.
dat$g_perfly<-(glycogen$amount*(150/1000))/(glycogen$flies)

###loop for lipids
lipids525<-data.frame("amount"=numeric(length=0), "line"=numeric(length=0), "diet"=character(length=0), "plate"=character(length=0), "flies"=character(length=0), stringsAsFactors=FALSE)

for (plateid in unique(dat$PLATE))
{
  #std.p<-subset(p_std, PLATE== plateid)
  std.l525<-avgstdl525
  star<-lm(abs~amount, data=l525_std)
  plate.s<-subset(dat,PLATE== plateid)
  workwork<-(plate.s$LIP525-star$coefficients[1])/star$coefficients[2]
  lipids525.s<-data.frame("amount"=workwork, "line"=plate.s$RIL, "diet"=plate.s$Treat, "plate"=plate.s$PLATE, "flies"=plate.s$FLIES, stringsAsFactors=FALSE)
  lipids525<-rbind(lipids525, lipids525.s)
}

#calculate amount per fly.
dat$lperfly<-(lipids525$amount*(100/1400))/(lipids525$flies)


###loop for protiens
protiens<-data.frame("amount"=numeric(length=0), "line"=numeric(length=0), "diet"=character(length=0), "plate"=character(length=0), "flies"=character(length=0), stringsAsFactors=FALSE)

for (plateid in unique(dat$PLATE))
{
  std.p<-subset(p_std, PLATE== plateid)
  star<-lm(abs~amount, data=p_std)
  plate.s<-subset(dat,PLATE== plateid)
  workwork<-(plate.s$PROT595-star$coefficients[1])/star$coefficients[2]
  protiens.s<-data.frame("amount"=workwork, "line"=plate.s$RIL, "diet"=plate.s$Treat, "plate"=plate.s$PLATE, "flies"=plate.s$FLIES, stringsAsFactors=FALSE)
  protiens<-rbind(protiens, protiens.s)
}

###calculate protien amount per fly.
dat$p_perfly<-(protiens$amount*(2/180))/(protiens$flies)

which(dat$RIL>2000)

##remove the 22000 from 'RIL' here
dat$RIL[which(dat$RIL>2000)]<-dat$RIL[which(dat$RIL>2000)]-22000
which(dat$RIL>2000)

write.table(dat, "../ProcessedData/Processed_perflyconc.txt", sep="\t")

tests<-count(dat, Line)
biological_replicates<-subset(tests, n>2, select= c("Line", "n"))