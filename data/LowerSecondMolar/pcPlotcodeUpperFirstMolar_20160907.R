library(vegan)
library(geomorph)
library(shapes)
library(devtools)
install_github("vbonhomme/Momocs")
library(Momocs)

### Plotting function to plot convex hulls
### Filename: Plot_ConvexHull.R
### Notes:
############################################################################

# INPUTS:
# xcoords: x-coordinates of point data
# ycoords: y-coordinates of point data
# lcolor: line color

# OUTPUTS:
# convex hull around data points in a particular color (specified by lcolor)

# FUNCTION:
Plot_ConvexHull<-function(xcoord, ycoord, lcolor){
  hpts <- chull(x = xcoord, y = ycoord)
  hpts <- c(hpts, hpts[1])
  polygon(xcoord[hpts], ycoord[hpts], col = rgb(0,0,0,0), border = lcolor, lwd = 2)
  
}  
# END OF FUNCTION

evaluate_harmonics <- function(coe, point_count){#this function takes original EF coefficient data and evaluates the number of harmonics that have to be used to catch 95-99% of the variance/power of the outlines; first harmonic is usually excluded. To include it, remove the [-1]'s.
  stopifnot(is.data.frame(coe))
  coef <- coe[,3:ncol(coe)] # remove non-coefficient columns (may change depending on what other variables are inserted besides areas and ID)
  co <- coef^2
  tharm <- apply(co, 2, sum)
  power <- apply(matrix(tharm, (point_count/2),4),1, sum) #group tharm into 30 harmonics (rows) each 	with 4 coeffs representing sum over all An, Bn, Cn, Dn, and sum by rows
  cumpower <- round(cumsum(power[-1])/sum(power[-1]),3)[1:15]# take out first harmonic and then calculate cum power
  vharm <- apply(coef, 2, var)
  variation <- apply(matrix(vharm, (point_count/2), 4),1,sum)
  cumvar <- round(cumsum(variation[-1])/sum(variation[-1]),3)[1:15]# cumulative variance excluding the first harmonic
  prelim.table <- rbind(cumpower, cumvar)
  rownames(prelim.table) <- c('Cumulative Power','Cumulative Variance')
  colnames(prelim.table) <- c('2nd','3rd','4th','5th','6th','7th','8th','9th','10th','11th','12th','13th','14th','15th','16th')
  return(prelim.table)
}	

evaluate_harmonics(outCoo$coo,100)


#Only needs to be run once:
#Shell script
#ImageMagick
folderGold <- "/Users/gregorymatthews/Dropbox/brophyTeeth/DeciduousTeethProject/"
#Convert to .jpg
#system(paste("cd ",folderGold,'; for i in *.png ; do convert "$i" "${i%.*}.jpg" ; done',sep=" "))
#Convert to gray scale
#system(paste("cd ",folderGold,'; for i in *.jpg; do convert "$i" -colorspace Gray "$i"; done',sep=" "))
#Trim to center
#system(paste("cd ",folderGold,'; for i in *.jpg; do convert -trim "$i" "$i"; done',sep=" "))


#Import MTurk results
#folders<-list.files("/Users/gregorymatthews/Dropbox/brophyTeeth/MechanicalTurk-May-110imagescopy/")

#for (j in folders){print(j)
#Pick out one folder
#folderTemp <- paste("/Users/gregorymatthews/Dropbox/brophyTeeth/MechanicalTurk-May-110imagescopy/",j,sep="")
#fff<-list.files(folderTemp,full=TRUE)

#Shell script
#ImageMagick
#system(paste("cd ",folderTemp,'; for i in *.png ; do convert "$i" -background white -alpha remove "$i" ; done',sep=" "))
#Convert to .jpg
#system(paste("cd ",folderTemp,'; for i in *.png ; do convert "$i" "${i%.*}.jpg" ; done',sep=" "))
#Convert to gray scale
#system(paste("cd ",folderTemp,'; for i in *.jpg; do convert "$i" -colorspace Gray "$i"; done',sep=" "))
#trim
#system(paste("cd ",folderTemp,'; for i in *.jpg; do convert -trim "$i" "$i"; done',sep=" "))
#}


#subsetFiles <- read.csv("/Users/gregorymatthews/Dropbox/brophyTeeth/DeciduousTeethProject/UpperFirstMolar/Udm2 TEST w no paleo.csv")
#subsetFiles$Udm2 <- as.character(subsetFiles$Udm2)
#subsetFiles$X <- as.character(subsetFiles$X)

key <- read.csv("/Users/gregorymatthews/Dropbox/brophyTeeth/DeciduousTeethProject/UpperFirstMolar/Udm1.csv")[-1,1:3]
names(key)<-c("SpecimenNumber","ID","Group")
key$ID<-as.character(key$ID)
key$ID <- gsub(" ","",key$ID)
key$ID <- gsub("_","",key$ID)
key$ID <- toupper(key$ID)
####################################################################
##Juliet's files 
####################################################################
g<-1
#Get the file names of Juliet's b/w images.
goldFiles<-list.files("/Users/gregorymatthews/Dropbox/brophyTeeth/DeciduousTeethProject/UpperFirstMolar/")
goldFiles <- goldFiles[grep("jpg",goldFiles)]
#goldFiles <- goldFiles[goldFiles%in%subsetFiles$X]

#I get an error with this tooth: "DSCN2879"
#I have no idea why. So I removed it for now.
#goldFiles <- goldFiles[-grep("DSCN2879",goldFiles)]
imgsList<-list()
for (g in c(1:length(goldFiles))){print(g)
  #for (g in 1:10){print(g)
  #Read in the files from Juliet
  imgID<-goldFiles[g]
  id<-strsplit(imgID,"[.]")[[1]][1]
  
  
  #Contains the image done by Juliet
  imgsList[[imgID]]<-import_jpg(paste("/Users/gregorymatthews/Dropbox/brophyTeeth/DeciduousTeethProject/UpperFirstMolar/",imgID,sep=""))[[1]]
  #temp<-import_jpg(paste("/Users/gregorymatthews/Dropbox/brophyTeeth/DeciduousTeethProject/UpperFirstMolar/",imgID,sep=""))[[1]]
  #harms<-efourier(temp,nb.h=50)
  #imgsList[[imgID]]<-efourier_shape(harms$an,harms$bn,harms$cn,harms$dn,nb.h=50,nb.pts=100)
}

outCoo<-Out(imgsList)
test<-efourier(outCoo,nb.h=20,norm=FALSE)

#proc<-fgProcrustes(outCoo)
#png("/Users/gregorymatthews/Dropbox/brophyTeeth/DeciduousTeethProject/UpperFirstMolarProcrustesStack.png")
#stack(proc,title="Upper Second Molar - Procrustes Stack")
#dev.off()
#Do EFA

#test<-efourier(Out(proc$coo),nb.h=15,norm=TRUE)

names(imgsList) <- toupper(names(imgsList))
names(imgsList) <- gsub(" ","",names(imgsList))
names(imgsList) <- gsub("_","",names(imgsList))

group<-rep(NA,nrow(test))
for (i in 1:length(group)){
  group[i] <- as.character(key$Group[key$ID==names(imgsList)[i]])
}

harmId<-c(paste0("A",c(1:8)),paste0("B",c(1:8)),paste0("C",c(1:8)),paste0("D",c(1:8)))

x<-test$coe[,harmId]
x<-x[complete.cases(x),]
#xScale<-apply(x,2,scale)
#cor(xScale)

pc<-princomp(x[,!colnames(x)%in%c("A1","B1","C1","D1")])


png("/Users/gregorymatthews/Dropbox/brophyTeeth/DeciduousTeethProject/UDM1_pcPLot_withLabels_TEST_20160907_ALLDATA.png",h=10,w=10,units="in",res=300)
col<-rep("black",length(group))
col[group=="Recent"]<-"gray"
col[group=="Hnaledi"]<-"red"
col[group=="Aafarensis"]<-"blue"
col[group=="Aafricanus"]<-"lightblue"
col[group=="EarlyHomosapiens"]<-"green"
col[group=="Homosp"]<-"lightgreen"
col[group=="Neandertals"]<-"pink"
col[group=="Probustus"]<-"purple"
col[group=="UP"]<-"orange"
col[group=="Omo"]<-"gold"
col[group=="Pboisei"]<-"purple2"

plot(-pc$scores[,2],pc$scores[,1],pch=16,col="white",cex=2,main="Principal Component Plot",sub="UDM1",xlab="PC1",ylab="PC2",xlim=c(-30,65))
legend(35,0,c("Hnaledi","Recent","Aafarensis","Aafricanus","Homosp","EarlyHomoSapiens","Neandertals","Probustus","UP","Omo","Pboisei"),
       c("red","gray","blue","lightblue","lightgreen","green","pink","purple","orange","gold","purple2"))
text(-pc$scores[,2],pc$scores[,1],substring(names(imgsList),1,8),col="black",cex=0.5)
for (i in 1:length(imgsList)){
  points(.30*scale(imgsList[[i]][,2])-pc$scores[i,2],.30*scale(imgsList[[i]][,1])+pc$scores[i,1],type="l",col=col[i])
}

xrange <- range(-pc$scores[,2])
yrange <- range(pc$scores[,1])

spec <-names(table(key$Group)[table(key$Group)>1])
for (sp in spec){
  temp<-pc$scores[group==sp,1:2]
  # Plot it up!
  par(tck = 0.02, mgp = c(1.7, 0.3, 0))
  #plot(pc$scores[,1], pc$scores[,2], type = "p", pch = 1, col = "black", xlim = c(xrange), ylim = c(yrange))
  Plot_ConvexHull(xcoord = -temp[,2], ycoord = temp[,1], lcolor = col[group==sp])
}
#http://127.0.0.1:42290/graphics/plot_zoom_png?width=528&height=667
dev.off()
