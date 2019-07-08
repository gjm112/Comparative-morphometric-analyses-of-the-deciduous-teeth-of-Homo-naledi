library(geomorph)
library(EBImage)
library(shapes)
library(devtools)
install_github("vbonhomme/Momocs")
library(Momocs)

# FUNCTION:
Plot_ConvexHull<-function(xcoord, ycoord, lcolor){
  hpts <- chull(x = xcoord, y = ycoord)
  hpts <- c(hpts, hpts[1])
  polygon(xcoord[hpts], ycoord[hpts], col = lcolor, border = lcolor, lwd = 7)
  
}  
# END OF FUNCTION

#Only needs to be run once:
#Shell script
#ImageMagick
folderGold <- "/Users/gregorymatthews/Dropbox/brophyTeeth/DeciduousHomininGimpTest"
#Convert to .jpg
system(paste("cd ",folderGold,'; for i in *.png ; do convert "$i" "${i%.*}.jpg" ; done',sep=" "))
#Convert to gray scale
system(paste("cd ",folderGold,'; for i in *.jpg; do convert "$i" -colorspace Gray "$i"; done',sep=" "))
#Trim to center
system(paste("cd ",folderGold,'; for i in *.jpg; do convert -trim "$i" "$i"; done',sep=" "))

####################################################################
##Hominin files 
####################################################################
g<-1
#Get the file names of Juliet's b/w images.
goldFiles<-list.files("/Users/gregorymatthews/Dropbox/brophyTeeth/DeciduousHomininGimpTest/")
goldFiles <- goldFiles[grep("jpg",goldFiles)]

homList<-list()
toothList<-list()

for (g in 1:length(goldFiles)){print(g)
  imgID<-goldFiles[g]
  
  imgsGold<-import_jpg(paste("/Users/gregorymatthews/Dropbox/brophyTeeth/DeciduousHomininGimpTest/",imgID,sep=""))
  toothList[[g]]<-imgsGold
  plot(imgsGold[[1]],main=g)
    efa<-efourier(imgsGold,nb.h=15,norm=TRUE)
    tempGold<-c(imgID,unlist(as.numeric(c(efa$an,efa$bn,efa$cn,efa$dn))))
    homList[[g]] <- tempGold
    }


datHom<-data.frame(do.call(rbind,homList))
names(datHom)[1]<-"file"
names(datHom)[-1]<-paste0("h",1:60)



############################################################
##PC code
############################################################
x<-datHom
x<-x[complete.cases(x),]
x[,-1]<-apply(x[,-1],2,as.numeric)
#xScale<-apply(x[,-1],2,scale)
xScale<-x[,-1]
rownames(xScale)<-as.character(x[,1])
cor(xScale)
pc<-princomp(xScale)

col<-rgb(0.5,0.5,0.5,0.5)
#PC Plot
png("/Users/gregorymatthews/Dropbox/brophyTeeth/DeciduousHomininGimpTest/pcPlot_homininGimpTest_20160718.png",h=10,w=10,units="in",res=200)
plot(pc$scores[,1],pc$scores[,2],pch=16,col=col,cex=2,main="Principal Component Plot",sub="Test",xlab="PC1",ylab="PC2")

xrange <- range(pc$scores[,1])
yrange <- range(pc$scores[,2])

spec <-names(table(x$Species)[table(x$Species)>1])
for (sp in spec){
  temp<-pc$scores[x$Species==sp,1:2]
  # Plot it up!
  par(tck = 0.02, mgp = c(1.7, 0.3, 0))
  #plot(pc$scores[,1], pc$scores[,2], type = "p", pch = 1, col = "black", xlim = c(xrange), ylim = c(yrange))
  Plot_ConvexHull(xcoord = temp[,1], ycoord = temp[,2], lcolor = col[x$Species==sp])
}

text(pc$scores[,1],pc$scores[,2],(rownames(xScale)),pos=1,cex=.5)
dev.off()


