library(vegan)
library(geomorph)
library(shapes)
library(devtools)
#install_github("vbonhomme/Momocs")
library(Momocs)

#I have slightly modified these function to prevent rotation.  
efourier.Out <- function(x, nb.h, smooth.it = 0, norm = TRUE, start = FALSE, verbose=TRUE, ...) {
  Out <- x
  # validates
  Out %<>% validate()
  q <- floor(min(sapply(Out$coo, nrow)/2))
  if (missing(nb.h)) {
    #nb.h <- ifelse(q >= 32, 32, q)
    nb.h <- calibrate_harmonicpower(Out, thresh = 99, verbose=FALSE, plot=FALSE)$minh
    if (verbose) message("'nb.h' not provided and set to ", nb.h, " (99% harmonic power)")
  }
  if (nb.h > q) {
    nb.h <- q  # should not be 1 #todo
    message("at least one outline has no more than ", q * 2,
            " coordinates. 'nb.h' has been set to ", q, " harmonics")
  }
  coo <- Out$coo
  col.n <- paste0(rep(LETTERS[1:4], each = nb.h), rep(1:nb.h,
                                                      times = 4))
  coe <- matrix(ncol = 4 * nb.h, nrow = length(coo), dimnames = list(names(coo),
                                                                     col.n))
  for (i in seq(along = coo)) {
    # todo: vectorize ?
    ef <- efourier(coo[[i]], nb.h = nb.h, smooth.it = smooth.it,
                   verbose = TRUE)
    if (norm) {
      ef <- efourier_norm(ef, start = start)
      coe[i, ] <- c(ef$an, ef$bn, ef$cn, ef$dn)
    }
  }
  
  coe[abs(coe) < 1e-12] <- 0  #not elegant but round normalized values to 0
  res <- OutCoe(coe = coe, fac = Out$fac, method = "efourier", norm = norm)
  res$cuts <- ncol(res$coe)
  return(res)
}


efourier_norm <- function(ef, start = FALSE) {
  A1 <- ef$an[1]
  B1 <- ef$bn[1]
  C1 <- ef$cn[1]
  D1 <- ef$dn[1]
  nb.h <- length(ef$an)
  theta <- 0.5 * atan(2 * (A1 * B1 + C1 * D1)/(A1^2 + C1^2 -
                                                 B1^2 - D1^2))%%pi
  phaseshift <- matrix(c(cos(theta), sin(theta), -sin(theta),
                         cos(theta)), 2, 2)
  M2 <- matrix(c(A1, C1, B1, D1), 2, 2) #%*% phaseshift
  v <- apply(M2^2, 2, sum)
  if (v[1] < v[2]) {
    theta <- theta + pi/2
  }
  theta <- (theta + pi/2)%%pi - pi/2
  Aa <- A1 * cos(theta) + B1 * sin(theta)
  Cc <- C1 * cos(theta) + D1 * sin(theta)
  scale <- sqrt(Aa^2 + Cc^2)
  psi <- atan(Cc/Aa)%%pi
  if (Aa < 0) {
    psi <- psi + pi
  }
  size <- 1/scale
  rotation <- matrix(c(cos(psi), -sin(psi), sin(psi), cos(psi)),
                     2, 2)
  A <- B <- C <- D <- numeric(nb.h)
  if (start) {
    theta <- 0
  }
  for (i in 1:nb.h) {
    mat <- size *
      matrix(c(ef$an[i], ef$cn[i],
               ef$bn[i], ef$dn[i]), 2, 2)
    A[i] <- mat[1, 1]
    B[i] <- mat[1, 2]
    C[i] <- mat[2, 1]
    D[i] <- mat[2, 2]
    lnef <- c(A[i], B[i], C[i], D[i])
  }
  list(an = A, bn = B, cn = C, dn = D, size = scale, theta = theta,
       psi = psi, ao = ef$ao, co = ef$co, lnef = lnef)
}

# path<-"/Users/gregorymatthews/Dropbox/brophyTeeth/DeciduousTeethProject/LowerFirstMolar/"
# file1<-"Ldm1.csv"
# toothname<-"LDM1"

##################################################
#This function reads in all the teeth and computes
#the average of each group of teeth
##################################################

deciduous_teeth <- function(path, file1, toothname){
#Only needs to be run once:
#Shell script
#ImageMagick
#folderGold <- "/Users/gregorymatthews/Dropbox/brophyTeeth/DeciduousTeethProject/"
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
  
#Read in the key file.  
key <- read.csv(paste(path, file1, sep=""))[-1,1:3]
names(key) <- c("SpecimenNumber", "ID", "Group")
key$ID <- as.character(key$ID)
key$ID <- gsub(" ", "", key$ID)
key$ID <- gsub("_", "", key$ID)
key$ID <- toupper(key$ID)

####################################################################
##Reading in files traced by the expert (Juliet Brophy)
####################################################################
g <- 1
#Get the file names of Juliet's b/w images.
goldFiles <- list.files(path)
goldFiles <- goldFiles[grep("jpg", goldFiles)]
#goldFiles <- goldFiles[goldFiles%in%subsetFiles$X]

#I get an error with this tooth: "DSCN2879"
#I have no idea why. So I removed it for now.
#goldFiles <- goldFiles[-grep("DSCN2879",goldFiles)]
imgsList <- list()
for (g in c(1:length(goldFiles))){
  #Read in the files from Juliet
  imgID <- goldFiles[g]
  id <- strsplit(imgID, "[.]")[[1]][1]
  #Contains the image done by Juliet
  imgsList[[imgID]] <- import_jpg(paste(path, imgID,sep=""))[[1]]
}

outCoo <- Out(imgsList)
test <- efourier(outCoo, nb.h=20, norm=TRUE)


#Finding average shape by group
numPoints <- 150
ptsArray <- array(dim=c(numPoints, 2, dim(test$coe)[1]))
for (j in 1:dim(test$coe)[1]){
  a <- grep("A", colnames(test$coe))
  b <- grep("B", colnames(test$coe))
  c <- grep("C", colnames(test$coe))
  d <- grep("D", colnames(test$coe))
  ptsArray[ , ,j] <- efourier_i(list(an = test$coe[j,a], 
                                     bn = test$coe[j,b], 
                                     cn = test$coe[j,c],
                                     dn = test$coe[j,d],
                                     a0 = 0,
                                     c0 = 0),
                                     nb.pts = numPoints)
  
}

ptsArrayNames <- gsub(" ", "", toupper(names(test)))
ptsArrayNames <- gsub("_", "", ptsArrayNames)



png(paste0("/Users/gregorymatthews/Dropbox/brophyTeeth/DeciduousTeethProject/", toothname, "_meanShapeBySpecies.jpg"), res=300, units="in", w=10, h=4)
grps <- sort(as.character(unique(key$Group)))
grps <- grps[grps != ""]
par(mfrow = c(2, 5))
par(mar = c(1, 1, 1, 1))
for (g in grps){
  ids <- which(ptsArrayNames%in%key$ID[key$Group == g])
  if(length(ids) > 1){
    gpa<-procGPA(ptsArray[ , ,ids])
    plot(gpa$mshape, main = g , xlab = "", ylab = "", type = "l", 
         xlim = c(-1.5, 1.5, ylim = c(-1.5, 1.5), frame.plot = "F", xaxt = "n", 
         yaxt = "n", sub = paste0("n=",length(ids)), asp = 1)
  }
  if(length(ids) == 1){
    plot(ptsArray[ , ,ids], main = g, xlab = "", ylab = "", type = "l", 
         xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5), frame.plot = "F", xaxt = "n", yaxt = "n",
         sub = "n=1", asp = 1)
  }
}
dev.off()

}


# path<-"/Users/gregorymatthews/Dropbox/brophyTeeth/DeciduousTeethProject/LowerFirstMolar/"
# file1<-"Ldm1.csv"
# toothname<-"LDM1"
####################################################
#Now run each different tooth
####################################################
deciduous_teeth("~/Dropbox/brophyTeeth/DeciduousTeethProject/LowerFirstMolar/", "Ldm1.csv", "LDM1")
deciduous_teeth("~/Dropbox/brophyTeeth/DeciduousTeethProject/LowerSecondMolar/", "Ldm2.csv", "LDM1")
deciduous_teeth("~/Dropbox/brophyTeeth/DeciduousTeethProject/LowerThirdMolar/", "Ldm3.csv", "LDM1")
deciduous_teeth("~/Dropbox/brophyTeeth/DeciduousTeethProject/UpperFirstMolar/", "Udm1.csv", "UDM1")
deciduous_teeth("~/Dropbox/brophyTeeth/DeciduousTeethProject/UpperSecondMolar/", "Udm2.csv", "UDM2")
deciduous_teeth("~/Dropbox/brophyTeeth/DeciduousTeethProject/UpperThirdMolar/", "Udm3.csv", "UDM3")




