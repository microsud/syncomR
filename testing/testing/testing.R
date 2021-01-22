
data <- dat.stab
dat.stab$trefId
data <- NULL
nmds.dt <- plot_nmds(dat.stab)#, time.col= 'Time_hr')

plot_nmds <- function(stab.out){
  titleSize = 0.95
  numberOfGates <- stab.out$numberOfGates
  trefId <- stab.out$trefId
  referenceState <- stab.out$referenceState
  experimentStart <- stab.out$experimentStart
  experimentEnd <- stab.out$experimentEnd
  maxEuclideanId <- stab.out$maxEuclideanId
  maxCanberraId <- stab.out$maxCanberraId
  data <- stab.out$data

  if (overrideMaxEuclidean == -1) {
    maxEuclideanId <- order(data$euclidean, decreasing=TRUE)[1]
  } else {
    maxEuclideanId <- match(overrideMaxEuclidean, data[[1]])
  }
  if (overrideMaxCanberra == -1) {
    maxCanberraId <- order(data$canberra, decreasing=TRUE)[1]
  } else {
    maxCanberraId <- match(overrideMaxCanberra, data[[1]])
  }

  mds.out <- metaMDS(rbind(data[,2:(numberOfGates+1)], referenceState), distance="bray", autotransform=FALSE, zerodist="add")
  plot(mds.out, type="n", main=expression('Community evolve apart s'[ref]), font.main=1, cex.main = titleSize)
  lines(mds.out$points[c((trefId+1):(length(mds.out$points[,1])-1)),], col = "black")
  # marking points in reference space
  if (trefId > 1) {
    points(mds.out$points[c(1:trefId),], col = "grey", pch = 21, bg="white")
  } else {
    points(mds.out$points[c(trefId, trefId),], col = "grey", pch = 21, bg="white")
  }
  # mark reference state
  points(mds.out$points[c(length(mds.out$points[,1]), length(mds.out$points[,1])),], col ="red", pch = 21, cex = 1.3, bg = "white")
  # mark evolution
  points(mds.out$points[c((trefId+1):(length(mds.out$points[,1])-1)),], bg = "white", pch = 21)
  # mark direction
  arrows(mds.out$points[length(mds.out$points[,1]),1],mds.out$points[length(mds.out$points[,1]),2],mds.out$points[(trefId+1),1],mds.out$points[(trefId+1),2], length = 0.08, angle = 8, lwd=1)
  # mark most distant states
  points(mds.out$points[c(maxEuclideanId, maxEuclideanId),], col = "deepskyblue3", pch = 24, cex = 1.2, bg = "white")
  points(mds.out$points[c(maxCanberraId, maxCanberraId),], col = "brown3", pch = 24, cex = 1.2, bg = "white")
  # mark final state
  points(mds.out$points[c(length(mds.out$points[,1])-1, length(mds.out$points[,1])-1),], col = "black", pch = 16, cex=1.3)

  return(mds.out)
}

head(nmds.dt$points)
nmd.dat <- data.frame(MDS1=nmds.dt$point[,1],MDS2=nmds.dt$point[,2])
head(nmd.dat)
#stems<-colSums(dune) #total abundances for each species
spps <- data.frame(scores(mds.out, display = "species")) #dataframe of species scoes for plotting
spps$species <- row.names(spps) # making a column with species names
#spps$colsums <- stems #adding the colSums from above
spps<-spps[!is.na(spps$NMDS1) & !is.na(spps$NMDS2),] #removes NAs

#spps.colmedian <- median(spps$colsums) #create an object that is the median of the abundance of the measured species
#spps.colmean <- mean(spps$colsums) #creates a mean instead if you wish to use
#spps2 <- subset(spps,spps$colsums > spps.colmean) #select the most abundant species. Could discard fewer by going something like - spps$colsums>(spps.colmedian/2) instead

spps$species <- factor(spps$species) #otherwise factor doesn't drop unused levels and it will throw an error
nmds.dt <- plot_nmds(dat.stab)
p <- ggplot(nmd.dat, aes(MDS1, MDS2)) + geom_point(shape= 21,size=4) + xlim(c(-0.5, 0.5)) +
  ylim(c(-0.5, 0.5)) + theme_classic() +
  geom_path(data = nmd.dat[c((trefId+1):(length(nmd.dat[,1])-1)),],
            aes(MDS1,MDS2), arrow = arrow())
#refphase
p <-p + geom_point(data=nmd.dat[c(1:trefId),], aes(MDS1,MDS2), fill = "grey60",size=4,  shape= 21)
p
#refstate
p <- p + geom_point(data=nmd.dat[c(length(nmd.dat[,1]), length(nmd.dat[,1])),],
                    aes(MDS1,MDS2), fill = "red",size=4,  shape= 21)
p
# evolution mds.out$points
p <- p + geom_point(data = nmd.dat[c((trefId+1):(length(nmd.dat[,1])-1)),],
                    aes(MDS1,MDS2), fill = "yellow", shape= 21, size=4)
p
#mark most distant states
p <- p + geom_point(data = nmd.dat[c(maxEuclideanId, maxEuclideanId),],
                    aes(MDS1,MDS2), fill = "steelblue", shape= 24, size=4)
p
p <- p + geom_point(data = nmd.dat[c(maxCanberraId, maxCanberraId),],
                    aes(MDS1,MDS2), fill = "brown3", shape= 24, size=4)
# mark final state
p <- p + geom_point(data = nmd.dat[c(length(nmd.dat[,1])-1, length(nmd.dat[,2])-1),],
                    aes(MDS1,MDS2), fill = "black", shape= 21, size=4)
#p <- p + geom_path(nmd.dat, aes(MDS1,MDS2), colour="black")

p




points(mds.out$points[c((trefId+1):(length(mds.out$points[,1])-1)),], bg = "white", pch = 21)
# mark direction
arrows(mds.out$points[length(mds.out$points[,1]),1],mds.out$points[length(mds.out$points[,1]),2],mds.out$points[(trefId+1),1],mds.out$points[(trefId+1),2], length = 0.08, angle = 8, lwd=1)
# mark most distant states
points(mds.out$points[c(maxEuclideanId, maxEuclideanId),], col = "deepskyblue3", pch = 24, cex = 1.2, bg = "white")
points(mds.out$points[c(maxCanberraId, maxCanberraId),], col = "brown3", pch = 24, cex = 1.2, bg = "white")
# mark final state
points(mds.out$points[c(length(mds.out$points[,1])-1, length(mds.out$points[,1])-1),], col = "black", pch = 16, cex=1.3)

data <- as.data.frame(nmds.dt$points[c((trefId+1):(length(nmds.dt$points[,1])-1)),])
head(data)
# points
nmds.dt$points[c((trefId+1):(length(nmds.dt$points[,1])-1)),]
# lines
nmds.dt$points[c((trefId+1):(length(nmds.dt$points[,1])-1)),]
nmds.dt$points[c(1:trefId),]

nmds.dt$points[c(trefId, trefId),]
# mark reference state
mds.out$points[c(length(mds.out$points[,1]), length(mds.out$points[,1])),]
# mark evolution
mds.out$points[c((trefId+1):(length(mds.out$points[,1])-1)),]
# mark direction
#mds.out$points[length(mds.out$points[,1]),1],mds.out$points[length(mds.out$points[,1]),2],mds.out$points[(trefId+1),1],mds.out$points[(trefId+1),2]
# mark most distant states
mds.out$points[c(maxEuclideanId, maxEuclideanId),]
mds.out$points[c(maxCanberraId, maxCanberraId),]
# mark final state
mds.out$points[c(length(mds.out$points[,1])-1, length(mds.out$points[,1])-1),]


if (trefId > 1) {
  points(mds.out$points[c(1:trefId),], col = "grey", pch = 21, bg="white")
} else {
  points(mds.out$points[c(trefId, trefId),], col = "grey", pch = 21, bg="white")
}


NMDS=data.frame(x=nmds.dt$point[,1],y=nmds.dt$point[,2],Country=as.factor(grouping_info[,1]),Latrine=as.factor(grouping_info[,2]),Depth=as.factor(grouping_info[,3]))



