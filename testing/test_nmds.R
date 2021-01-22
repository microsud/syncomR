plot_nmds_stability <- function(stab.in,
                                refphase = "grey60",
                                refstate = "red",
                                evolution.points = "yellow",
                                max.distant.states = "steelblue",
                                final.state = "black"){
  stab.in <- dat.stab
  numberOfGates <- stab.in$numberOfGates
  trefId <- stab.in$trefId
  referenceState <- stab.in$referenceState
  experimentStart <- stab.in$experimentStart
  experimentEnd <- stab.in$experimentEnd
  maxEuclideanId <- stab.in$maxEuclideanId
  maxCanberraId <- stab.in$maxCanberraId
  data <- stab.in$data


  mds.out <- metaMDS(rbind(data[,2:(numberOfGates+1)], referenceState), distance="bray", autotransform=FALSE, zerodist="add")
  nmd.dat <- data.frame(MDS1=mds.out$point[,1],MDS2=mds.out$point[,2])
    #mds.out$points[length(mds.out$points[,1]),1],mds.out$points[length(mds.out$points[,1]),2]

  #nmd.dat[length(nmd.dat$MDS1),1:2] #from
  #nmd.dat[length(nmd.dat$MDS1),2] #from
  #nmd.dat[(trefId+1),1:2] #to
  #nmd.dat[(trefId+1),2] #to
  #choose refstate to refstate+1 for start arrow
  arow <- rbind(nmd.dat[length(nmd.dat$MDS1),1:2], nmd.dat[(trefId+1),1:2])
  p <- ggplot(nmd.dat, aes(MDS1, MDS2)) + geom_point(shape= 21,size=4) + xlim(c(-0.5, 0.5)) +
    ylim(c(-0.5, 0.5)) + theme_syncom() + geom_hline(aes(yintercept = 0), alpha = 0.1) +
    geom_vline(aes(xintercept = 0), alpha = 0.1) +
    geom_path(data = nmd.dat[c((trefId+1):(length(nmd.dat[,1])-1)),],
              aes(MDS1,MDS2), arrow = arrow(angle = 15, length = unit(0.1, "inches"),
                                            ends = "last", type = "closed")) +
    geom_path(data = arow,
              aes(MDS1,MDS2), linejoin= "mitre", arrow = arrow(angle = 15, length = unit(0.1, "inches"),
                                            ends = "last", type = "closed"),size = 0.5)
  #refphase
  p <-p + geom_point(data=nmd.dat[c(1:trefId),], aes(MDS1,MDS2), fill = "grey60",size=4,  shape= 21)

  p#refstate
  p <- p + geom_point(data=nmd.dat[c(length(nmd.dat[,1]), length(nmd.dat[,1])),],
                      aes(MDS1,MDS2), fill = "red",size=4,  shape= 21)

  # evolution mds.out$points
  p <- p + geom_point(data = nmd.dat[c((trefId+1):(length(nmd.dat[,1])-1)),],
                      aes(MDS1,MDS2), fill = "yellow", shape= 21, size=4)

  #mark most distant states
  p <- p + geom_point(data = nmd.dat[c(maxEuclideanId, maxEuclideanId),],
                      aes(MDS1,MDS2), fill = "steelblue", shape= 24, size=4)

  p <- p + geom_point(data = nmd.dat[c(maxCanberraId, maxCanberraId),],
                      aes(MDS1,MDS2), fill = "brown3", shape= 24, size=4)
  # mark final state
  p <- p + geom_point(data = nmd.dat[c(length(nmd.dat[,1])-1, length(nmd.dat[,2])-1),],
                      aes(MDS1,MDS2), fill = "black", shape= 21, size=4)
  p <- p + labs(title = "",subtitle = expression(paste("Community evolve apart ", s[ref])))
  #p <- p + geom_path(nmd.dat, aes(MDS1,MDS2), colour="black")


  return(list("nmds_data"=mds.out,"plot" = p))
}




