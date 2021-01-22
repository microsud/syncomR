# using mean rel abund
ps.rar.b5 <- rarefy_even_depth(ps1.b5)
#otu.rar.b5 <- abundances(ps.rar)
otu.tb <- taxa_time_table(ps.rar.b5, normalize = F, time.col= 'Time_hr', remove.zero = F)
active.max.growths <- maxgrowth(OTUs.active.persistent + .1)
mxg <- maxgrowth(t(otu.tb)+.1)
maxgrowth <- function(sbs){

  # create new matrix to calculate growth rates
  y.s <- matrix(ncol = ncol(sbs), nrow = nrow(sbs)-1)
  for(i in 1:(nrow(sbs)-1)){
    y.s[i,] = log(sbs[i+1,] / sbs[i,])
  }

  # for the gr's just calculated (ys) where is the max
  # for each OTU (i.e., each col) in the matrix
  max.growth.time <- apply(y.s, 2, which.max)

  # now, what is the gr at that max id'ed above
  max.growths <- apply(y.s, 2, max)

  # here, export a df of max gr and time of max gr for each OTU in the input df
  return(data.frame(
    rate = max.growths, sample = max.growth.time))
}
