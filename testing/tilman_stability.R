
library(syncomR)
data("SyncomGMM")
head(SyncomGMM)
#A Clark, 14 November 2013
#Example of Lehman/Tilman stability
#Not error checked - make sure this makes sense before using it

#Make fake data
dat<-data.frame(plot=rep(1:10, each=20), species=letters[1:5], times=rep(1:4, each=5), biomass=runif(200))
#plot is plot ID, species is species ID, times is sampling date (or ID), biomass is biomass
#NOTE - if a species has "zero" abundance, you will need to mark it as such.
#I.e., number of samples within a particular plot must be equal for all species
ps.df <- microbiomeutilities::phy_to_ldf(SyncomFiltData, "compositional")
colnames(ps.df)
ps.df.sub <- ps.df[,c("StudyIdentifier", "Species", "Time_hr_num", "Abundance")]
colnames(ps.df.sub) <- c("plot", "species", "times", "biomass")
dat <- ps.df.sub
#sum of biomass per plot
sum_abund<-rowSums(tapply(dat$biomass, list(dat$plot, dat$species), mean))

#sum of variance for each species in each plot
sum_var<-rowSums(tapply(dat$biomass, list(dat$plot, dat$species), var))

#function to sum of covariance for each possible combination of species in each plot
find_covar_sum<-function(x) {
  nspecies<-length(unique(x$species))
  nsamples<-length(unique(x$times))

  comparisons<-outer(1:nspecies, 1:nspecies, ">")


  #Make a matrix with species biomass for the plot in each column, and rows ordered by sampling event
  sp_mat<-matrix(nrow=nsamples, ncol=nspecies, data=x$biomass[order(paste(x$species, x$times))])

  return(var(sp_mat)[comparisons]) #Record only unique, non-diagonal values
}

#Run function for each plot
nsites<-length(unique(dat$plot))
sum_covar<-numeric(nsites)
n<-1
for(i in sort(unique(dat$plot))) {
  x<-dat[dat$plot==i,]
  sum_covar[n]<-sum(find_covar_sum(x))
  names(sum_covar)[n]<-i
  n<-n+1
}

#Formula: by plot, sum(biomass of species)/sqrt(sum(var(biomass of species))+sum(cov(biomass of species)))
#Sensu Lehman and Tilman, 2000. American Naturalist, Vol. 156, No. 5 (November 2000), pp. 534-552
Lehman_Tilman_Stability<-sum_abund/sqrt(sum_var+sum_covar)

Lehman_Tilman_Stability
