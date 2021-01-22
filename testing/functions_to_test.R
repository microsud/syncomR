library(progress)
# using mean rel abund
ps.rar.b5 <- rarefy_even_depth(ps1.b5)
#otu.rar.b5 <- abundances(ps.rar)
otu.tb <- taxa_time_table(ps.rar.b5, normalize = F, time.col= 'Time_hr', remove.zero = F)
OTUs.total.persistent <- t(otu.tb)
OTUs.total.nozero <- replace(OTUs.total.persistent, which(OTUs.total.persistent == 0), .5)
nfd.total <- calc.neg.freq.dep(OTUs.total.nozero)

nfd.total %>% filter(-fd > 0) -> nfd.total
nfd.total.cov <- cov(log(nfd.total$eq.freq), log(-nfd.total$fd))
nfd.total %>% ggplot(aes(x = eq.freq, y = -fd)) +
  geom_point() +
  scale_y_log10() +
  scale_x_log10() +
  geom_smooth(method = "lm") +
  ylab("Negative Frequency Dependence") +
  xlab("Equilibrium Frequency")

# Generate null distributions
iternum <- 50
tot.nd <- get.null.distr(comm = OTUs.total.nozero, iters = iternum)
#act.nd <- get.null.distr(comm = OTUs.active.nozero, iters = iternum)
# Compare observed with null
hist(tot.nd, breaks = 20); abline(v = nfd.total.cov, col = "red")
(tot.nfd.pval <- rank(c(nfd.total.cov, tot.nd))[1] / length(c(nfd.total.cov, tot.nd)))
nfd.obs <- rbind.data.frame(cbind.data.frame(subset = "Total", nfd = nfd.total.cov))
nfd.nulldists <- rbind.data.frame(cbind.data.frame(subset = "Total", nfd = tot.nd))
nfd.nulldists %>%
  ggplot(aes(y = nfd, x = subset)) +
  #geom_jitter(color = "grey", alpha = 0.1) +
  geom_violin(alpha = 0.5, show.legend = F, fill = "grey") +
  geom_point(data = nfd.obs, aes(x = subset, y = nfd),
             color = "black", size = 5, show.legend = F) +
  scale_fill_manual(values = strain.colors[c(1,4)]) +
  coord_flip() +
  ylab("Covariance between Equilibrium Frequency\n and Negative Frequency Dependence") +
  xlab("")

# using mean rel abund
calc.neg.freq.dep <- function(sbs){
  y.s <- matrix(ncol = ncol(sbs), nrow = nrow(sbs)-1)
  for(i in 1:(nrow(sbs)-1)){
    y.s[i,] = log(sbs[i+1,] / sbs[i,])
  }
  eq.freq <- vector(length = ncol(y.s))
  fd <- vector(length = ncol(y.s))
  sbs.rel <- decostand(sbs[-nrow(sbs),], method = "total")
  for(i in 1:ncol(y.s)){
    mod <- lm(y.s[,i] ~ sbs.rel[,i])
    eq.freq[i] <- mean(sbs.rel[,i])
    fd[i] <- coef(mod)[2]
  }
  # plot(log10(eq.freq), log10(-fd),
  #      ylab = "log10 Negative Freq. Dependence",
  #      xlab = "log10 Equilibrium Frequency")
  return(na.omit(cbind.data.frame(eq.freq, fd)))
}
shuffle.ts <- function(species.ts){
  species.ts[sample(1:length(species.ts))]
}
get.null.distr <- function(comm, iters = 5000){
  pb <- progress_bar$new(total = iters)
  nfd.cov <- rep(NA, iters)
  for(i in 1:iters){
    nullcom <- apply(comm, MARGIN = 2, FUN = shuffle.ts) # randomize sampling order
    nfd <- calc.neg.freq.dep(nullcom)
    nfd <- nfd %>% filter(-fd > 0)
    nfd.cov[i] <- cov(log(nfd$eq.freq), log(-nfd$fd))
    pb$tick()
  }
  return(nfd.cov)
}

