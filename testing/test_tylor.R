library(syncomR)
data("SyncomFiltData")
ps1.b5 <- subset_samples(SyncomFiltData, StudyIdentifier== "Bioreactor A")
otu.b5.rel <- taxa_time_table(ps1.b5, normalize = TRUE,
                              time.col= 'Time_hr',
                              remove.zero = TRUE)
library(seqtime)
f5.taylor <- taylor(otu.b5.rel, type = "taylor", boot = 1000, interval = "prediction",
                    lower.conf = 0.025, upper.conf = 0.975, pseudo = 0.000000001,
                    col = "black", header = "Bioreactor A", label = F, plot = TRUE)

plot(logmeans,logvars, xlim=xlim, ylim=ylim,xlab="Log(mean)", ylab="Log(variance)", main=paste("Taylor's law",header,"\np-value:",round(pval,3),", adjusted R2:",round(sum$adj.r.squared,2),", slope:",round(linreg$coefficients[2],2)),type="p",pch="+",col=col)
abline(linreg,bty="n",col="red")
lines(reg.data$logmeans, interval[,2], col="blue", lty=2)
lines(reg.data$logmeans, interval[,3], col="blue", lty=2)


colnames(reg.data) <- c("logvars.reg", "logmeans.reg")
logs <- cbind(logmeans, logvars)
df.tay <- cbind(reg.data, logs,interval)

p <- ggplot(df.tay, aes(logmeans,logvars)) +
  geom_point() +
  stat_smooth(method = lm)
# 3. Add prediction intervals
p + geom_line(aes(y = lwr), color = "red", linetype = "dashed")+
  geom_line(aes(y = upr), color = "red", linetype = "dashed")
