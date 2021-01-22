theme_set(theme_bw())
min_theme <- theme_update(
  #panel.border = element_blank(),
  panel.grid = element_blank(),
  panel.spacing = unit(0, "line"),
  #axis.ticks = element_blank(),
  legend.title = element_text(size = 16),
  legend.text = element_text(size = 14),
  axis.text = element_text(size = 14),
  axis.title = element_text(size = 14),
  strip.background = element_blank(),
  strip.text = element_text(size = 16),
  legend.key = element_blank()
)

strain.colors<-c(Akkermansia_muciniphila ='#e6194b', Bifidobacterium_adolescentis = '#3cb44b',
                 Collinsella_aerofaciens= '#ffe119', Bacteroides_ovatus ='#4363d8',
                 Bacteroides_xylanisolvens='#f58231', Agathobacter_rectalis ='#911eb4',
                 Anaerobutyricum_soehngenii ='#46f0f0', Eubacterium_siraeum ='#f032e6',
                 Blautia_hydrogenotrophica='#bcf60c', Coprococcus_catus='#fabebe',
                 Flavonifractor_plautii ='#008080', Roseburia_intestinalis ='#e6beff',
                 Faecalibacterium_prausnitzii ='#9a6324', Blautia_obeum= '#e4cd05',
                 Ruminococcus_bromii='#800000', Subdoligranulum_variabile='#aaffc3',
                 Acetate = '#808000', Butyrate = '#ffd8b1',
                 Formate = '#000075', Propionate = '#808080',
                 Lactate = '#ffffff', Succinate = '#000000')

library(syncomR)
library(ggplot2)
library(phyloseq)
library(microbiome)
# format_ps
data(SyncomRawCounts)
pseq <- format_ps(SyncomRawCounts, tax.level = "Species")

data(SyncomRawCounts)
copy_num <- read.table("inst/extras/mdbmm_copy_numbers.txt", header = T,
                      row.names = 1, stringsAsFactors = F, sep = "\t")

# normalize_copy_number
ps_cp_normalized <- normalize_copy_number(pseq,
                                          column_with_ids = "Genus_species",
                                          copy_num_tab = copy_num)

data(SyncomFiltData)
fasting_cols <- c("#b3de69", "#fb8072", "#80b1d3")

pl <- taxa_coverage(SyncomFiltData, coverage = 0.9999,
                    time.column = "Time_hr_num",
                    shape.variable = "Acetate_Feed",
                    color.variable = "Fasting",
                    #facet.var = "StudyIdentifier",
                    color.pal = fasting_cols,
                    y.breaks = seq(0, 16, 1),
                    y.limits = c(11, 16))
pl + facet_wrap(~ StudyIdentifier) + theme_syncom()


data(SyncomFiltData)
ps1.b5 <- subset_samples(SyncomFiltData, StudyIdentifier== "Bioreactor A")
#remove_T80 <- c("Ferm_1_5_80", "Ferm_1_6_80", "Ferm_1_8_80")
#ps1.sub <- prune_samples( !(sample_names(ps1.b5) %in% remove_T80), ps1.b5)
otu.tb <- taxa_time_table(ps1.b5, normalize = TRUE, time.col= 'Time_hr', remove.zero = TRUE)
head(otu.tb)



p <- temporal_turnover(ps1.b5, tree =NULL,time.col = 'Time_hr_num',
                       method = "canberra", compositional = TRUE)
p + theme_syncom()


fer_cols <- c(`Bioreactor A`= "#b2182b", `Bioreactor B`="#2166ac", `Bioreactor C` = "#35978f")
px <- temporal_diversity(SyncomFiltData,time.col = 'Time_hr_num', div.measure = "gini")
px + facet_grid(~StudyIdentifier) + theme_syncom()



data(SyncomFiltData)
SyncomFiltData.rel <- microbiome::transform(SyncomFiltData, "compositional")
SyncomFiltData.rel <- subset_samples(SyncomFiltData.rel, Time_hr_num <=200)
p <- plot_syncom_composition(SyncomFiltData.rel,
                         type = "bar",
                         time.col = "Time_hr_num",
                         taxa.level = "OTU",
                         sp.fill.pal = syncom_colors("BacterialSpecies"), facet.var = "StudyIdentifier")
print(p+theme_syncom(base_size = 6))

library(syncomR)
data(SyncomFiltData)
ps1.b5 <- subset_samples(SyncomFiltData, StudyIdentifier== "Bioreactor A")
#ps1.sub <- prune_samples( !(sample_names(ps1.b5) %in% remove_T80), ps1.b5)
#ps1.sub <- subset_samples(ps1.b5, Time_hr_num >= 140)

dat.stab <- stability_properties(ps1.b5, time.col = "Time_hr",
                                 experimentStart = 28,
                                 #experimentEnd = 460,
                                 tref = 152)
p_ba <- plot_nmds_stability(dat.stab)
p_ba$plot + theme_syncom() + xlab("NMDS1") + ylab("NMDS2")

#head(dat.stab)

plot_stability_properties(dat.stab, property = "dist.ot")




data(SyncomFiltData)
remove_T80 <- c("Ferm_1_5_80", "Ferm_1_6_80", "Ferm_1_8_80")
SyncomFiltData <- prune_samples( !(sample_names(SyncomFiltData) %in% remove_T80), SyncomFiltData)
fer_cols <- c(`Bioreactor A`= "#b2182b", `Bioreactor B`="#2166ac", `Bioreactor C` = "#35978f")
p<- plot_trajectory(SyncomFiltData, time.col = "Time_hr", taxa = "Akkermansia_muciniphila",
                    type = "Species", group.variable = "StudyIdentifier",
                    color.pal = fer_cols)
p + theme_classic()


p<- plot_trajectory(SyncomFiltData, time.col = "Time_hr",
                    type = "all", group.variable = "StudyIdentifier",
                    color.pal = fer_cols)
p + theme_classic()

data(SyncomGMM)
library(ggplot2)
p<- plot_modul_abundances(SyncomGMM,
                          tax.variable = c("Akkermansia_muciniphila", "Bacteroides_xylanisolvens", "Bacteroides_ovatus"),
                     mm.variable = c("starch degradation", "mucin degradation"),
                     color.pal = strain.colors,
                     nrow=2, ncol = 1)
p + theme_classic()

p<- plot_module_abundances(SyncomGMM,
                          tax.variable = c("Akkermansia_muciniphila", "Bacteroides_xylanisolvens", "Bacteroides_ovatus"),
                          mm.variable = "alanine degradation I",
                          color.pal = strain.colors,
                          nrow=2, ncol = 1)
p + theme_classic()

data(SyncomFiltData)
library(microbiome)
hd_tax <- find_hyperdominant_taxa(SyncomFiltData, cutoff= 0.5)
hd_tax


library(syncomR)
ps1.b5 <- subset_samples(SyncomFiltData, StudyIdentifier== "Bioreactor A")
#otu.tb <- taxa_time_table(ps1.b5, normalize = T, time.col= 'Time_hr', remove.zero = T)

p.a <- plot_growth_rate(ps1.b5, rarefy = TRUE, depth=min(sample_sums(ps1.b5)))
print(p.a + scale_color_manual(values = strain.colors)+ theme_syncom())

ps1.b6 <- subset_samples(SyncomFiltData, StudyIdentifier== "Bioreactor B")
#otu.tb <- taxa_time_table(ps1.b5, normalize = T, time.col= 'Time_hr', remove.zero = T)

p.b <- plot_growth_rate(ps1.b6, rarefy = TRUE, depth=min(sample_sums(ps1.b6)))
print(p.b + theme_syncom())

OTUs.total.nozero <- replace(otu.tb, which(otu.tb == 0), .5)
h <-calc.neg.freq.dep(otu.tb)

################

data(SyncomFiltData)
ps1.b5 <- subset_samples(SyncomFiltData, StudyIdentifier== "Bioreactor A")
otu.tb <- taxa_time_table(ps1.b5, normalize = T, time.col= 'Time_hr', remove.zero = TRUE)
head(otu.tb)

Density <- otu.tb[2,]  # create a sequence of numbers from 0 to 500, representing a range of population densities
summary(Density)
## CONSTANTS

b_max <- 0.5   # maximum reproduction (at low densities)
d_min <- 0.001  # minimum mortality (at low densities)

a <- 0.00005    # D-D terms
c <- 0.00001

b <- b_max - a*Density
d <- d_min + c*Density

plot(Density,b,type="l",col="green",lwd=2,main="Density-Dependence!")
points(Density,d,type="l",col="red",lwd=2)
legend("topleft",col=c("green","red"),lty=c(1,1),legend=c("per-capita birth rate","per-capita mortality"),bty="n")




syncom.rel <- SyncomFiltData
det <- c(0, 0.1, 0.5, 2, 5, 20)/100
prevalences <- seq(.05, 1, .05)
#ggplot(d) + geom_point(aes(x, y)) + scale_x_continuous(trans="log10", limits=c(NA,1))

plot_core(syncom.rel, prevalences = prevalences, detections = det, plot.type = "lineplot") + xlab("Relative Abundance (%)")

taxa_coverage



ps <- syncom.rel
psbA <- subset_samples(ps, StudyIdentifier =="Bioreactor A")
abund = taxa_time_table(psbA, normalize = F, time.col = "Time_hr_num", remove.zero = F)
#samples = sample_data(ps)
cutoff <- seq(0, 1, 0.1)
mymat = data.frame(matrix(nrow=nrow(abund))) # create a 10 x 10 matrix (of 10 rows and 10 columns)
dim(mymat)
for(j in cutoff) # for each column
  {
    diversity = as.data.frame(coverage(t(abund), j))
    mymat <- cbind(mymat, diversity) # assign values based on position: product of two indexes
    #mymat <- mymat[,-1]
    #colnames(mymat) <- c(as.character(cutoff))
}
mymat <- mymat[,-1]
colnames(mymat) <- c(as.character(cutoff))
mymat$Time <- rownames(mymat)
head(mymat)
mymat.m <- reshape2::melt(mymat)
mymat.m$Time <- as.numeric(mymat.m$Time)
#mymat.m$variable <- as.numeric(mymat.m$variable)
#mymat.m$value <- as.numeric(mymat.m$value)
colnames(mymat.m) <- c("Time", "Threshold", "Coverage")
head(mymat.m)
library(RColorBrewer)
p <- ggplot(mymat.m, aes(x = Time, y = Threshold))
p <- p + geom_tile(aes(fill=Coverage), size=1)
p  <- p + theme_syncom() + scale_fill_gradientn(colours = brewer.pal(9,"Spectral"))
p
p <- p + xlab()
p <- p + ylab()
p <- p + ggtitle()
list(plot = p, data = x)



df.m <- merge(meta(ps), mymat)
head(df.m)


df <- melt(df, "ID")
names(df) <- c("Abundance", "Prevalence", "Count")
df$Abundance <- as.numeric(as.character(df$Abundance))
df$Prevalence <- as.numeric(as.character(df$Prevalence))
df$Count <- as.numeric(as.character(df$Count))
p <- ggplot(mymat.m, aes(x = variable, y = value, #color = Prevalence,
                    group = variable))
p <- p + geom_line()
p <- p + geom_point()
p <- p + scale_x_log10()
p <- p + xlab()
p <- p + ylab()
p <- p + ggtitle()
list(plot = p, data = x)

head(df)

library(syncomR)
data("SynComRNAKEGG")
kegg <- SynComRNAKEGG

head(SynComRNAKEGG)
kegdf <- group_by_variable(SynComRNAKEGG, by = "KOID")
head(kegdf)
class(kegdf)

kegdf[1:5,1:5]




data(SyncomFiltData)
ps1.b5 <- subset_samples(SyncomFiltData, StudyIdentifier == "Bioreactor A")
ps1.sub <- subset_samples(ps1.b5, Time_hr_num >= 28)
dat.stab <- stability_properties(ps1.sub, time.col = "Time_hr", tref = 152)
plot_resilience(dat.stab, method = "euclidean")
plot_resilience(dat.stab, method = "canberra")











