# install.packages("extendedForest_1.6.2.1.tar.gz", repos = NULL, type = "source")
# install.packages("gradientForest_0.1-37.tar.gz", repos = NULL, type = "source")
# install.packages('clv')
# install.packages("RDAforest_2.6.9.tar.gz")

library(RDAforest)
library(vegan)
library(dplyr)
library(clv)
library(extendedForest)
library(gradientForest)
library(rnaturalearth)
library(rnaturalearthdata)
library(terra)
library(viridis)
library(data.table)
library(RColorBrewer)
library(ggplot2)

################## GETTING THE DATA IN ##################
geno.df <- fread("../../03.GEA/BANS.ds6x.pass-maf-0.05.SNP.above4x.nomiss.raw_noheader.raw")
geno <- geno.df[,c(-1:-6)]
geno <- as.matrix(geno)
rownames(geno) <- geno.df$IID


env.all <- fread("../../03.GEA/data/BANS.initial.worldclim.txt")
head(env.all)
env <- env.all[,7:ncol(env.all)]
head(env)

latlong <- env.all %>% dplyr::select(Lat, Long)

envc <- read.table("data/BANS_envc.tsv")
envf <- read.table("data/BANS_envf.tsv")

ecotype <- read.table("../../03.GEA/data/BANS_sample_list_pop.tsv", header=F) %>% dplyr::select(V1, V3) %>% rename(ind=V1, cluster=V3)

################## DATA EXPLORE ##################

# distances between individuals (hence the t function, to transponse)
cordist=1-cor(t(geno))
# ordination:
ord=capscale(cordist~1)

# write the variance explained plot
png("results/CA_variance_explained.png", width = 800, height = 600, res = 120)

plot(ord$CA$eig / sum(ord$CA$eig),
     xlab = "PC",
     ylab = "Proportion of variance explained",
     type = "b",      # optional: adds points connected by lines
     pch = 16,        # solid points
     col = "steelblue")

dev.off()

# plot population structure
so=data.frame(scores(ord,scaling=1,display="sites"))

so$cluster <- ecotype$cluster

region_order <- c( "Cluster_1", "Cluster_5", "Cluster_4", "Cluster_8", "Cluster_12", "Cluster_2", "Cluster_3", "Cluster_7", "Cluster_11", "Cluster_14", "Cluster_13", "Cluster_6", "Cluster_9", "Cluster_10" )

# Make sure cluster is a factor in the desired order
so$cluster <- factor(so$cluster, levels = region_order)

# Define clusters by region
West_clusters    <- c("Cluster_1", "Cluster_5")
Midwest_clusters <- c("Cluster_4", "Cluster_8", "Cluster_12", "Cluster_2", "Cluster_3", "Cluster_7")
East_clusters    <- c("Cluster_11", "Cluster_14", "Cluster_13", "Cluster_6", "Cluster_9", "Cluster_10")

# Assign colors by cluster
W.cols <- setNames(brewer.pal(n = length(West_clusters), "Reds"), West_clusters)
Midwest.cols <- setNames(brewer.pal(n = length(Midwest_clusters), "Blues"), Midwest_clusters)
E.cols <- setNames(brewer.pal(n = length(East_clusters), "Greens"), East_clusters)

# Combine all colors
cluster_colors <- c(W.cols, Midwest.cols, E.cols)

png("results/pop_struc.png", width = 800, height = 600, res = 120)
ggplot(so, aes(MDS1, MDS2, color = cluster)) +
  geom_point(size = 2) +
  coord_equal() +
  theme_bw() +
  scale_color_manual(values = cluster_colors)
dev.off()

# rsync -avzP ericacnr@colostate.edu@login.rc.colorado.edu:/scratch/alpine/ericacnr@colostate.edu/BANS/06.genomic_offset/RDAForest_offset/results/*.png \
# /Users/ericarobertson/Desktop/BANS_adaptive_units/analysis/06.genomic_offset/results

################## CHECKING IBD ################## 

GCD=gcd.dist(latlong)
latlon.gcd=GCD[[1]]
distGCD=GCD[[2]]

png("results/IBD_signal.png", width = 800, height = 600, res = 120)

plot(as.dist(cordist)~distGCD,pch=16,cex=0.6,col=rgb(0,0,0,alpha=0.2))

dev.off()

protest(capscale(distGCD~1),capscale(cordist~1))

## removing IBD signal
latlon.gcd=GCD[[1]]
ord1=capscale(cordist~1+Condition(as.matrix(latlon.gcd)))

################## CLEANING PREDICTORS ##################
pc=hclust(as.dist(1-cor(env)))

png("results/predictors_cor.png", width = 800, height = 600, res = 120)
plot(pc)
abline(h=0.1, col="red")
dev.off()

rsync -avzP ericacnr@colostate.edu@login.rc.colorado.edu:/scratch/alpine/ericacnr@colostate.edu/BANS/06.genomic_offset/RDAForest_offset/results/predictors_cor.png \
/Users/ericarobertson/Desktop/BANS_adaptive_units/analysis/06.genomic_offset/results

################## EXPLORATORY RDA ################## 
gf=makeGF(ord1,env,pcs2keep=c(1:40))

gf$result

eigen.var=(ord1$CA$eig/sum(ord1$CA$eig))[names(gf$result)]
# total variance explained by model
sum(eigen.var*gf$result)

tokeep=15
# computing properly scaled importances:
imps=data.frame(importance_RDAforest(gf,ord1))
# some data frame housekeeping...
names(imps)="R2"
imps$var=row.names(imps)
# reordering predictors by their importances:
imps$var=factor(imps$var,levels=imps$var[order(imps$R2)])

# plotting
png("results/predictor_importance.png", width = 800, height = 600, res = 120)

ggplot(imps,aes(var,R2))+geom_bar(stat="identity")+coord_flip()+theme_bw()
dev.off()

rsync -avzP ericacnr@colostate.edu@login.rc.colorado.edu:/scratch/alpine/ericacnr@colostate.edu/BANS/06.genomic_offset/RDAForest_offset/results/predictor_importance.png \
/Users/ericarobertson/Desktop/BANS_adaptive_units/analysis/06.genomic_offset/results
################## TURNOVER CURVES ################## 
png("results/turnover_curves.png")
plot_gf_turnovers(gf,imps$var[1:6])
dev.off()

rsync -avzP ericacnr@colostate.edu@login.rc.colorado.edu:/scratch/alpine/ericacnr@colostate.edu/BANS/06.genomic_offset/RDAForest_offset/results/turnover_curves.png \
/Users/ericarobertson/Desktop/BANS_adaptive_units/analysis/06.genomic_offset/results

################## VARIABLE SELECTION ################## 
rownames(cordist)=geno$IID
mm=mtrySelJack(Y=cordist,X=env,covariates=latlon.gcd,nreps=31, prop.positive.cutoff=0.5,top.pcs=tokeep)

mm$goodvars
save(mm, "mm.Rdata")

png("results/good_variables.png")
ggplot(mm$prop.positive,aes(var,prop.positive))+
  geom_bar(stat="identity")+
  coord_flip()+
  geom_hline(yintercept=0.5,col="red")
dev.off()

rsync -avzP ericacnr@colostate.edu@login.rc.colorado.edu:/scratch/alpine/ericacnr@colostate.edu/BANS/06.genomic_offset/RDAForest_offset/results/good_variables.png \
/Users/ericarobertson/Desktop/BANS_adaptive_units/analysis/06.genomic_offset/results
################## FORMING PREDICTIONS ################## 

env <- as.data.frame(env)
# present-day
oj=ordinationJackknife(Y=cordist,X=env[, mm$goodvars, drop = FALSE], newX=envc[,mm$goodvars,drop = FALSE]
                       ,covariates=latlon.gcd,nreps=25,top.pcs=tokeep,extra=0.1)
# future
ojf=ordinationJackknife(Y=cordist,X=env[, mm$goodvars, drop = FALSE], newX=envf[,mm$goodvars,drop = FALSE]
                        ,covariates=latlon.gcd,nreps=25,top.pcs=tokeep,extra=0.2)

png("results/boxplots_importance.png")
ggplot(oj$all.importances,aes(variable,importance))+geom_boxplot(outlier.shape = NA)+coord_flip()
dev.off()

rsync -avzP ericacnr@colostate.edu@login.rc.colorado.edu:/scratch/alpine/ericacnr@colostate.edu/BANS/06.genomic_offset/RDAForest_offset/results/boxplots_importance.png \
/Users/ericarobertson/Desktop/BANS_adaptive_units/analysis/06.genomic_offset/results

################## MAPS OF PREDICTED ADAPTIVE NEIGHBORHOODS ################## 
goods=oj$goodrows
# predictor data restricted to only those spots:
ras2=envc[which(goods),]
xy2=envc[goods,c("x","y")]
names(xy2)=c("lon","lat")
rfpreds=oj$predictions.direct
turnovers=oj$predictions.turnover
bests=names(oj$median.importance)[1:3]

### CLUSTERING INTO ADAPTIVE NEIGHBORHOODS

pa1=plot_adaptation(rfpreds,ras2[,bests],xy2,main="direct preds",
                    # options affecting PCA plot:
                    rangeExp=1.5,
                    scal=10,
                    jitscale=0.05,
                    # options affecting map and PCA colors:
                    color.scheme="001",
                    lighten=0.8,
                    # options affecting clustering:
                    cluster.guide = NULL,
                    nclust=15,
                    cluster.merge.level=0.333 # the default
)

ggsave("results/adaptation_plot.png", plot = pa1, width = 8, height = 6, dpi = 300)

rsync -avzP ericacnr@colostate.edu@login.rc.colorado.edu:/scratch/alpine/ericacnr@colostate.edu/BANS/06.genomic_offset/RDAForest_offset/results/adaptation_plot.png \
/Users/ericarobertson/Desktop/BANS_adaptive_units/analysis/06.genomic_offset/results

pa2=plot_adaptation(rfpreds,ras2[,bests],xy2,main="turnovers",
                    # options affecting PCA plot:
                    rangeExp=1.5,
                    scal=10,
                    jitscale=0.05,
                    # options affecting map and PCA colors:
                    color.scheme="001",
                    lighten=0.8,
                    # options affecting clustering:
                    cluster.guide = turnovers,
                    nclust=12,
                    cluster.merge.level=0.25
)
ggsave("results/less_noisy_adaptation_plot.png", plot = pa1, width = 8, height = 6, dpi = 300)

# rsync -avzP ericacnr@colostate.edu@login.rc.colorado.edu:/scratch/alpine/ericacnr@colostate.edu/BANS/06.genomic_offset/RDAForest_offset/results/less_noisy_adaptation_plot.png \
# /Users/ericarobertson/Desktop/BANS_adaptive_units/analysis/06.genomic_offset/results


################## OFFSET ################## 
OFFS=gen_offset_oj(X=oj,Y=ojf,sx=envc[,1:2],sy=envf[,1:2])

landscape=st_as_sf(ne_countries(scale="medium",continent="north america"))
# fancy way
pdf("results/BANS_gen_offset.pdf",height=5, width=7)
gg=ggplot()+geom_sf(data=landscape)+
  xlim(min(OFFS$x),max(OFFS$x))+
  ylim(min(OFFS$y), max(OFFS$y))+
  theme_minimal()+
  geom_raster(data=OFFS,aes(x=x,y=y,fill=offset))+scale_fill_viridis(option="inferno",direction = -1)
plot(gg)
dev.off()






