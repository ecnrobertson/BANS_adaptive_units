# #!/bin/bash
# #SBATCH --job-name=RDA
# #SBATCH --output=RDA%j.out
# #SBATCH --error=RDA%j.err
# #SBATCH -t 8:00:00
# #SBATCH -p amilan
# #SBATCH --nodes=1
# #SBATCH --ntasks-per-node 24
# #SBATCH --mem=90G
# #SBATCH --mail-type=ALL
# #SBATCH  --mail-user=ericacnr@colostate.edu
# 
# source ~/.bashrc
# 
# ############ Bash to make the .raw file ############
# conda activate GWAS2
# plink --bfile BANS.ds6x.maf.0.05.SNP.above4x.maxmiss.8.imputed4.1.ld25-10-0.5 -aec --recode A --out data/BANS.ds6x.maf.0.05.SNP.above4x.maxmiss.8.imputed4.1.ld25-10-0.5
# cut -d$'\t' -f7- data/BANS.ds6x.maf.0.05.SNP.above4x.maxmiss.8.imputed4.1.ld25-10-0.5.raw > data/BANS.ds6x.maf.0.05.SNP.above4x.maxmiss.8.imputed4.1.ld25-10-0.5_noheader.raw
# 
# ############ Load and Run RDA.R script ############
# conda activate R
# # Run your R script
# Rscript 03.2.RDA.R

############ Libraries ############ 
library(dplyr)       # select, rename, mutate, pipes (%>%)
library(readr)       # read_delim, read_csv
library(data.table)  # fread
library(sf)              # st_as_sf
library(sp)              # as(..., "Spatial")
library(geosphere)       # distm, distHaversine
library(adespatial)      # dbmem
library(rnaturalearth)   # ne_states
library(adegenet)   # read.PLINK, genlight, locNames, indNames
library(vegan)      # rda, ordiR2step, RsquareAdj, scores, anova.cca, eigenvals
library(ggplot2)    # ggplot
library(RColorBrewer) # brewer.pal

############ reading in the data we need ############ 
pops <- read.table(file = "data/BANS_sample_list_pop.tsv", sep="\t", header = FALSE) %>% rename(BGP_ID = V1, Pop = V3) %>% dplyr::select(BGP_ID, Pop)

genotypes <- fread("data/subsampled_noheader.raw")

indorder <- fread("BANS.ds6x.maf.0.05.SNP.above4x.maxmiss.8.imputed4.1.ld25-10-0.5.fam") %>% rename(BGP_ID = V1) %>% dplyr::select(BGP_ID)

env <- read_delim("data/BANS.initial.worldclim.landscape_edit.txt", delim = "\t") %>% select(-cfvo)

env <- left_join(indorder, env, by = "BGP_ID")

env$BGP_ID == indorder$BGP_ID

#BIO16 = Precipitation of Wettest Quarter
#bioclay = g/kg
#BIO8 = Mean Temperature of Wettest Quarter
#BIO10 = Mean Temperature of Warmest Quarter
#BIO11 = Mean Temperature of Coldest Quarter

env %>% as.data.frame() %>% dplyr::select(-"BGP_ID", bio16, clay, bio08, bio10, bio11) %>% write.table("data/BANS.climtop5.txt", row.names = F, quote = F, sep = "\t")

#env2 <- fread("data/BANS.climtop4.txt", sep = "\t")
print("data read in correctly")
############ PCA for pop structure variables ############ 
plink_pca <- read.table("BANS.ds6x.maf.0.05.SNP.above4x.maxmiss.8.imputed4.1.ld25-10-0.5.eigenvec")
n_pcs <- ncol(plink_pca) - 2  # number of PC columns
colnames(plink_pca) <- c("FID", "IID", paste0("PC", 1:n_pcs))
plink_pca <- plink_pca %>% dplyr::select(-FID) %>% rename(BGP_ID = IID)
print("read in PCA variables")
############ MEMS for spatial autocorrelation ############ 
coords <- read_csv("data/BANS_sample_data.csv") %>%
  rename(BGP_ID = "BGP ID") %>% 
  filter(BGP_ID %in% indorder$BGP_ID) %>% 
  dplyr::select(BGP_ID, Long, Lat)

# reorder coordinates to match PCA data
print("reorder coordinates to match PCA data")
coords <- coords[(match(indorder$BGP_ID, coords$BGP_ID)),]

coords_sf <- st_as_sf(coords, coords = c("Long", "Lat"), crs = 4326)

gd <- geosphere::distm(as(coords_sf, "Spatial"), fun = geosphere::distHaversine)
dist <- as.dist(gd)

# generate moran's eigenvector maps (MEMs) based on coordinates
print("generate MEMS based on coordinates")
mem <- adespatial::dbmem(dist, MEM.autocor = "all")

#moran_test <- moran.randtest(mem)

# rda to assess how much MEMs explain genetic variation (PC1–PC2)
mem_rda <- vegan::rda(plink_pca[,c("PC1", "PC2")] ~ ., as.data.frame(mem[,1:10]))
summary(mem_rda)

# create a null (intercept-only) rda model as a baseline
mem_rda_full <- vegan::rda(plink_pca[,c("PC1", "PC2")] ~ 1, as.data.frame(mem[,1:10]))
summary(mem_rda_full)

# stepwise forward selection of mem variables based on adjusted R-square and significance
selmem <- vegan::ordiR2step(mem_rda_full, mem_rda, Pin = .01)
summary(selmem)
# extract the names of the selected mem variables
selmem <- names(selmem$terminfo$ordered)

# subset the mem matrix to include only the mems that explain population structure
mem_popstr <- mem[,selmem]

saveRDS(mem_popstr, file = "BANS.mems.popstr_updated.rds")

coords_mems_pc <- cbind(coords, mem, plink_pca[,c("PC1", "PC2")])

map <- rnaturalearth::ne_states(country = c("United States of America", "Canada", "Mexico"), returnclass = 'sf')

# plot the populations you defined. use plotly to zoom around easy
# mem4_map<- ggplot(map) + geom_sf() + geom_point(data = coords_mems_pc, aes(x = Long, y = Lat, color = MEM4, BGP_ID = BGP_ID)) +
#   coord_sf(xlim = c(-125, -70), ylim = c(23, 58)) + scale_color_viridis_c()
# ggsave("mem4_map.pdf", mem4_map, width = 8, height = 6)
# 
# 
# pdf("coords_mems_pc_plot.pdf", width = 7, height = 7)
# plot(coords_mems_pc[, c("PC1", "PC2", paste0("MEM", 1:9))])
# dev.off()
print("MEMS identified correctly")
############ RDA models prep ############ 

plink_data <- read.PLINK("data/BANS.ds6x.maf.0.05.SNP.above4x.maxmiss.8.imputed4.1.ld25-10-0.5_noheader.raw")
print(nrow(plink_data))

# Convert genlight to genotype matrix
geno_matrix <- as.matrix(plink_data)

# Keep SNP names
snp_names <- locNames(plink_data)  

# Recreate genlight object with SNP names
plink_subset <- new("genlight", geno_matrix, loc.names = snp_names, ind.names = indNames(plink_data))
print(nrow(plink_subset))

# Save as CSV
geno_df <- as.data.frame(plink_subset)
nrow(geno_df)

env <- fread("data/BANS.climtop5.txt")
pop <- fread("data/BANS_sample_list_pop.tsv", header = F)  %>% dplyr::select(-V2) %>% rename(BGP_ID = V1, ClimGroup = V3)

bgp_ids <- data.frame(BGP_ID = rownames(geno_df))
pop <- pop %>% slice(match(bgp_ids$BGP_ID, pop$BGP_ID)) ## make ind same order as snpfile & env data
env <- cbind(pop, env)

##modify data types like brenna did...

# Convert BGP_ID to character (to avoid factor-related issues)
env[, BGP_ID := as.character(BGP_ID)]

# Convert ClimGroup to a factor if it's categorical
env[, ClimGroup := as.factor(ClimGroup)]

all(rownames(geno_df) == env$BGP_ID)

pred <- subset(env, select=c("bio02", "bio05", "bio07", "bio03"))

pred_mems <- cbind(pred, as.data.frame(mem_popstr))

pred_pc <- cbind(pred, plink_pca %>% dplyr::select(PC1, PC2))

#this is the null model?
bans.full_geno_df.ind.rda <- rda(geno_df ~ ., data=pred, scale=T)

print("RDA data loaded and null model fit")
############ Actual RDA Models ############ 

### 1. Fit RDAs ----
print("fit RDAs")
# Baseline (env only)
rda_env      <- rda(geno_df ~ bio02 + bio05 + bio07 + bio03, data = pred, scale = TRUE)

# Control for MEMs (full, uncorrelated, single MEM1)
rda_env_mems <- rda(geno_df ~ bio02 + bio05 + bio07 + bio03 +
                      Condition(MEM1 + MEM5 + MEM4 + MEM3 + MEM7),
                    data = pred_mems, scale = TRUE)

rda_env_mems_uncor <- rda(geno_df ~ bio02 + bio05 + bio07 + bio03 +
                            Condition(MEM4 + MEM3 + MEM7),
                          data = pred_mems, scale = TRUE)

rda_env_mem1 <- rda(geno_df ~ bio02 + bio05 + bio07 + bio03 +
                      Condition(MEM1),
                    data = pred_mems, scale = TRUE)

# Control for PCs
rda_env_pc1   <- rda(geno_df ~ bio02 + bio05 + bio07 + bio03 +
                       Condition(PC1),
                     data = pred_pc, scale = TRUE)

rda_env_pc12  <- rda(geno_df ~ bio02 + bio05 + bio07 + bio03 +
                       Condition(PC1 + PC2),
                     data = pred_pc, scale = TRUE)


### 2. Plots ----
print("plot RDA results")
# Save all plots into one PDF
pdf("RDA_models_updated.pdf", width = 10, height = 12)

# Arrange 2 plots per page (2 rows, 1 column)
par(mfrow = c(2, 1), mar = c(5, 5, 4, 2))

### Baseline
plot(rda_env, scaling = 3, main = "Baseline RDA (Axis 1 vs 2)")
plot(rda_env, choices = c(1, 3), scaling = 3, main = "Baseline RDA (Axis 1 vs 3)")

### MEMs
plot(rda_env_mems, scaling = 3, main = "RDA with MEMs (Axis 1 vs 2)")
plot(rda_env_mems, choices = c(1, 3), scaling = 3, main = "RDA with MEMs (Axis 1 vs 3)")

plot(rda_env_mems_uncor, scaling = 3, main = "RDA with Uncorrelated MEMs (Axis 1 vs 2)")
plot(rda_env_mem1, scaling = 3, main = "RDA with MEM1 only")

### PCs
plot(rda_env_pc1, scaling = 3, main = "RDA with PC1 (Axis 1 vs 2)")
plot(rda_env_pc1, choices = c(1, 3), scaling = 3, main = "RDA with PC1 (Axis 1 vs 3)")
plot(rda_env_pc12, scaling = 3, main = "RDA with PC1 + PC2 (Axis 1 vs 2)")

dev.off()

############ RDA Model Comparison ############ 
### Helper function to extract key stats from an RDA ----
rda_summary <- function(model, name) {
  r2   <- RsquareAdj(model)$adj.r.squared
  eig  <- summary(eigenvals(model, model = "constrained"))[1]             # variance on axis 1
  data.frame(Model = name, AdjR2 = r2, Axis1_var = eig)
}

### Build comparison table ----
results <- rbind(
  rda_summary(rda_env,        "Env only"),
  rda_summary(rda_env_mems,   "Env + MEMs"),
  rda_summary(rda_env_mems_uncor, "Env + MEMs (uncor)"),
  rda_summary(rda_env_mem1,   "Env + MEM1"),
  rda_summary(rda_env_pc1,    "Env + PC1"),
  rda_summary(rda_env_pc12,   "Env + PC1+PC2")
)

write.csv(results,
          file = "RDA_model_summaries_update.csv",
          row.names = FALSE)
print("model comparison table done")
############ Plotting RDA Results ############
region_order <- c(
  "Cluster_1", "Cluster_5", 
  "Cluster_4", "Cluster_8", "Cluster_12", "Cluster_2", "Cluster_3", "Cluster_7",
  "Cluster_11", "Cluster_14", "Cluster_13", "Cluster_6", "Cluster_9", "Cluster_10"
)

eco <- factor(env$ClimGroup, levels = region_order)

bg <- colorRampPalette(brewer.pal(min(12, length(unique(eco))), "Paired"))(length(unique(eco)))

# Define your clusters grouped into regions
West_clusters    <- c("Cluster_1", "Cluster_5")
Midwest_clusters <- c("Cluster_4", "Cluster_8", "Cluster_12", "Cluster_2", "Cluster_3", "Cluster_7")
East_clusters    <- c("Cluster_11", "Cluster_14", "Cluster_13", "Cluster_6", "Cluster_9", "Cluster_10")

# Assign colors by region
W.cols <- setNames(brewer.pal(n = length(West_clusters), "Reds"), West_clusters)
Midwest.cols <- setNames(brewer.pal(n = length(Midwest_clusters), "Blues"), Midwest_clusters)
E.cols <- setNames(brewer.pal(n = length(East_clusters), "Greens"), East_clusters)

# Combine into one named vector
bg <- c(W.cols, Midwest.cols, E.cols)

# Reorder bg to match your factor order
bg <- bg[region_order]

# Check with barplot
barplot(rep(1, length(bg)), col = bg, border = NA, names.arg = names(bg), las = 2)

par(mar=c(5, 5, 4, 8))
plot(rda_env_mems_uncor, type="n", scaling=3)  # Set up empty plot
#points(rda_env_mems_uncor, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg[eco]) # Individuals
points(rda_env_mems_uncor, display="sites", pch=21, cex=1.3,
       col="gray32", scaling=3, bg=bg[as.character(eco)])
text(rda_env_mems_uncor, scaling=3, display="bp", col="black", cex=1)  # Predictors
#legend("bottomright", legend=levels(eco), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)
#legend("bottomright", legend=levels(eco), bty="n", col="gray32", pch=21, cex=0.6, pt.bg=bg)
dev.off()

print("plot of RDA results with best model")
############ Getting Outlier SNPs ############
load.rda <- scores(rda_env_mems_uncor, choices=c(1:3), display="species")

hist(load.rda[,1], main="Loadings on RDA1")
hist(load.rda[,2], main="Loadings on RDA2")
hist(load.rda[,3], main="Loadings on RDA3")

outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)   #find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]             #locus names in these tails
}

cand1 <- outliers(load.rda[,1],3) # 75998
cand2 <- outliers(load.rda[,2],3) # 63347
cand3 <- outliers(load.rda[,3],3) # 63725

##Find total # candidates
ncand <- length(cand1) + length(cand2) + length(cand3)
print("number of candidate snps")
print(ncand)

df.cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
df.cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
df.cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))

colnames(df.cand1) <- colnames(df.cand2) <- colnames(df.cand3) <- c("axis","snp","loading")

df.cand <- rbind(df.cand1, df.cand2, df.cand3)
df.cand$snp <- as.character(df.cand$snp)

df.cand$snp <- gsub("-", ".", df.cand$snp)

# pred <- subset(env, select=c(bio02 + bio05 + bio07 + bio03))

##Add environmental correlations to candidate snps
foo <- matrix(nrow=(ncand), ncol=4)  #4 columns for 4 predictors
colnames(foo) <- c("bio02", "bio05", "bio07", "bio03")
pred2 <- pred[,1:4]

for (i in 1:length(df.cand$snp)) {
  nam <- df.cand[i,2]
  snp.gen <- geno_df[[nam]]
  foo[i,] <- apply(pred2,2,function(x) cor(x,snp.gen))
}

cand <- cbind.data.frame(df.cand,foo)

length(cand$snp[duplicated(cand$snp)])
foo <- cbind(cand$axis, duplicated(cand$snp)) 
table(foo[foo[,1]==1,2]) # no duplicates on axis 1
table(foo[foo[,1]==2,2]) # 18202 duplicates on axis 2
table(foo[foo[,1]==3,2]) # 18787 duplicates on axis 3
cand <- cand[!duplicated(cand$snp),] # remove duplicate detection

cols <- as.numeric(ncol(cand))
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,(cols+1)] <- names(which.max(abs(bar[4:cols]))) # gives the variable
  cand[i,(cols+2)] <- max(abs(bar[4:cols]))              # gives the correlation
}

colnames(cand)[cols+1] <- "predictor"
colnames(cand)[cols+2] <- "correlation"

table(cand$predictor)

#rda_env_mems_uncor
write.table(cand,"ind.rda.env_memp.uncor.cand.snps_update.txt",row.names=F,sep = "\t", quote=F)
print("cand snp results table output")

sel <- cand$snp
env <- cand$predictor
env[env=="bio02"] <- '#1b9e77'
env[env=="bio05"] <- '#d95f02'
env[env=="bio07"] <- "purple4"
env[env=="bio03"] <-  '#e7298a'

col.pred <- rownames(rda_env_mems_uncor$CCA$v) # pull the SNP names

for (i in 1:length(sel)) {           # color code candidate SNPs
  foo <- match(sel[i],col.pred)
  col.pred[foo] <- env[i]
}

col.pred[grep("scaffold",col.pred)] <- '#f1eef6' # non-candidate SNPs
col.pred[!grepl("^#", col.pred)] <- "#f1eef6"

col.pred[!col.pred %in% env] <- "gray80"
empty <- col.pred
empty.outline <- ifelse(empty %in% c("gray80"), "gray50", "black")

empty.outline <- ifelse(empty == "#00FF0000", "#00FF0000", "gray32")

bg <- c('#1b9e77','#d95f02',"purple4",'#e7298a')

print("plotting adaptive snps")
plot(rda_env_mems_uncor, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1))
#points(cawa2.rda, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3)
points(rda_env_mems_uncor, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3)
text(rda_env_mems_uncor, scaling=3, display="bp", col="#0868ac", cex=1)
legend("bottomright", legend=c("bio02","bio05","bio07","bio03"), 
       bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)



