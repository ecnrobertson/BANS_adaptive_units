# #!/bin/bash
# #SBATCH --job-name=RDA
# #SBATCH --output=RDA.%j.out
# #SBATCH --error=RDA.%j.err
# #SBATCH -t 8:00:00
# #SBATCH -p amilan
# #SBATCH --qos=normalm
# #SBATCH --nodes=1
# #SBATCH --ntasks-per-node 24
# #SBATCH --mem=90G
# #SBATCH --mail-type=ALL
# #SBATCH  --mail-user=ericacnr@colostate.edu
# 
# source ~/.bashrc
# 
# ############ Bash to make the .raw file ############
# bfile="/scratch/alpine/ericacnr@colostate.edu/BANS/02.5.neutral_pop_struc/pruned_allindv"
# conda activate GWAS2
# plink --bfile $bfile -aec --recode A --out data/BANS.all.imputed.pruned
# cut -d$'\t' -f7- data/BANS.all.imputed.pruned.raw > data/BANS.all.imputed.pruned_noheader.raw

# ############ Load and Run RDA.R script ############
# conda activate R
# # Run your R script
# Rscript 03.2.RDA.R

############################# RDA.R ##################################
##### File Checklist #####
# rsync -avzP /Users/ericarobertson/Desktop/BANS_adaptive_units/analysis/03.GEA/bg_colors.rds ericacnr@colostate.edu@login.rc.colorado.edu:/scratch/alpine/ericacnr@colostate.edu/BANS/03.GEA/data/
# rsync -avzP /Users/ericarobertson/Desktop/BANS_adaptive_units/data/BANS_all_sample_data.csv ericacnr@colostate.edu@login.rc.colorado.edu:/scratch/alpine/ericacnr@colostate.edu/BANS/03.GEA/data/
# rsync -avzP /Users/ericarobertson/Desktop/BANS_adaptive_units/analysis/02.5.delineate_ESUs/bg_colors.rds  ericacnr@colostate.edu@login.rc.colorado.edu:/scratch/alpine/ericacnr@colostate.edu/BANS/03.GEA/data/
# In data/
# BANS_sample_list_pop.tsv
# BANS.all.imputed.pruned_noheader.raw
# pruned_allindv.fam
# BANS_topRDA_vars.csv
# BANS_all_sample_data.csv
# bg_colors.rds

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

############ PCA for pop structure variables ############ 
plink_pca <- read.table("../02.5.neutral_pop_struc/PC1_allindv.ld25-10-0.5.eigenvec", header = F) %>% select(-V1)
n_pcs <- ncol(plink_pca) - 2  # number of PC columns
colnames(plink_pca) <- c("FID", "IID", paste0("PC", 1:n_pcs))
plink_pca <- plink_pca %>% dplyr::select(-FID) %>% rename(BGP_ID = IID)
print("read in PCA variables")

############ MEMS for spatial autocorrelation ############ 
indorder <- fread("data/pruned_allindv.fam") %>% rename(BGP_ID = V1) %>% dplyr::select(BGP_ID)
coords <- read_csv("data/BANS_all_sample_data.csv") %>%
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
selmem <- vegan::ordiR2step(mem_rda_full, mem_rda, R2permutations=1000,Pin=0.01,R2scope=T)
summary(selmem)
# extract the names of the selected mem variables
selmem <- names(selmem$terminfo$ordered)

# subset the mem matrix to include only the mems that explain population structure
mem_popstr <- mem[,selmem]

saveRDS(mem_popstr, file = "data/BANS.mems.popstr_updated.rds")

coords_mems_pc <- cbind(coords, mem, plink_pca[,c("PC1", "PC2")])

print("MEMS identified correctly")
############ RDA models prep ############ 

plink_data <- read.PLINK("data/BANS.all.imputed.pruned_noheader.raw")
print(nrow(plink_data))

# Convert genlight to genotype matrix
geno_matrix <- as.matrix(plink_data)

# Keep SNP names
snp_names <- locNames(plink_data)  

# Recreate genlight object with SNP names
plink_subset <- new("genlight", geno_matrix, loc.names = snp_names, ind.names = indNames(plink_data))
print(nrow(plink_subset))

# Move to Dataframe format
geno_df <- as.data.frame(plink_subset)
rownames(geno_df) <- indNames(plink_data)
print(nrow(geno_df)==nrow(plink_subset))

pops <- read.table(file = "data/BANS_sample_list_pop.tsv", sep="\t", header = FALSE) %>% rename(BGP_ID = V1, Pop = V3) %>% dplyr::select(BGP_ID, Pop)
indorder <- fread("data/pruned_allindv.fam") %>% rename(BGP_ID = V1) %>% dplyr::select(BGP_ID)

env <- fread("data/BANS_topRDA_vars.csv")
env <- left_join(indorder, env, by = "BGP_ID")

print(env$BGP_ID==rownames(geno_df))
bgp_ids <- data.frame(BGP_ID = rownames(geno_df))

env <- cbind(pop, env)

# Convert BGP_ID to character (to avoid factor-related issues)
env[, BGP_ID := as.character(BGP_ID)]

# Convert ClimGroup to a factor if it's categorical
env[, Pop := as.factor(Pop)]

pred <- env[,8:17]

pred_mems <- cbind(pred, as.data.frame(mem_popstr))

pred_pc_mems_env <- cbind(pred_mems, plink_pca %>% dplyr::select(PC1, PC2))


print("RDA data loaded")
############ Actual RDA Models ############ 

### 1. Fit RDAs ----
print("fit RDAs")
# Baseline (env only)
rda_env <- rda(geno_df ~ Temperate.or.Subpolar.Shrubland + Subpolar.Taiga.Needleleaf.Forest + bio05 + Temperate.or.Subpolar.Grassland + Wetland + bio15 + clay + Temperate.or.Subpolar.Broadleaf.Deciduous.Forest + sand + bio02, data = pred_pc_mems_env, scale = TRUE)
saveRDS(rda_env, file = "results/RDA_output/rda_env.rds")
print("rda_env")
# Control for MEMs (full, uncorrelated, single MEM1)
rda_env_mems_uncor <- rda(geno_df ~ Temperate.or.Subpolar.Shrubland + Subpolar.Taiga.Needleleaf.Forest + bio05 + Temperate.or.Subpolar.Grassland + Wetland + bio15 + clay + Temperate.or.Subpolar.Broadleaf.Deciduous.Forest + sand + bio02 +
                            Condition(MEM2 + MEM4 + MEM8 + MEM6),
                          data = pred_pc_mems_env, scale = TRUE)
saveRDS(rda_env_mems_uncor, file = "results/RDA_output/rda_env_mems_uncor.rds")
print("rda_env_mems_uncor")
rda_env_mem2 <- rda(geno_df ~ Temperate.or.Subpolar.Shrubland + Subpolar.Taiga.Needleleaf.Forest + bio05 + Temperate.or.Subpolar.Grassland + Wetland + bio15 + clay + Temperate.or.Subpolar.Broadleaf.Deciduous.Forest + sand + bio02 +
                      Condition(MEM2),
                    data = pred_pc_mems_env, scale = TRUE)

# Control for PCs
rda_env_pc1 <- rda(geno_df ~ Temperate.or.Subpolar.Shrubland + Subpolar.Taiga.Needleleaf.Forest + bio05 + Temperate.or.Subpolar.Grassland + Wetland + bio15 + clay + Temperate.or.Subpolar.Broadleaf.Deciduous.Forest + sand + bio02 +
                     Condition(PC1),
                   data = pred_pc_mems_env, scale = TRUE)

rda_env_pc12 <- rda(geno_df ~ Temperate.or.Subpolar.Shrubland + Subpolar.Taiga.Needleleaf.Forest + bio05 + Temperate.or.Subpolar.Grassland + Wetland + bio15 + clay + Temperate.or.Subpolar.Broadleaf.Deciduous.Forest + sand + bio02 +
                      Condition(PC1 + PC2),
                    data = pred_pc_mems_env, scale = TRUE)

############ RDA Model Comparison ############ 
### Helper function to extract key stats from an RDA ----
rda_summary <- function(model, name) {
  r2   <- RsquareAdj(model)$adj.r.squared
  pval <- anova.cca(model, parallel = getOption("mc.cores"))$`Pr(>F)`[1]  # model-level test
  eigvals <- eigenvals(model, model = "constrained")
  axis1_var <- eigvals[1] / sum(eigvals) # variance on axis 1
  data.frame(Model = name, AdjR2 = r2, Pval = pval, Axis1_var = axis1_var)
}

### Build comparison table ----
results <- rbind(
  rda_summary(rda_env,        "Env only"),
  rda_summary(rda_env_mems_uncor, "Env + MEMs (uncor)"),
  rda_summary(rda_env_mem1,   "Env + MEM1"),
  rda_summary(rda_env_pc1,    "Env + PC1"),
  rda_summary(rda_env_pc12,   "Env + PC1+PC2")
)

write.csv(results,
          file = "RDA_model_summaries.csv",
          row.names = FALSE)
print("model comparison table done")

############ Plotting RDA Results ############
clim_groups <- unique(env$Pop)
eco <- factor(env$Pop, levels = clim_groups)  # Ensure factor levels match
bg <- readRDS("data/bg_colors.rds")

rda_obj <- readRDS("results/RDA_output/rda_env_mems_uncor.rds")
env_df  <- readRDS("results/RDA_output/env.rds")

region_order <- c(
  "Cluster_1", "Cluster_5", "Cluster_16", "Cluster_4", "Cluster_9", "Cluster_12", "Cluster_2", "Cluster_20",
  "Cluster_3", "Cluster_8", "Cluster_17", "Cluster_21", "Cluster_18", "Cluster_19",
  "Cluster_6", "Cluster_14", "Cluster_11", "Cluster_10", "Cluster_13", "Cluster_15", "Cluster_7"
)

bg <- readRDS("../02.5.delineate_ESUs/bg_colors.rds")
bg <- bg[region_order]  # enforce consistent order

env_df  <- env_df[,-1:-3]                  # must contain your metadata (and Long/Lat if mapping)
group_col <- "Pop"              # change to "ClimGroup" if that's what you want
scaling_val <- 3

# 1) extract scores
site_scores <- scores(rda_obj, display = "sites", scaling = scaling_val)
bp_scores   <- scores(rda_obj, display = "bp",    scaling = scaling_val)

sites_df <- as.data.frame(site_scores) %>%
  tibble::rownames_to_column(var="IndID")

bp_df <- as.data.frame(bp_scores) %>%
  tibble::rownames_to_column("Var")

# 2) Join group labels safely by ID
env2 <- env_df %>%
  select(BGP_ID, all_of(group_col)) %>%
  rename(PopID = all_of(group_col), IndID=BGP_ID)   # always call it PopID in the plotting df

sites_df <- left_join(sites_df, env2, by = "IndID")

# 3) PDF biplot (colored by ClimGroup)
eco <- factor(env_df$Pop, levels = region_order)

eig <- summary(rda_obj)$cont$importance
prop <- eig["Proportion Explained", ]
pct <- round(100 * prop, 1)
xlab <- paste0("RDA1 (", pct[1], "%)")
ylab <- paste0("RDA2 (", pct[2], "%)")

pdf("bans.allind.allsnp.rda_env_mems_uncor.pdf", width = 10, height = 8)
par(mar = c(5, 5, 4, 8))
plot(rda_obj,
  type = "n",
  scaling = 3,
  xlab = xlab,
  ylab = ylab)
# sites
points(rda_obj,
  display = "sites",
  pch = 21,
  cex = 1.2,
  col = "gray32",
  bg  = bg[eco],
  scaling = 3)
# bp vectors: draw arrows + labels (clearer than text-only)
bp <- scores(rda_obj, display = "bp", scaling = 3)
# compute a sensible multiplier
mul <- ordiArrowMul(bp)
arrows(0, 0,
  bp[, 1] * mul,
  bp[, 2] * mul,
  length = 0.08,
  col = "black")

text(bp[, 1] * mul,
  bp[, 2] * mul,
  labels = rownames(bp),
  col = "black",
  cex = 0.9,
  pos = 3)
dev.off()

print("plot of RDA results with best model")

############ Getting Outlier SNPs ############
load.rda <- scores(rda_env_mems_uncor, choices=c(1:3), display="species")

outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)   #find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]             #locus names in these tails
}

cand1 <- outliers(load.rda[,1],3)
cand2 <- outliers(load.rda[,2],3)
cand3 <- outliers(load.rda[,3],3)

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

##Add environmental correlations to candidate snps
foo <- matrix(nrow=(ncand), ncol=ncol(pred))  #4 columns for 4 predictors
colnames(foo) <- colnames(pred)

for (i in 1:length(df.cand$snp)) {
  nam <- df.cand[i,2]
  snp.gen <- geno_df[[nam]]
  foo[i,] <- apply(pred,2,function(x) cor(x,snp.gen))
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



