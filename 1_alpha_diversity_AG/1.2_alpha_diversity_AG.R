library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(DivNet)
library(doParallel)
library(doSNOW)
library(gridExtra)

#############
# load data #
#############

## load microbiome data
load('./design_AG/agdata.rda')

## load sample/matched_data
load('./design_AG/dat_matched_smoke_ag.RData')

## load W matrix for randomization test
load("./design_AG/W_paired_smoke_ag.Rdata")

############################
# create a phyloseq object # 
############################

## order covariates data frame
sample_df <- matched_df[order(rownames(matched_df)),]
sample_df$W <- as.factor(sample_df$W)

# create a phyloseq object
ps <- phyloseq(otu_table(agdata),
               sample_data(sample_df),
               tax_table(agdata))
ps

## locate the species that are totally absent in the matched data
empty_species <- rowSums(otu_table(ps))
length(which(empty_species == 0))
# remove them
ps_prune <- prune_taxa(empty_species > 1, ps)

#############################################################
### 1. ESTIMATE THE TOTAL ALPHA DIVERSITY FOR EACH SAMPLE ###
#############################################################

## agglomerate data to phylum level
ps_phyl <- tax_glom(ps_prune, taxrank="Phylum", NArm = FALSE)

# calculate plug-in shannon index
shannon_plug <- sample_shannon(ps_phyl)

# pick which taxon is to be the base ()
base_abundant_taxa <- rownames(tax_table(ps_phyl)[which(tax_table(ps_phyl)[,"Phylum"] == "Firmicutes")])

## estimate shannon with divnet function  
divnet_phylum <- divnet(ps_phyl, 
                        base = base_abundant_taxa,
                        ncores = 4)
divnet_phylum

## retrieve statistics for inference
sample_data(ps_prune)[,"DivNet_W"] <- summary(divnet_phylum$shannon)$estimate
sample_data(ps_prune)[,"DN_error_W"] <- summary(divnet_phylum$shannon)$error
sample_data(ps_prune)[,"lower"] <- summary(divnet_phylum$shannon)$lower
sample_data(ps_prune)[,"upper"] <- summary(divnet_phylum$shannon)$upper
sample_data(ps_prune)[,"sample_counts"] <- sample_sums(ps_prune)
sample_data(ps_prune)[,"shannon_plug"] <- summary(shannon_plug)$estimate
sample_data(ps_prune)[,"id"] <- rownames(sample_data(ps_prune))

# plot raw data
ggplot(sample_data(ps_prune), aes(x = id)) + 
  geom_point(aes(y = DivNet_W), colour = "red") + 
  geom_errorbar(aes(ymin=lower, ymax=upper)) +
  geom_point(aes(y = shannon_plug), colour = "blue", alpha = .4) +
  xlab('Samples') + ylab("Shannon estimate")

# code run on cluster (computationally intensive)
# load('/Users/alicesommer/Desktop/DivNet_cluster/ps_DivNet_smoke.RData')

g_PM <- ggplot(sample_data(ps_prune), aes(color = factor(W), y = DivNet_W)) +
  geom_boxplot(alpha = .5) + ylab('DivNet shannon index') +
  scale_x_discrete(name = "") +
  scale_colour_manual(values = c("darkgreen","blue4"), limits=c("1","0"), name ="Long-term PM2.5", labels = c("Low","High")) +
  # scale_colour_manual(values = c("darkgreen","blue4"), limits=c("1","0"), name ="Smoking", labels = c("No","Yes")) +
  theme(legend.position = "top", legend.key.size =  unit(0.1, "in")) +
  guides(color=guide_legend(nrow=2,byrow=TRUE))  

#######################################################
### 2. RANDOMIZATION-TEST WITH BETTA TEST STATISTIC ###
#######################################################

# import modified function for betta funtions
source("./misc/model_betta.R")

x <- cbind(1, 
           # sample_data(ps_prune)$u3csex
           sample_data(ps_prune)$W )

reg <- betta(sample_data(ps_prune)$DivNet_W,
             sample_data(ps_prune)$DN_error_W , X = x)
reg$table
estim_obs <- reg$table[2,1]

## W matrix 
W_paired_smoke <- W_paired_smoke[order(rownames(W_paired_smoke)), ]
dim(W_paired_smoke)

# set the number of randomizations
nrep <- ncol(W_paired_smoke)/2000

# create a matrix where the t_rand will be saved
t_array <- NULL

for(i in 1:nrep){
  print(i)
  x = cbind(1, 
            # sample_data(ps_prune)$u3csex
            W_paired_smoke[,i] )
  
  reg = betta(sample_data(ps_prune)$DivNet_W,
              sample_data(ps_prune)$DN_error_W , X = x)
  
  # fill t_array
  t_array[i] = reg$table[2,1] 
}

## calculate p_value
# p_value <- mean(t_array >= t_array[1])
p_value <- mean(t_array >= estim_obs)
p_value

## plot distribution of test-statistic
hist(t_array, breaks = 30, main = "", xlab = "shannon beta (Smoking)")
abline(v = estim_obs, col = 'red', lwd = 2, lty = 2)



