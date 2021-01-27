library(phyloseq); packageVersion("phyloseq")
library(ggplot2)
library(gridExtra)
library(compositions)

#############
# load data #
#############

## load microbiome data
load('./design_AG/agdata.rda')

## load sample/matched_data
load('./design_AG/dat_matched_smoke_ag.RData')

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
ps_prune <- prune_taxa(empty_species != 0, ps)

## agglomerate to Genus ##
ps_Genus <- tax_glom(ps_prune, taxrank = "Genus", NArm = FALSE)
vec_taxa_Gen <- apply(otu_table(ps_Genus), 1, function(x) sum(x > 0, na.rm = TRUE))
table(vec_taxa_Gen)
# 5% prevalence filtering
ps_Genus_prune <- prune_taxa(vec_taxa_Gen >= 25, ps_Genus)
ps_Genus_prune

ps <- ps_Genus_prune

p = nrow(otu_table(ps))
n = ncol(otu_table(ps))
  
## prevalence
taxa_prev <- apply(unname(otu_table(ps)), 1, function(x) sum(x > 0, na.rm = TRUE)/n)
  
## variance
ps_ait <- transform_sample_counts(ps, function(x) {log((x + 1)/sum(x))})
taxa_var_ait <- apply(unname(otu_table(ps_ait)), 1, function(x) var(x))

#########################
### SELECT REFERENCES ###
#########################

plot(taxa_prev, taxa_var_ait)

condition <- taxa_var_ait < 3 & taxa_prev > .9
condition_ref <- which(condition)

# check if at least one reference taxa is in all samples
reference_values <- apply(otu_table(ps)[condition_ref,], 2, sum)
table(reference_values)[1:3]

plot(taxa_prev, taxa_var_ait, col = factor(condition))

# save(condition_ref, file = "./4_differential_abundance_AG/selected_ref_genus_AG.RData")

