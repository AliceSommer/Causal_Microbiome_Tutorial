library(dacomp)
library(phyloseq); packageVersion("phyloseq")

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
ps_prune <- prune_taxa(empty_species != 0, ps)

## agglomerate to Genus ##
ps_Genus <- tax_glom(ps_prune, taxrank = "Genus", NArm = FALSE)
vec_taxa_Gen <- apply(otu_table(ps_Genus), 1, function(x) sum(x > 0, na.rm = TRUE))
table(vec_taxa_Gen)
# 5% prevalence filtering
ps_Genus_prune <- prune_taxa(vec_taxa_Gen >= 25, ps_Genus)
ps_Genus_prune

ps <- ps_Genus_prune

# reference sets selected
load('./4_differential_abundance_AG/selected_ref_genus_AG.RData')

nrep = 100

verbose = T

Y = sample_data(ps)$W # research variable

# Permutated "W" matrix
W_paired <- W_paired_smoke[order(rownames(W_paired_smoke)), ]
Y_perm <- W_paired_smoke[,1:nrep] # the first column of the Y_perm matrix is W_obs
  
# Microbiome Data
X = t(otu_table(ps))
X = as(X, "matrix")
  
# reference set
ind_reference_taxa = condition_ref
  
########################################
#### TEST WITH DACOMP-RATIO METHOD ####
########################################
  
p_test = ncol(X)
n_test = nrow(X)
stats_matrix_test = matrix(NA, ncol = p_test, nrow = nrep)
ratio_matrix_test = matrix(NA, nrow = n_test, ncol = 1)
  
# Compute reference values
reference_values = apply(X[,ind_reference_taxa], 1, sum)
  
print(paste('ref. sanity check:', sum(reference_values == 0)))
  
#iterate over taxa and test
for(t in 1:p_test){
    
  if(verbose)
    if(t%% ceiling(p_test/10) == 1)
      cat(paste0('Testing taxon : ',t,'/',p_test,' \n\r'))
    
  nom_test = X[,t]
  dnom_test = reference_values
    
  #no need to test reference taxa
  if(t %in% ind_reference_taxa){
    print(t)
    # nom_test = apply(X_test[,ind_reference_taxa[-which(ind_reference_taxa == t)]], 1, sum)
    next
  }
    
  ratio_matrix_test[,1] = nom_test/(dnom_test+nom_test) ##### references not present in all samples
    
  stats_matrix_test[,t] = dacomp:::Compute.resample.test(ratio_matrix_test, Y_perm, 
                                                         statistic = DACOMP.TEST.NAME.LOG_FOLD_DIFFERENCE_IN_MEANS)
    
}
  

# computes p-values:
p.values.ratio.normalization = apply(stats_matrix_test, 2, function(x) mean(x >= x[1]))
  
head(sort(p.values.ratio.normalization),20)
which(p.values.ratio.normalization <= 0.02)
unname(tax_table(ps)[which(p.values.ratio.normalization <= 0.02),])

# "Raoultella"
# "Anaerostipes"
# "Sarcandra"

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Multiple comparison adjustment #
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  
### Lee, Forastiere, Miratix, Pillai (2017) method  for multiple comparison adjustment ### 
  
# STEP 1 to 3: recorded in "stats_matrix"
dim(stats_matrix_test) # the first row is the observed
  
# AT STEP 4: consider at each rep. T_rep as T_obs to calculate nrep p-values 
# for the hypothetical test statistics
  
# hyp_matrix <- stats_matrix_test # remove first row (obs.)
hyp_p_value <- matrix(NA, ncol = p_test, nrow = dim(stats_matrix_test)[1])
  
# based on value (hyp_obs) of each row
for (r in 1:nrep){
  if(verbose)
    if(r%% ceiling(nrep/10) == 1)
      cat(paste0('Testing rep : ',r,'/',nrep,' \n\r'))
  # calc. hypothetical p_value on each column of the matrix 
  hyp_p_value[r,] <- apply(stats_matrix_test, 2, function(x) mean(x >= x[r]))  
}
  
# for each rep. take the min. p_value
min_p_nrep <- apply(hyp_p_value, 1, function(x) min(x, na.rm = TRUE))
  
# calculate the proportion of min_p_nrep that is sm/eq. p_value (for obs.)
p_value_adj <- sapply(p.values.ratio.normalization, function(x) mean(min_p_nrep <= x))
head(sort(p_value_adj))

# p_adj_rejections <- which(p_value_adj <= 0.2)  
# unname(tax_table(ps)[p_adj_rejections,])
  




