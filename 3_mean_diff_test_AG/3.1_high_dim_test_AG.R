library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")

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

# ## agglomerate to Species ##
# ps_Species <- tax_glom(ps, taxrank = "Species", NArm = FALSE)
# vec_taxa_Sp <- apply(otu_table(ps_Species), 2, function(x) sum(x > 0, na.rm = TRUE))
# length(which(vec_taxa_Sp > samp_perc))
# ps_Species <- prune_taxa(vec_taxa_Sp >= 0, ps_Species)
# ps_Species

## agglomerate to Genus ##
ps_Genus <- tax_glom(ps_prune, taxrank = "Genus", NArm = FALSE)
vec_taxa_Gen <- apply(otu_table(ps_Genus), 1, function(x) sum(x > 0, na.rm = TRUE))
table(vec_taxa_Gen)
ps_Genus_prune <- prune_taxa(vec_taxa_Gen >= 0, ps_Genus)
ps_Genus_prune

ps <- ps_Genus_prune

###########################################
### 1. TEST-STATISTIC FROM CAO, LIN, LI ###
###########################################
# Yuanpei Cao,  Wei Lin,  Hongzhe Li (2017) - Two-sample tests of high-dimensional means for compositional data 

# Filter the data
n<-dim(otu_table(ps))[2]
p<-dim(otu_table(ps))[1]

# Add 0.5 on zeros
otu_table(ps)[otu_table(ps)==0] <- 0.5

# Separate the sample into two groups
ps_control = subset_samples(ps, sample_data(ps)$W == 0)
ps_treated = subset_samples(ps, sample_data(ps)$W == 1)

X_c = t(otu_table(ps_control)); dim(X_c)
X_c = as(X_c, "matrix")

X_t = t(otu_table(ps_treated)); dim(X_t)
X_t = as(X_t, "matrix")

nx <- dim(X_c)[1]
ny <- dim(X_t)[1]

x <- X_c/(rowSums(X_c)%*%matrix(1,1,p))
y <- X_t/(rowSums(X_t)%*%matrix(1,1,p))

log_x<-log(x)
log_y<-log(y)

clog_TX <- log_x-1/p*rowSums(log_x)%*%matrix(1,1,p)
clog_TY <- log_y-1/p*rowSums(log_y)%*%matrix(1,1,p)

# Mean
# x_mean <- colSums(x)/nx
# y_mean <- colSums(y)/ny
# 
# log_x_mean <- colSums(log_x)/nx
# log_y_mean <- colSums(log_y)/ny

clr_x_mean <- colSums(clog_TX)/nx
clr_y_mean <- colSums(clog_TY)/ny

# Variance
# x_var <- diag(var(x))*(nx-1)/nx
# y_var <- diag(var(y))*(ny-1)/ny
# x_stat_var <- (x_var*nx + y_var*ny)/(nx*ny)
# 
# log_x_var <- diag(var(log_x))*(nx-1)/nx
# log_y_var <- diag(var(log_y))*(ny-1)/ny
# log_x_stat_var <- (log_x_var*nx + log_y_var*ny)/(nx*ny)

clr_x_var <- diag(var(clog_TX))*(nx-1)/nx
clr_y_var <- diag(var(clog_TY))*(ny-1)/ny
clr_x_stat_var <- (clr_x_var*nx + clr_y_var*ny)/(nx*ny)

# p-values
# x_stat <- max((x_mean/sqrt(x_stat_var) - y_mean/sqrt(x_stat_var))^2)
# p_x <- 1-exp(-1/sqrt(pi)*exp(-(x_stat-(2*log(p)-log(log(p))))/2))
# log_x_stat <- max((log_x_mean/sqrt(log_x_stat_var) - log_y_mean/sqrt(log_x_stat_var))^2)
# p_logx <- 1-exp(-1/sqrt(pi)*exp(-(log_x_stat-(2*log(p)-log(log(p))))/2))
clr_x_stat <- max((clr_x_mean/sqrt(clr_x_stat_var) - clr_y_mean/sqrt(clr_x_stat_var))^2)
p_clrx <- 1-exp(-1/sqrt(pi)*exp(-(clr_x_stat-(2*log(p)-log(log(p))))/2))
clr_x_stat; p_clrx

#######################################################
### 2. RANDOMIZATION-TEST MEAN DIFF. TEST STATISTIC ###
#######################################################

matrix_otu_bind <- rbind(clog_TX, clog_TY)
dim(matrix_otu_bind)

# set the number of randomizations
W_paired <- W_paired_smoke[order(rownames(W_paired_smoke)), ]
nrep <- ncol(W_paired)/100

# create vector to be save 
# first is the obs. test stat
Tarray = clr_x_stat

for (i in 1:nrep){
  W_suffle <- W_paired[,i]
  
  clog_TX_rand <- matrix_otu_bind[W_suffle == 0,]
  clog_TY_rand <- matrix_otu_bind[W_suffle == 1,]
  
  clr_x_mean_rand <- colSums(clog_TX_rand)/nx
  clr_y_mean_rand <- colSums(clog_TY_rand)/ny
  
  clr_x_var_rand <- diag(var(clog_TX_rand))*(nx-1)/nx
  clr_y_var_rand <- diag(var(clog_TY_rand))*(ny-1)/ny
  clr_x_stat_var_rand <- (clr_x_var_rand*nx + clr_y_var_rand*ny)/(nx*ny)
  
  clr_x_stat_rand <- max((clr_x_mean_rand/sqrt(clr_x_stat_var_rand) - clr_y_mean_rand/sqrt(clr_x_stat_var_rand))^2)
  
  Tarray = c(Tarray, clr_x_stat_rand)
  
  print(i)
}

hist(Tarray, main = NULL, xlab = NULL, breaks = 100)
abline(v=clr_x_stat, lty = 2, lwd = 2)
pval = mean(Tarray >= clr_x_stat)
pval

# estimate: 50.0806
# p-value: 0.000999001
