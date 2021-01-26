library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(breakaway)
library(dplyr)
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

################################################
### 1. ESTIMATE THE RICHNESS FOR EACH SAMPLE ###
################################################

# calculate plug-in richness
rich <- sample_richness(ps_prune)

## estimate richness with breakaway function  
ba <- breakaway(ps_prune)
ba[[1]]
# plot(ba, ps_prune, color = "W")

## retrieve statistics for inference
sample_data(ps_prune)[,"breakaway_W"] <- summary(ba)$estimate
sample_data(ps_prune)[,"ba_error_W"] <- summary(ba)$error
sample_data(ps_prune)[,"lower"] <- summary(ba)$lower
sample_data(ps_prune)[,"upper"] <- summary(ba)$upper
sample_data(ps_prune)[,"sample_counts"] <- sample_sums(ps_prune)
sample_data(ps_prune)[,"richness"] <- summary(rich)$estimate
sample_data(ps_prune)[,"id"] <- rownames(sample_data(ps_prune))

ggplot(sample_data(ps_prune), aes(x = id)) + 
  geom_point(aes(y = breakaway_W), colour = "red") + 
  # geom_errorbar(aes(ymin=lower, ymax=upper)) +
  geom_point(aes(y = richness), colour = "blue", alpha = .5) +
  xlab('Samples') + ylab("BA estimate")

g_PM <- ggplot(sample_data(ps_prune_out), aes(color = factor(W), y = breakaway_W)) +
  geom_boxplot(alpha = .5) + ylab('breakaway richness measure') +
  scale_x_discrete(name = "") +
  scale_colour_manual(values = c("darkgreen","blue4"), limits=c("1","0"), name ="Smoking", labels = c("No","Yes")) +
  theme(legend.position = "top", legend.key.size =  unit(0.1, "in")) +
  guides(color=guide_legend(nrow=2,byrow=TRUE))

#######################################################
### 2. RANDOMIZATION-TEST WITH BETTA TEST STATISTIC ###
#######################################################
x <- cbind(1, 
          # sample_data(ps_prune)$sex,
          sample_data(ps_prune)$W)
head(sample_data(ps_prune)$breakaway_W)
head(sample_data(ps_prune)$ba_error_W)

reg <- betta(sample_data(ps_prune)$breakaway_W,
             sample_data(ps_prune)$ba_error_W, 
             X = x)
reg$table
estim_obs <- reg$table[2,1]

## W matrix 
dim(W_paired_smoke)
W_paired_smoke <- W_paired_smoke[order(rownames(W_paired_smoke)), ]

# set the number of randomizations
nrep <- ncol(W_paired_smoke)/1000
# nrep <- ncol(W_paired)/100

# create a matrix where the t_rand will be saved
t_array <- NULL

for(i in 2:nrep){
  print(i)
  x = cbind(1, 
            # sample_data(ps_prune)$u3tcigsmk, 
            # sample_data(ps_prune)$u3csex
            W_paired_smoke[,i])
            # W_paired[,i]) 
            
  
  reg = betta(sample_data(ps_prune)$breakaway_W,
              sample_data(ps_prune)$ba_error_W, X = x)
  
  # fill t_array
  t_array[i] = reg$table[2,1] 
}

## calculate p_value
p_value <- mean(t_array >= estim_obs, na.rm = TRUE)
p_value

## plot distribution of test-statistic
hist(t_array, breaks = 30, main = "", xlab = "breakaway beta (Smoking)")
abline(v = estim_obs, col = 'red', lwd = 2, lty = 2)

