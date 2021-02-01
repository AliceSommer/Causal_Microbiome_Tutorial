
# install develop2 version of NetCoMi
# devtools::install_github("stefpeschel/NetCoMi", ref = "develop2", force = TRUE)

library(NetCoMi)
library(metagMisc)
library(phyloseq)
library(igraph)
library(plyr)
library(ForceAtlas2)
library(pals)
library(RColorBrewer)
library(corrplot)

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
# necessary for NetCoMi analyses 
matched_df$W <- as.factor(matched_df$W)
sample_df <- matched_df[order(matched_df$W, matched_df$pair_nb),]
otu_table(agdata) <- otu_table(agdata)[,rownames(sample_df)]

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

## agglomerate to Family ##
ps_Fam <- tax_glom(ps_prune, taxrank = "Family", NArm = FALSE)
vec_taxa_Fam <- apply(otu_table(ps_Fam), 1, function(x) sum(x > 0, na.rm = TRUE))
table(vec_taxa_Fam)
ps_Fam_prune <- prune_taxa(vec_taxa_Fam >= 0, ps_Fam)
ps_Fam_prune

ps <- ps_Fam_prune

# check order of matching
cols <- c("sex", "pair_nb", "W")
tail(sample_data(ps)[,cols])

###########
# NetCoMi # 
###########

gut_split <- metagMisc::phyloseq_sep_variable(ps, "W")

net_W <- netConstruct(gut_split$`0`, gut_split$`1`, verbose = 2,
                           filtTax = "highestVar",
                           filtTaxPar = list(highestVar = 50),
                           measure = "spieceasi",
                           measurePar = list(method = "glasso",
                                             nlambda=20,
                                             pulsar.params=list(rep.num=20)),
                           normMethod = "none", zeroMethod = "none",
                           sparsMethod = "none", seed = 123456, 
                           matchDesign = c(1,1),
                           jointPrepro = FALSE )


props_W <- netAnalyze(net_W, clustMethod = "cluster_fast_greedy", 
                      connectivity = FALSE,
                      avDissIgnoreInf = TRUE,
                      weightClustCoef = FALSE )

plot(props_W, sameLayout = TRUE, layoutGroup = 1,
     nodeSize = "eigenvector", cexNodes = 1.5, cexLabels = 1.8,
     groupNames = c("Non-Smokers", "Daily Smokers"),
     hubBorderCol  = "gray40")

# this matrix can be stored externally
# and computed in parallel with permGroupMat
assoPerm <- createAssoPerm(props_W, nPerm = 2, 
                               computeAsso = TRUE,
                               storeCountsPerm = FALSE,
                               fileStoreAssoPerm = './5_networks_AG/assoPerm',
                               seed = 123456, append = FALSE)

comp_W <- netCompare(props_W, 
                     permTest = TRUE, nPerm = 2,
                     storeAssoPerm = FALSE,
                     fileLoadAssoPerm = "./5_networks_AG/assoPerm",
                     storeCountsPerm = FALSE, seed = 123456,
                     returnPermProps = TRUE)

summary(comp_W, showCentr = c("degree", "eigen"), numbTaxa = 5)

diff_net <- diffnet(net_W, nPerm = 2, fileLoadAssoPerm = "./5_networks_AG/assoPerm")
