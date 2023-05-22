#################################################################################################
#                                                                                               #
# Supplementary Material Data D4 for Bernstein et al.: Phylogenomics of Fresh and Formalin      #
# Specimens Resolves the Systematics of Old World Mud Snakes (Serpentes: Homalopsidae) and      #
# Expands Biogeographic Inference                                                               #              
#                                                                                               #   
#################################################################################################

### Set Working Directory ###

setwd("C:/Users/hugo/OneDrive - Indian Institute of Science/Indian Institute of Science (IISc)/Kartik Shanker Lab/Work/Homalopsidae/X1_Origin and Dispersal_Paper/Analysis/X8_APE/X0_R Working Directory/X4_Final_Genomics")


### Load Required Packages ###

rm(list=ls()); #clear R history
library(ape)
library(phytools)
library(geiger)
library(rlist)
library(phangorn)


### Import Nexus Tree File ###

tree <- read.nexus(file="Homalopsidae_Final_Genomics_Dated.tre")
tree


### Making consensus tree ### # (Do this in case you're using a combined tree rather than a consensus tree) #

#treemcc=maxCladeCred(tree,rooted=T)
#write.nexus(treemcc, file = "contree.nex", translate = TRUE)
#contree<-read.nexus(file="contree.nex")


## Case 1 - Terrestrial/Fossorial vs.Brackish vs. Freshwater - 0 = Fossorial; 1 = Freshwater; 2 = Brackish Water  ##

data <- read.csv(file="Supplementary_DataD5_APE_ASR_Case1.csv",row.names=1,header=TRUE)
habitat1 <- setNames(as.factor(data[,"Habitat"]),
                     rownames(data))
check <- name.check(tree,data)

fitace_homalopsidae_C1 <- ace(habitat1, tree, type = "discrete", method = "ML", CI = TRUE, model =  "ARD",
                              scaled = FALSE, kappa = 1, corStruct = NULL, ip = 0.1,
                              use.expm = FALSE, use.eigen = TRUE, marginal = FALSE)
save(fitace_homalopsidae_C1,file = "homalopsidae_ape_genomics_25012022_C1.rdata")

fitace_homalopsidae_C1 # to view details
round(fitace_homalopsidae_C1$lik.anc,3) # to round up the likelihood value at each internal node and viewing
dotTree(tree,habitat1,colors=setNames(c("brown","blue","green"),c("0","1","2")),ftype="i",fsize=0.7) # Putting pie on the tips representing the trait values
nodelabels(node=1:tree$Nnode+Ntip(tree),pie=fitace_homalopsidae_C1$lik.anc,piecol=c("brown","blue","green"),cex=0.3) # putting likelihood pies on the internal nodes
##_________________________________________________________##


## Case 2 - Terrestrial/Fossorial vs.Aquatic - 0 = Fossorial; 1 = Aquatic  ##

data <- read.csv(file="Supplementary_DataD6_APE_ASR_Case2.csv",row.names=1,header=TRUE)
habitat2 <- setNames(as.factor(data[,"Habitat"]),
                     rownames(data))
check <- name.check(tree,data)

fitace_homalopsidae_C2 <- ace(habitat2, tree, type = "discrete", method = "ML", CI = TRUE, model =  "ARD",
                              scaled = FALSE, kappa = 1, corStruct = NULL, ip = 0.1,
                              use.expm = FALSE, use.eigen = TRUE, marginal = FALSE)
save(fitace_homalopsidae_C2,file = "homalopsidae_ape_genomics_25012022_C2.rdata")

fitace_homalopsidae_C2 # to view details
round(fitace_homalopsidae_C2$lik.anc,3) # to round up the likelihood value at each internal node and viewing
dotTree(tree,habitat2,colors=setNames(c("brown","blue"),c("0","1")),ftype="i",fsize=0.7) # Putting pie on the tips representing the trait values
nodelabels(node=1:tree$Nnode+Ntip(tree),pie=fitace_homalopsidae_C2$lik.anc,piecol=c("brown","blue"),cex=0.3) # putting likelihood pies on the internal nodes

##_________________________________________________________##


## Case 3 - Terrestrial/Fossorial vs.Brackish vs. Freshwater WITH MISSING DATA REMOVED - 0 = Fossorial; 1 = Freshwater; 2 = Brackish Water  ##

data <- read.csv(file="Supplementary_DataD7_APE_ASR_Case3.csv",row.names=1,header=TRUE)
habitat3 <- setNames(as.factor(data[,"Habitat"]),
                     rownames(data))
check <- name.check(tree,data)

tree3 <-drop.tip(tree,tip=check$tree_not_data) # (Use this if dropping missing data present in tree but not in .csv file) #

fitace_homalopsidae_C3 <- ace(habitat3, tree3, type = "discrete", method = "ML", CI = TRUE, model =  "ARD",
                              scaled = FALSE, kappa = 1, corStruct = NULL, ip = 0.1,
                              use.expm = FALSE, use.eigen = TRUE, marginal = FALSE)
save(fitace_homalopsidae_C3,file = "homalopsidae_ape_genomics_25012022_C3.rdata")

fitace_homalopsidae_C3 # to view details
round(fitace_homalopsidae_C3$lik.anc,3) # to round up the likelihood value at each internal node and viewing
dotTree(tree3,habitat3,colors=setNames(c("brown","blue","green"),c("0","1","2")),ftype="i",fsize=0.7) # Putting pie on the tips representing the trait values
nodelabels(node=1:tree3$Nnode+Ntip(tree3),pie=fitace_homalopsidae_C3$lik.anc,piecol=c("brown","blue","green"),cex=0.3) # putting likelihood pies on the internal nodes

##_________________________________________________________##

## Case 4 - Terrestrial/Fossorial vs.Brackish vs. Freshwater WITH MISSING DATA & OUTGROUPS REMOVED - 0 = Fossorial; 1 = Freshwater; 2 = Brackish Water  ##

data <- read.csv(file="Supplementary_DataD8_APE_ASR_Case4.csv",row.names=1,header=TRUE)
habitat5 <- setNames(as.factor(data[,"Habitat"]),
                     rownames(data))
check <- name.check(tree,data)

tree5 <-drop.tip(tree,tip=check$tree_not_data) # (Use this if dropping missing data present in tree but not in .csv file) #

fitace_homalopsidae_C5 <- ace(habitat5, tree5, type = "discrete", method = "ML", CI = TRUE, model =  "ARD",
                              scaled = FALSE, kappa = 1, corStruct = NULL, ip = 0.1,
                              use.expm = FALSE, use.eigen = TRUE, marginal = FALSE)
save(fitace_homalopsidae_C5,file = "homalopsidae_ape_genomics_25012022_C4.rdata")

fitace_homalopsidae_C5 # to view details
round(fitace_homalopsidae_C5$lik.anc,3) # to round up the likelihood value at each internal node and viewing
dotTree(tree5,habitat5,colors=setNames(c("brown","blue","green"),c("0","1","2")),ftype="i",fsize=0.7) # Putting pie on the tips representing the trait values
nodelabels(node=1:tree5$Nnode+Ntip(tree5),pie=fitace_homalopsidae_C5$lik.anc,piecol=c("brown","blue","green"),cex=0.3) # putting likelihood pies on the internal nodes

##_________________________________________________________##



## Case 5 - Terrestrial/Fossorial vs.Aquatic WITH MISSING DATA & OUTGROUPS REMOVED - 0 = Fossorial; 1 = Aquatic;  ##

data <- read.csv(file="Supplementary_DataD9_APE_ASR_Case5.csv",row.names=1,header=TRUE)
habitat6 <- setNames(as.factor(data[,"Habitat"]),
                     rownames(data))
check <- name.check(tree,data)

tree6 <-drop.tip(tree,tip=check$tree_not_data) # (Use this if dropping missing data present in tree but not in .csv file) #

fitace_homalopsidae_C6 <- ace(habitat6, tree6, type = "discrete", method = "ML", CI = TRUE, model =  "ARD",
                              scaled = FALSE, kappa = 1, corStruct = NULL, ip = 0.1,
                              use.expm = FALSE, use.eigen = TRUE, marginal = FALSE)
save(fitace_homalopsidae_C6,file = "homalopsidae_ape_genomics_25012022_C5.rdata")

fitace_homalopsidae_C6 # to view details
round(fitace_homalopsidae_C6$lik.anc,3) # to round up the likelihood value at each internal node and viewing
dotTree(tree6,habitat6,colors=setNames(c("brown","blue"),c("0","1")),ftype="i",fsize=0.7) # Putting pie on the tips representing the trait values
nodelabels(node=1:tree6$Nnode+Ntip(tree6),pie=fitace_homalopsidae_C6$lik.anc,piecol=c("brown","blue"),cex=0.3) # putting likelihood pies on the internal nodes

##_________________________________________________________##
