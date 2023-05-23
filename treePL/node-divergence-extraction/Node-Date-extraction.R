#set working directory
setwd("D:/Documents/Homalopsidae_Phylogenomics/R_Burbrink2020_TreeDates/Dated_Trees_Fast_Bootstrapped[1]/Dated_Trees_Fast_Bootstrapped/")
getwd()

#Load packages
library(phytools)
library(ape)
library(HDInterval)

##read the dated trees into R
list.files()->l



##remove some r script stuff at the top
l[-c(1:2)]->l
lapply(l,function(x) read.tree(x))->t



##get node number and make sure its all the same (this was bootstrapped on the same tree). I am getting the stem date between a Homalopsid and Cyclocorus
unlist(lapply(c(1:length(t)),function(x) getMRCA(t[[x]],c("Elapidae_Micrurus_fulvius_I8101", "Colubridae_Simophis_rhinostoma_I12212"))))->node


lapply(t,function(x) branching.times(x))->bt
lapply(c(1:length(bt)),function(x) which(row.names(as.data.frame(bt[[x]]))==node[[x]]))->n2
lapply(c(1:length(bt)),function(x) data.frame(bt[[x]])[n2[[x]],])->new_date
unlist(new_date)->ND



summary(ND)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 37.37   47.25   47.64   47.35   48.03   50.23 

quantile(ND,0.95)


quantile(ND,c(0.025,0.975))
# 2.5%    97.5% 
#   37.59323 49.32168 

plot(density(ND))

#Getting the 95% HPD
hdi(ND) 

##distribution of dates look s a little choppy but certainly centered on ~44-50 my with lots of stuff in the tails (35 and 55 my)

citation("phytools")
citation("ape")
citation("HDInterval")
