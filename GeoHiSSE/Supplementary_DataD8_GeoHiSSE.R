#################################################################################################
#                                                                                               #
# Supplementary Material Data D1 for Bernstein et al.: GeoHiSSE for Homalopsidae Phylogenomics  #
#                                                                                               #   
#################################################################################################

# install required packages (commented out if already done)
# install.packages("diversitree")
# library(devtools)
# install_github(repo = "thej022214/hisse", ref = "master")
# install.packages("devtools")



#Load libraries
suppressWarnings(library(hisse))
suppressWarnings(library(diversitree))


# Set the working directory
setwd("D:/Documents/Homalopsidae_Phylogenomics/GeoHisse/temp_new/")

#This code is based on the  Tutorial "Running GeoHiSSE" by Daniel Caetano and Jeremy M. Beaulieu

###Read in the time-calibrated phylogeny
phy <- read.tree("Supplementary_DataD9_time-calibrated-tree.tre")
phy$node.labels <- NULL
sim.dat <- read.csv("Supplementary_DataD10_homa_aq_geohisse.csv", stringsAsFactors=FALSE)

#trim the tips you don't need (outgroups + the fangless homalopsids)
phy<-drop.tip(phy, "Bothrops_moojeni_SRS2115278")
phy<-drop.tip(phy, "Bothrops_pauloensis_SRS2115295")
phy<-drop.tip(phy, "Micrurus_brasiliensis_SRS2115293")
phy<-drop.tip(phy, "Chironius_exoletus_SRS2115290")
phy<-drop.tip(phy, "Philodryas_olfersii_SRS2115301")
phy<-drop.tip(phy, "Brachyorrhos_albus_FMNH134323")

###Setting up the Models: see Supplementary Material Text S1 for further explanation on models.

#Set the sampling fraction
samp.frac.fresh <- (18/26) #our sampling/a total of 26 freshwater homalopsids
samp.frac.brack <- (9/17)  #our sampling/a total of 17 brackish homalopsids
samp.frac.both  <- (3/3)   #our sampling/a total of 3  freshwater/brackish homalopsids

samp.frac <- c(samp.frac.fresh,samp.frac.brack,samp.frac.both)

## Model 1 - Dispersal parameters vary only, no range-dependent diversification.
turnover <- c(1,1,0)
eps <- c(1,1)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=0)
trans.rate.mod <- ParEqual(trans.rate, c(1,2))
mod1 <- GeoHiSSE(phy = phy, data = sim.dat, f=samp.frac,
                 turnover=turnover, eps=eps,
                 hidden.states=FALSE, trans.rate=trans.rate.mod,
                 turnover.upper=100, trans.upper=10)

# To conduct a canonical GeoSSE model, where range evolution affects diversification, we add back the turnover
# rate for the widespread range, such that there are three turnover parameters and two extinction fraction
# parameters estimated:
## Model 2. Canonical GeoSSE model, range effect on diversification
turnover <- c(1,2,3)
eps <- c(1,1)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=0)
trans.rate.mod <- ParEqual(trans.rate, c(1,2))
mod2 <- GeoHiSSE(phy = phy, data = sim.dat, f=samp.frac,
                 turnover=turnover, eps=eps,
                 hidden.states=FALSE, trans.rate=trans.rate.mod,
                 turnover.upper=100, trans.upper=10)

# Models 3 and 4 below each have 2 hidden states. In this case the models will be more complex. First, we
# will show how to set up a range-independent model of diversification. Remember, as with mod1 above, if
# diversification is independent of range-evolution we must remove the turnover rate for the widespread range.
# Again, internally, GeoSSE() will recognize this (if assume.cladogenetic=TRUE) and simply set s01 to be
# equal the rate for s00 and s11.
## Model 3. GeoHiSSE model with 1 hidden trait, no range-dependent diversification.
## Note below how parameters vary among hidden classes but are the same within each
## hidden class.
turnover <- c(1,1,0,2,2,0)
eps <- c(1,1,1,1)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=1, make.null=TRUE)
trans.rate.mod <- ParEqual(trans.rate, c(1,2))
mod3 <- GeoHiSSE(phy = phy, data = sim.dat, f=samp.frac,
                 turnover=turnover, eps=eps,
                 hidden.states=TRUE, trans.rate=trans.rate.mod,
                 turnover.upper=100, trans.upper=10)

# Finally, if we want to fit a GeoHiSSE model we would do the following:
## Model 4. GeoHiSSE model with 1 hidden trait, range-dependent diversification.
turnover <- c(1,2,3,4,5,6)
eps <- c(1,1,1,1)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=1)
trans.rate.mod <- ParEqual(trans.rate, c(1,2,3,4))  # I added in ,3,4 here -SMH
mod4 <- GeoHiSSE(phy = phy, data = sim.dat, f=samp.frac,
                 turnover=turnover, eps=eps,
                 hidden.states=TRUE, trans.rate=trans.rate.mod,
                 turnover.upper=100, trans.upper=10)


# Model 4a: Make a state-independent version of mod4
turnover <- c(1,1,0,2,2,0,3,3,0,4,4,0,5,5,0,6,6,0)
eps <- rep(1, 12)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=5, make.null=TRUE)
trans.rate.mod <- ParEqual(trans.rate, c(1,2))
mod4a <- GeoHiSSE(phy = phy, data = sim.dat, f=samp.frac,
                 turnover=turnover, eps=eps,
                 hidden.states=TRUE, trans.rate=trans.rate.mod,
                 turnover.upper=100, trans.upper=10)

# We will also show how to fit a complementary set of models that remove the cladogenetic effect entirely, such
# that all changes occur along branches (i.e., anagenetic change). This requires the removal of the turnover
# rate for lineages in the widespread range and ensuring that range contraction is distinct from the extinction
# of endemics:
## Model 5. MuSSE-like model with no hidden trait, no cladogenetic effects. simplest state dependent model
turnover <- c(1,2,0)
eps <- c(1,1)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=0, make.null=FALSE, 
                                    separate.extirpation = TRUE)
trans.rate.mod <- ParEqual(trans.rate, c(1,2))                         
trans.rate.mod <- ParEqual(trans.rate.mod, c(2,3))
mod5 <- GeoHiSSE(phy = phy, data = sim.dat, f=samp.frac,
                 turnover=turnover, eps=eps,
                 hidden.states=FALSE, trans.rate=trans.rate.mod,
                 turnover.upper=100, trans.upper=10,
                 assume.cladogenetic = FALSE)

# Model 6 - a model like model 1, but including separate.extirpation=TRUE
turnover <- c(1,1,0)
eps <- c(1,1)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=0, separate.extirpation = TRUE)
trans.rate.mod <- ParEqual(trans.rate, c(1,2))
trans.rate.mod <- ParEqual(trans.rate.mod, c(2,3))
mod6 <- GeoHiSSE(phy = phy, data = sim.dat, f=samp.frac,
                 turnover=turnover, eps=eps,
                 hidden.states=FALSE, trans.rate=trans.rate.mod,
                 turnover.upper=100, trans.upper=10)


# Model 7 - a separate.extirpation version of model 3
turnover <- c(1,1,0,2,2,0)
eps <- c(1,1,1,1)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=1, make.null=TRUE, separate.extirpation = TRUE)
trans.rate.mod <- ParEqual(trans.rate, c(3,4))
trans.rate.mod <- ParEqual(trans.rate.mod, c(1,2))
mod7 <- GeoHiSSE(phy = phy, data = sim.dat, f=samp.frac,
                 turnover=turnover, eps=eps,
                 hidden.states=TRUE, trans.rate=trans.rate.mod,
                 turnover.upper=100, trans.upper=10)

## Model 8 - a separate.extirpation version of model 4
turnover <- c(1,2,0,3,4,0)
eps <- c(1,1,1,1)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=1, separate.extirpation = TRUE)
trans.rate.mod <- ParEqual(trans.rate, c(1,2))
trans.rate.mod <- ParEqual(trans.rate.mod, c(2,3))
trans.rate.mod <- ParEqual(trans.rate.mod, c(3,4))
trans.rate.mod <- ParEqual(trans.rate.mod, c(4,5))
trans.rate.mod <- ParEqual(trans.rate.mod, c(2,4))
trans.rate.mod <- ParEqual(trans.rate.mod, c(1,3))
mod8 <- GeoHiSSE(phy = phy, data = sim.dat, f=samp.frac,
                 turnover=turnover, eps=eps,
                 hidden.states=TRUE, trans.rate=trans.rate.mod,
                 turnover.upper=100, trans.upper=10)

## model 9 - a separate.extirpation=TRUE of model 4a
turnover <- c(1,1,0,2,2,0,3,3,0,4,4,0)
eps <- rep(1, 8)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=3, make.null=TRUE, separate.extirpation = TRUE)
trans.rate.mod <- ParEqual(trans.rate, c(1,2))
trans.rate.mod <- ParEqual(trans.rate.mod, c(2,3))
mod9 <- GeoHiSSE(phy = phy, data = sim.dat, f=samp.frac,
                  turnover=turnover, eps=eps,
                  hidden.states=TRUE, trans.rate=trans.rate.mod,
                  turnover.upper=100, trans.upper=10)

######### Now let's start making variants of the separate.extirpation=TRUE models that have more transition rates

# Model 10 - same as 5 but with 4 trans rates
turnover <- c(1,2,0)
eps <- c(1,1)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=0, make.null=FALSE, separate.extirpation = TRUE)
mod10 <- GeoHiSSE(phy = phy, data = sim.dat, f=samp.frac,
                 turnover=turnover, eps=eps,
                 hidden.states=FALSE, trans.rate=trans.rate,
                 turnover.upper=100, trans.upper=10,
                 assume.cladogenetic = FALSE)

# Model 11 - same as mod 6 but with 4 trans rates
turnover <- c(1,1,0)
eps <- c(1,1)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=0, separate.extirpation = TRUE)
mod11 <- GeoHiSSE(phy = phy, data = sim.dat, f=samp.frac,
                 turnover=turnover, eps=eps,
                 hidden.states=FALSE, trans.rate=trans.rate,
                 turnover.upper=100, trans.upper=10)

# Model 12 - same as 7 but with asymmetric trans rates
turnover <- c(1,1,0,2,2,0)
eps <- c(1,1,1,1)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=1, make.null=TRUE, separate.extirpation = TRUE)
mod12 <- GeoHiSSE(phy = phy, data = sim.dat, f=samp.frac,
                 turnover=turnover, eps=eps,
                 hidden.states=TRUE, trans.rate=trans.rate,
                 turnover.upper=100, trans.upper=10)

# Model 13 - same as 8 but with asymmetric trans rates
turnover <- c(1,2,0,3,4,0)
eps <- c(1,1,1,1)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=1, separate.extirpation = TRUE)
trans.rate.mod <- ParEqual(trans.rate, c(3,7))
trans.rate.mod <- ParEqual(trans.rate.mod, c(4,7))
trans.rate.mod <- ParEqual(trans.rate.mod, c(1,5))
trans.rate.mod <- ParEqual(trans.rate.mod, c(2,5))
mod13 <- GeoHiSSE(phy = phy, data = sim.dat, f=samp.frac,
                 turnover=turnover, eps=eps,
                 hidden.states=TRUE, trans.rate=trans.rate.mod,
                 turnover.upper=100, trans.upper=10)

# Model 14 - same as 9 but with asymmetric trans rates
turnover <- c(1,1,0,2,2,0,3,3,0,4,4,0)
eps <- rep(1, 8)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=3, make.null=TRUE, separate.extirpation = TRUE)
mod14 <- GeoHiSSE(phy = phy, data = sim.dat, f=samp.frac,
                 turnover=turnover, eps=eps,
                 hidden.states=TRUE, trans.rate=trans.rate,
                 turnover.upper=100, trans.upper=10)

# Compare the AIC weights across all models
# Akaike weights are important to evaluate the relative importance of each of the models to explain the variation
# observed in the data. This quantity takes into account pennalties associated to the number of free parametes.
# Models with higher weight show better fit to the data and, as a result, have more weight when performing
# model averaging (see below).

weights <- GetAICWeights(list(model1 = mod1, model2 = mod2, model3 = mod3, model4 = mod4, model4a = mod4a, model5=mod5,
  model6=mod6, model7 = mod7, model8=mod8, model9=mod9, model10 = mod10, model11 = mod11, model12 = mod12,
  model13 =mod13, model14 = mod14), criterion="AIC")

sort(weights)

## As the number of models in the set grows, naming each model in the set can become hard.
## So one can use a list (created by some automated code) as an input also:
list.geohisse <- list(model1 = mod1, model2 = mod2, model3 = mod3, model4 = mod4, model4a = mod4a, model5 = mod5, model6 = mod6, model7 = mod7, model8 = mod8,
                      model9 = mod9, model10 = mod10, model11 = mod11, model12 = mod12, model13 = mod13, model14 = mod14)
GetAICWeights(list.geohisse, criterion="AIC")
AICWeights <- GetAICWeights(list.geohisse, criterion="AIC")
#write out the output
write.csv(AICWeights, file = "AICWeights.csv")

### Model averaging and plotting.

# Now we can model average the results. Note that this step will reflect the Akaike model weights that we
# computed above.

recon.mod1 <- MarginReconGeoSSE(phy = mod1$phy, data = mod1$data, f = mod1$f,
                                pars = mod1$solution, hidden.states = 1,
                                root.type = mod1$root.type, root.p = mod1$root.p,
                                AIC = mod1$AIC, n.cores = 1)
recon.mod2 <- MarginReconGeoSSE(phy = mod2$phy, data = mod2$data, f = mod2$f,
                                pars = mod2$solution, hidden.states = 1,
                                root.type = mod2$root.type, root.p = mod2$root.p,
                                AIC = mod2$AIC, n.cores = 1)
recon.mod3 <- MarginReconGeoSSE(phy = mod3$phy, data = mod3$data, f = mod3$f,
                                pars = mod3$solution, hidden.states = 2,
                                root.type = mod3$root.type, root.p = mod3$root.p,
                                AIC = mod3$AIC, n.cores = 1)
recon.mod4 <- MarginReconGeoSSE(phy = mod4$phy, data = mod4$data, f = mod4$f,
                                pars = mod4$solution, hidden.states = 2,
                                root.type = mod4$root.type, root.p = mod4$root.p,
                                AIC = mod4$AIC, n.cores = 1)
recon.mod4a <- MarginReconGeoSSE(phy = mod4a$phy, data = mod4a$data, f = mod4a$f,
                                pars = mod4a$solution, hidden.states = 2,
                                root.type = mod4a$root.type, root.p = mod4a$root.p,
                                AIC = mod4a$AIC, n.cores = 1)
recon.mod5 <- MarginReconGeoSSE(phy = mod5$phy, data = mod5$data, f = mod5$f,
                                pars = mod5$solution, hidden.states = 2,
                                root.type = mod5$root.type, root.p = mod5$root.p,
                                AIC = mod5$AIC, n.cores = 1)
recon.mod6 <- MarginReconGeoSSE(phy = mod6$phy, data = mod6$data, f = mod6$f,
                                pars = mod6$solution, hidden.states = 2,
                                root.type = mod6$root.type, root.p = mod6$root.p,
                                AIC = mod6$AIC, n.cores = 1)
recon.mod7 <- MarginReconGeoSSE(phy = mod7$phy, data = mod7$data, f = mod7$f,
                                pars = mod7$solution, hidden.states = 2,
                                root.type = mod7$root.type, root.p = mod7$root.p,
                                AIC = mod7$AIC, n.cores = 1)
recon.mod8 <- MarginReconGeoSSE(phy = mod8$phy, data = mod8$data, f = mod8$f,
                                pars = mod8$solution, hidden.states = 2,
                                root.type = mod8$root.type, root.p = mod8$root.p,
                                AIC = mod8$AIC, n.cores = 1)
recon.mod9 <- MarginReconGeoSSE(phy = mod9$phy, data = mod9$data, f = mod9$f,
                                pars = mod9$solution, hidden.states = 2,
                                root.type = mod9$root.type, root.p = mod9$root.p,
                                AIC = mod9$AIC, n.cores = 1)
recon.mod10 <- MarginReconGeoSSE(phy = mod10$phy, data = mod10$data, f = mod10$f,
                                pars = mod10$solution, hidden.states = 2,
                                root.type = mod10$root.type, root.p = mod10$root.p,
                                AIC = mod10$AIC, n.cores = 1)
recon.mod11 <- MarginReconGeoSSE(phy = mod11$phy, data = mod11$data, f = mod11$f,
                                pars = mod11$solution, hidden.states = 2,
                                root.type = mod11$root.type, root.p = mod11$root.p,
                                AIC = mod11$AIC, n.cores = 1)
recon.mod12 <- MarginReconGeoSSE(phy = mod12$phy, data = mod12$data, f = mod12$f,
                                pars = mod12$solution, hidden.states = 2,
                                root.type = mod12$root.type, root.p = mod12$root.p,
                                AIC = mod12$AIC, n.cores = 1)
recon.mod13 <- MarginReconGeoSSE(phy = mod13$phy, data = mod13$data, f = mod13$f,
                                pars = mod13$solution, hidden.states = 2,
                                root.type = mod13$root.type, root.p = mod13$root.p,
                                AIC = mod13$AIC, n.cores = 1)
recon.mod14 <- MarginReconGeoSSE(phy = mod14$phy, data = mod14$data, f = mod14$f,
                                pars = mod14$solution, hidden.states = 2,
                                root.type = mod14$root.type, root.p = mod14$root.p,
                                AIC = mod14$AIC, n.cores = 1)

# Now that we have the AIC associated with each model and their reconstruction across the nodes of the tree
# we can compute the model average:

recon.models <- list(recon.mod1, recon.mod2, recon.mod3, recon.mod4, recon.mod4a, recon.mod5, recon.mod6, recon.mod7, recon.mod8, recon.mod9,
                     recon.mod10, recon.mod11, recon.mod12, recon.mod13, recon.mod14)
model.ave.rates <- GetModelAveRates(x = recon.models, type = "tips")
#write out the output
write.csv(model.ave.rates, file = "model.ave.rates.csv")

# The result of the reconstruction is a matrix with the parameter estimates for each of the tips species averaged
# over all models. Note that for the GeoSSE model there is no "extinction" parameter associated with
# widespread (01) lineages. Also not that one can change the type of model averaging (between tips, nodes,
# and both) when calling the GetModelAveRates function.

# Finally, we can plot the use the resulting data matrix to make a plot of the results.

plot.geohisse.states(x = recon.models, rate.param = "net.div", type = "phylogram",
                     show.tip.label = FALSE, legend = TRUE,  legend.kernel = "hist", state.colors = c("light blue", "royal blue", "black"), edge.width = 6)

###References
# Caetano, D.S., B.C. O'Meara, and J.M. Beaulieu. 2018. Hidden state models improve state-dependent
# diversification approaches, including biogeographic models. Evolution, 72:2308-2324.



#######################################################################################
                                                                                      #
# Below is the code for running GeoHiSSE on the cyt-b tree (fresh+degraded specimens) #
                                                                                      #
#######################################################################################


#install required packages (commented out if already done)
# install.packages("diversitree")
# library(devtools)
# install_github(repo = "thej022214/hisse", ref = "master")
# install.packages("devtools")



#Load libraries
suppressWarnings(library(hisse))
suppressWarnings(library(diversitree))


# Set the working directory
setwd("D:/Documents/Homalopsidae_Phylogenomics/GeoHisse_FreshForm-CYTB/temp_new/")

#Start with the Tutorial:Running GeoHiSSE by Daniel Caetano and Jeremy M. Beaulieu

###Simulating a range-independent process

# We will simulate a phylogenetic tree using neutral geographical ranges, and incorporate two different rates of
# diversification. Thus, the correct process here is: "rates of diversification vary independently of the geographic
# ranges".

# Generate a list with the parameters of the model:
# pars <- SimulateGeoHiSSE(hidden.traits = 1, return.GeoHiSSE_pars = TRUE)
# pars

# The object pars is a list with all the parameter values for this model in the correct order and format, but all
# values are 0. Thus, we need to populate these parameters with numbers in order to perform the simulation.
# pars$model.pars[,1] <- c(0.1, 0.1, 0.1, 0.03, 0.03, 0.05, 0.05)
# pars$model.pars[,2] <- c(0.2, 0.2, 0.2, 0.03, 0.03, 0.05, 0.05)
# pars$q.01[1,2] <- pars$q.01[2,1] <- 0.005
# pars$q.0[1,2] <- pars$q.0[2,1] <- 0.005
# pars$q.1[1,2] <- pars$q.1[2,1] <- 0.005
# pars

# Now we can use the parameters with the same function we applied before SimulateGeoHiSSE to generate
# both the data and the phylogeny.
# Here we will set the seed for the simulation, so the outcome of the simulation is always the same. Note that
# you can change the seed or skip this lines to generate a different, random, dataset.
#set.seed(42)
#sim.geohisse <- SimulateGeoHiSSE(pars=pars, hidden.traits = 1, x0 = "01A", max.taxa = 500)
#sim.dat <- data.frame(taxon=sim.geohisse$data[,1], ranges=as.numeric(sim.geohisse$data[,2]))
phy <- read.tree("Supplementary_DataD11_time-calibrated-tree.tre")
plot(phy)
phy$node.labels <- NULL
sim.dat <- read.csv("Supplementary_DataD12_homa_aq_geohisse.csv", stringsAsFactors=FALSE)
View(sim.dat)
#trim the tips you don't need
phy<-drop.tip(phy, "Bothrops_moojeni")
phy<-drop.tip(phy, "Calamophis_ruudelangi")
phy<-drop.tip(phy, "Brachyorrhos_wallacei")
phy<-drop.tip(phy, "Brachyorrhos_albus")
phy<-drop.tip(phy, "Brachyorrhos_raffrayi")
phy<-drop.tip(phy, "Brachyorrhos_gastrotaenius")
plot(phy)

###Setting up the Models
# We will fit a total of four models. Two models with a range-indendent diversification process and two other
# models in which the range have an effect on the diversification rate of the lineages (each with either one or
# two rate classes).

# Note that the function to estimate the parameters of the model is commented out below. Just uncomment
# and run to perform the estimate of the models. Here we will load results from a previous estimate. READ THIS***************************!!!!!!!!!

# Models 1 and 2 below do not include hidden classes. For model 1 (mod1), we are assuming equal rates
# regardless of biogeographic region. This requires a particular set up of the turnover rate with respect to the
# widespread range, which involves removing it from the model. However, if assume.cladogenetic=TRUE, this
# does not mean that we are excluding it from the model. Internally, GeoHiSSE() will recognize this and set the
# speciation rate for, s01, to be equal to the rate of s00 and s11. Removing the turnover rate for the widespread
# range like this is required for any model where the diversification rates are independent of range evolution:

#note, we we have 1 homalopsis unknown out of the 47 (and H. hardwickii is assumed to be freshwater)

samp.frac.fresh <- (20/23) #out of 26 freshwater homalopsids
samp.frac.brack <- (9/17) #out of 17 brackish homalopsids
samp.frac.both  <- (6/6) #out of 3 'freshwater' freshwater/brackish homalopsids

samp.frac <- c(samp.frac.fresh,samp.frac.brack,samp.frac.both)




#####if choosing one samp frac for all states: samp.frac <-nrow(states)/#of homalopsids of rear-fanged group/clade of interest
#These are sampling fractions; you're counting for species that aren't in the tree; important parameter as if you have incomplete sampling of phylogeny and are doing diversification analyses, taxa missing form the tree look like extinctinos or a non-speciation (speciation will be underestimated and extinctino over estimated)
## Get fraction of fresh species in the tree - there are 3645 snake species on Reptile Database as of Dec. 7, 2017 **********must change the # of taxa to how many rear-fanged homalopsids there are
#non_arb_frac<-length(states[states[,2]==1,2])/(samp.frac) #change length(all_arb_species) to total # of fresh homalopsids (nrow(states)/(# of rear-fanged homalopsids) for each state) # number of species coded as non-arboreal (state 1)/total number of snakes minus total # we have identified as arboreal
## Get fraction of primarily brackish species that are in the tree
#arb_frac<-length(states[states[,2]==2,2])/(samp.frac)  # change length(all_arb_species) total # of brackish homalopsids;  number of species in tree coded as primarily arboreal (state 2) divided by total number we identified as primarily arboreal
## Get fraction of non-aquatic species that are in the tree
#semi_arb_frac<-length(states[states[,2]==0,2])/(samp.frac)  # change length(all_arb_species) number of non-aquatic in tree coded as semi arboreal (state 0) divided by total number we identified as semi arboreal
## Set up the sampling fraction
##samp_frac<-c(semi_arb_frac, non_arb_frac, arb_frac) #vector of 3 diff. values; sampling fraction for states 0, 1, and 2 (the order here is important); must have a vector of three here whether you choose the same or different numbers
#samp_frac<-c(samp.frac, samp.frac, samp.frac)

## Model 1 - Dispersal parameters vary only, no range-dependent diversification.
turnover <- c(1,1,0)
eps <- c(1,1)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=0)
trans.rate.mod <- ParEqual(trans.rate, c(1,2))
mod1 <- GeoHiSSE(phy = phy, data = sim.dat, f=samp.frac,
                 turnover=turnover, eps=eps,
                 hidden.states=FALSE, trans.rate=trans.rate.mod,
                 turnover.upper=100, trans.upper=10)

# To conduct a canonical GeoSSE model, where range evolution affects diversification, we add back the turnover
# rate for the widespread range, such that there are three turnover parameters and two extinction fraction
# parameters estimated:
## Model 2. Canonical GeoSSE model, range effect on diversification
turnover <- c(1,2,3)
eps <- c(1,1)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=0)
trans.rate.mod <- ParEqual(trans.rate, c(1,2))
mod2 <- GeoHiSSE(phy = phy, data = sim.dat, f=samp.frac,
                 turnover=turnover, eps=eps,
                 hidden.states=FALSE, trans.rate=trans.rate.mod,
                 turnover.upper=100, trans.upper=10)
View(data)
# Models 3 and 4 below each have 2 hidden states. In this case the models will be more complex. First, we
# will show how to set up a range-independent model of diversification. Remember, as with mod1 above, if
# diversification is independent of range-evolution we must remove the turnover rate for the widespread range.
# Again, internally, GeoSSE() will recognize this (if assume.cladogenetic=TRUE) and simply set s01 to be
# equal the rate for s00 and s11.
## Model 3. GeoHiSSE model with 1 hidden trait, no range-dependent diversification.
## Note below how parameters vary among hidden classes but are the same within each
## hidden class.
turnover <- c(1,1,0,2,2,0)
eps <- c(1,1,1,1)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=1, make.null=TRUE)
trans.rate.mod <- ParEqual(trans.rate, c(1,2))
mod3 <- GeoHiSSE(phy = phy, data = sim.dat, f=samp.frac,
                 turnover=turnover, eps=eps,
                 hidden.states=TRUE, trans.rate=trans.rate.mod,
                 turnover.upper=100, trans.upper=10)

# Finally, if we want to fit a GeoHiSSE model we would do the following:
## Model 4. GeoHiSSE model with 1 hidden trait, range-dependent diversification.
turnover <- c(1,2,3,4,5,6)
eps <- c(1,1,1,1)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=1)
trans.rate.mod <- ParEqual(trans.rate, c(1,2,3,4))  # I added in ,3,4 here -SMH
mod4 <- GeoHiSSE(phy = phy, data = sim.dat, f=samp.frac,
                 turnover=turnover, eps=eps,
                 hidden.states=TRUE, trans.rate=trans.rate.mod,
                 turnover.upper=100, trans.upper=10)


# Make a state-independent version of mod4 - this is model 4a
turnover <- c(1,1,0,2,2,0,3,3,0,4,4,0,5,5,0,6,6,0)
eps <- rep(1, 12)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=5, make.null=TRUE)
trans.rate.mod <- ParEqual(trans.rate, c(1,2))
mod4a <- GeoHiSSE(phy = phy, data = sim.dat, f=samp.frac,
                  turnover=turnover, eps=eps,
                  hidden.states=TRUE, trans.rate=trans.rate.mod,
                  turnover.upper=100, trans.upper=10)


# just wanna do a quick comparison of these so far
# GetAICWeights(list(model1 = mod1, model2 = mod2, model3 = mod3, model4 = mod4, model4a = mod4a), criterion="AIC")

# So far, it looks like model 1 is best fit, followed by model 2





# We will also show how to fit a complementary set of models that remove the cladogenetic effect entirely, such
# that all changes occur along branches (i.e., anagenetic change). This requires the removal of the turnover
# rate for lineages in the widespread range and ensuring that range contraction is distinct from the extinction
# of endemics:
## Model 5. MuSSE-like model with no hidden trait, no cladogenetic effects. simplest state dependent model
turnover <- c(1,2,0)
eps <- c(1,1)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=0, make.null=FALSE, 
                                    separate.extirpation = TRUE)
trans.rate.mod <- ParEqual(trans.rate, c(1,2))                         
trans.rate.mod <- ParEqual(trans.rate.mod, c(2,3))
mod5 <- GeoHiSSE(phy = phy, data = sim.dat, f=samp.frac,
                 turnover=turnover, eps=eps,
                 hidden.states=FALSE, trans.rate=trans.rate.mod,
                 turnover.upper=100, trans.upper=10,
                 assume.cladogenetic = FALSE)


## I just wanna see how this compares to some of the other models
# GetAICWeights(list(model1 = mod1, model2 = mod2, model3 = mod3, model4 = mod4, model4a = mod4a, model5=mod5), criterion="AIC")

# So far, this fits the best

# Model 6 - a model like model 1, but including separate.extirpation=TRUE
turnover <- c(1,1,0)
eps <- c(1,1)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=0, separate.extirpation = TRUE)
trans.rate.mod <- ParEqual(trans.rate, c(1,2))
trans.rate.mod <- ParEqual(trans.rate.mod, c(2,3))
mod6 <- GeoHiSSE(phy = phy, data = sim.dat, f=samp.frac,
                 turnover=turnover, eps=eps,
                 hidden.states=FALSE, trans.rate=trans.rate.mod,
                 turnover.upper=100, trans.upper=10)


# Model 7 - a separate.extirpation version of model 3
turnover <- c(1,1,0,2,2,0)
eps <- c(1,1,1,1)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=1, make.null=TRUE, separate.extirpation = TRUE)
trans.rate.mod <- ParEqual(trans.rate, c(3,4))
trans.rate.mod <- ParEqual(trans.rate.mod, c(1,2))
mod7 <- GeoHiSSE(phy = phy, data = sim.dat, f=samp.frac,
                 turnover=turnover, eps=eps,
                 hidden.states=TRUE, trans.rate=trans.rate.mod,
                 turnover.upper=100, trans.upper=10)


## Model 8 - a separate.extirpation version of model 4
turnover <- c(1,2,0,3,4,0)
eps <- c(1,1,1,1)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=1, separate.extirpation = TRUE)
trans.rate.mod <- ParEqual(trans.rate, c(1,2))
trans.rate.mod <- ParEqual(trans.rate.mod, c(2,3))
trans.rate.mod <- ParEqual(trans.rate.mod, c(3,4))
trans.rate.mod <- ParEqual(trans.rate.mod, c(4,5))
trans.rate.mod <- ParEqual(trans.rate.mod, c(2,4))
trans.rate.mod <- ParEqual(trans.rate.mod, c(1,3))
mod8 <- GeoHiSSE(phy = phy, data = sim.dat, f=samp.frac,
                 turnover=turnover, eps=eps,
                 hidden.states=TRUE, trans.rate=trans.rate.mod,
                 turnover.upper=100, trans.upper=10)


## model 9 - a separate.extirpation=TRUE of model 4a
turnover <- c(1,1,0,2,2,0,3,3,0,4,4,0)
eps <- rep(1, 8)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=3, make.null=TRUE, separate.extirpation = TRUE)
trans.rate.mod <- ParEqual(trans.rate, c(1,2))
trans.rate.mod <- ParEqual(trans.rate.mod, c(2,3))
mod9 <- GeoHiSSE(phy = phy, data = sim.dat, f=samp.frac,
                 turnover=turnover, eps=eps,
                 hidden.states=TRUE, trans.rate=trans.rate.mod,
                 turnover.upper=100, trans.upper=10)





# another quick comparison of these so far
# GetAICWeights(list(model1 = mod1, model2 = mod2, model3 = mod3, model4 = mod4, model4a = mod4a, model5=mod5,
#   model6=mod6, model7 = mod7, model8=mod8, model9=mod9), criterion="AIC")



######### Alright, now let's start making variants of the separate.extirpation=TRUE models that have more transition rates

# Model 10 - same as 5 but with 4 trans rates
turnover <- c(1,2,0)
eps <- c(1,1)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=0, make.null=FALSE, separate.extirpation = TRUE)
mod10 <- GeoHiSSE(phy = phy, data = sim.dat, f=samp.frac,
                  turnover=turnover, eps=eps,
                  hidden.states=FALSE, trans.rate=trans.rate,
                  turnover.upper=100, trans.upper=10,
                  assume.cladogenetic = FALSE)


# Model 11 - same as mod 6 but with 4 trans rates
turnover <- c(1,1,0)
eps <- c(1,1)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=0, separate.extirpation = TRUE)
mod11 <- GeoHiSSE(phy = phy, data = sim.dat, f=samp.frac,
                  turnover=turnover, eps=eps,
                  hidden.states=FALSE, trans.rate=trans.rate,
                  turnover.upper=100, trans.upper=10)




# Model 12 - same as 7 but with asymmetric trans rates
turnover <- c(1,1,0,2,2,0)
eps <- c(1,1,1,1)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=1, make.null=TRUE, separate.extirpation = TRUE)
mod12 <- GeoHiSSE(phy = phy, data = sim.dat, f=samp.frac,
                  turnover=turnover, eps=eps,
                  hidden.states=TRUE, trans.rate=trans.rate,
                  turnover.upper=100, trans.upper=10)



# Model 13 - same as 8 but with asymmetric trans rates
turnover <- c(1,2,0,3,4,0)
eps <- c(1,1,1,1)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=1, separate.extirpation = TRUE)
trans.rate.mod <- ParEqual(trans.rate, c(3,7))
trans.rate.mod <- ParEqual(trans.rate.mod, c(4,7))
trans.rate.mod <- ParEqual(trans.rate.mod, c(1,5))
trans.rate.mod <- ParEqual(trans.rate.mod, c(2,5))
mod13 <- GeoHiSSE(phy = phy, data = sim.dat, f=samp.frac,
                  turnover=turnover, eps=eps,
                  hidden.states=TRUE, trans.rate=trans.rate.mod,
                  turnover.upper=100, trans.upper=10)


# Model 14 - same as 9 but with asymmetric trans rates
turnover <- c(1,1,0,2,2,0,3,3,0,4,4,0)
eps <- rep(1, 8)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=3, make.null=TRUE, separate.extirpation = TRUE)
mod14 <- GeoHiSSE(phy = phy, data = sim.dat, f=samp.frac,
                  turnover=turnover, eps=eps,
                  hidden.states=TRUE, trans.rate=trans.rate,
                  turnover.upper=100, trans.upper=10)




# Let's see how these all shake out
weights <- GetAICWeights(list(model1 = mod1, model2 = mod2, model3 = mod3, model4 = mod4, model4a = mod4a, model5=mod5,
                              model6=mod6, model7 = mod7, model8=mod8, model9=mod9, model10 = mod10, model11 = mod11, model12 = mod12,
                              model13 =mod13, model14 = mod14), criterion="AIC")

sort(weights)


# Looks like models 5 & 6 are still by far the best, cool. Means we don't need to explore any of the even more complex possibilities








#################################################################################################################################################



# An explicit three-state MuSSE/MuHiSSE model can, and probably should, be included in the set of models.
# This can be done by using the MuHiSSE() function. The details for doing so can be found in the Running a
# Multistate HiSSE model vignette.

### Computing Akaike Weights
# Akaike weights are important to evaluate the relative importance of each of the models to explain the variation
# observed in the data. This quantity takes into account pennalties associated to the number of free parametes.
# Models with higher weight show better fit to the data and, as a result, have more weight when performing
# model averaging (see below).

# To compute model weight we can use one of the functions of the package. This will work with both HiSSE
# and GeoHiSSE objects.

load("geohisse_new_vignette.Rsave")
GetAICWeights(list(model1 = mod1, model2 = mod2, model3 = mod3, model4 = mod4, model4a = mod4a, model5 = mod5, model6 = mod6, model7 = mod7, model8 = mod8,
                   model9 = mod9, model10 = mod10, model11 = mod11, model12 = mod12, model13 = mod13, model14 = mod14), criterion="AIC")

## As the number of models in the set grows, naming each model in the set can become hard.
## So one can use a list (created by some automated code) as an input also:
list.geohisse <- list(model1 = mod1, model2 = mod2, model3 = mod3, model4 = mod4, model4a = mod4a, model5 = mod5, model6 = mod6, model7 = mod7, model8 = mod8,
                      model9 = mod9, model10 = mod10, model11 = mod11, model12 = mod12, model13 = mod13, model14 = mod14)
GetAICWeights(list.geohisse, criterion="AIC")
AICWeights <- GetAICWeights(list.geohisse, criterion="AIC")
#write out the output
write.csv(AICWeights, file = "AICWeights.csv")

### Model averaging and plotting.

# Now we can model average the results. Note that this step will reflect the Akaike model weights that we
# computed above.

# For this we need first to perform a marginal reconstruction for each of the models in the set. This will
# reconstruct the hidden states at the nodes of the phylogeny. Then we can use this information to compute
# the model average for the rates.

# These can take a while to run. We will load the results of previous analyses. Uncomment the code below to
# perform the reconstructions.

recon.mod1 <- MarginReconGeoSSE(phy = mod1$phy, data = mod1$data, f = mod1$f,
                                pars = mod1$solution, hidden.states = 1,
                                root.type = mod1$root.type, root.p = mod1$root.p,
                                AIC = mod1$AIC, n.cores = 1)
recon.mod2 <- MarginReconGeoSSE(phy = mod2$phy, data = mod2$data, f = mod2$f,
                                pars = mod2$solution, hidden.states = 1,
                                root.type = mod2$root.type, root.p = mod2$root.p,
                                AIC = mod2$AIC, n.cores = 1)
recon.mod3 <- MarginReconGeoSSE(phy = mod3$phy, data = mod3$data, f = mod3$f,
                                pars = mod3$solution, hidden.states = 2,
                                root.type = mod3$root.type, root.p = mod3$root.p,
                                AIC = mod3$AIC, n.cores = 1)
recon.mod4 <- MarginReconGeoSSE(phy = mod4$phy, data = mod4$data, f = mod4$f,
                                pars = mod4$solution, hidden.states = 2,
                                root.type = mod4$root.type, root.p = mod4$root.p,
                                AIC = mod4$AIC, n.cores = 1)
recon.mod4a <- MarginReconGeoSSE(phy = mod4a$phy, data = mod4a$data, f = mod4a$f,
                                 pars = mod4a$solution, hidden.states = 2,
                                 root.type = mod4a$root.type, root.p = mod4a$root.p,
                                 AIC = mod4a$AIC, n.cores = 1)
recon.mod5 <- MarginReconGeoSSE(phy = mod5$phy, data = mod5$data, f = mod5$f,
                                pars = mod5$solution, hidden.states = 2,
                                root.type = mod5$root.type, root.p = mod5$root.p,
                                AIC = mod5$AIC, n.cores = 1)
recon.mod6 <- MarginReconGeoSSE(phy = mod6$phy, data = mod6$data, f = mod6$f,
                                pars = mod6$solution, hidden.states = 2,
                                root.type = mod6$root.type, root.p = mod6$root.p,
                                AIC = mod6$AIC, n.cores = 1)
recon.mod7 <- MarginReconGeoSSE(phy = mod7$phy, data = mod7$data, f = mod7$f,
                                pars = mod7$solution, hidden.states = 2,
                                root.type = mod7$root.type, root.p = mod7$root.p,
                                AIC = mod7$AIC, n.cores = 1)
recon.mod8 <- MarginReconGeoSSE(phy = mod8$phy, data = mod8$data, f = mod8$f,
                                pars = mod8$solution, hidden.states = 2,
                                root.type = mod8$root.type, root.p = mod8$root.p,
                                AIC = mod8$AIC, n.cores = 1)
recon.mod9 <- MarginReconGeoSSE(phy = mod9$phy, data = mod9$data, f = mod9$f,
                                pars = mod9$solution, hidden.states = 2,
                                root.type = mod9$root.type, root.p = mod9$root.p,
                                AIC = mod9$AIC, n.cores = 1)
recon.mod10 <- MarginReconGeoSSE(phy = mod10$phy, data = mod10$data, f = mod10$f,
                                 pars = mod10$solution, hidden.states = 2,
                                 root.type = mod10$root.type, root.p = mod10$root.p,
                                 AIC = mod10$AIC, n.cores = 1)
recon.mod11 <- MarginReconGeoSSE(phy = mod11$phy, data = mod11$data, f = mod11$f,
                                 pars = mod11$solution, hidden.states = 2,
                                 root.type = mod11$root.type, root.p = mod11$root.p,
                                 AIC = mod11$AIC, n.cores = 1)
recon.mod12 <- MarginReconGeoSSE(phy = mod12$phy, data = mod12$data, f = mod12$f,
                                 pars = mod12$solution, hidden.states = 2,
                                 root.type = mod12$root.type, root.p = mod12$root.p,
                                 AIC = mod12$AIC, n.cores = 1)
recon.mod13 <- MarginReconGeoSSE(phy = mod13$phy, data = mod13$data, f = mod13$f,
                                 pars = mod13$solution, hidden.states = 2,
                                 root.type = mod13$root.type, root.p = mod13$root.p,
                                 AIC = mod13$AIC, n.cores = 1)
recon.mod14 <- MarginReconGeoSSE(phy = mod14$phy, data = mod14$data, f = mod14$f,
                                 pars = mod14$solution, hidden.states = 2,
                                 root.type = mod14$root.type, root.p = mod14$root.p,
                                 AIC = mod14$AIC, n.cores = 1)



## Load previous results:
load("geohisse_recons_new_vignette.Rsave")

# Now that we have the AIC associated with each model and their reconstruction across the nodes of the tree
# we can compute the model average:

recon.models <- list(recon.mod1, recon.mod2, recon.mod3, recon.mod4, recon.mod4a, recon.mod5, recon.mod6, recon.mod7, recon.mod8, recon.mod9,
                     recon.mod10, recon.mod11, recon.mod12, recon.mod13, recon.mod14)
model.ave.rates <- GetModelAveRates(x = recon.models, type = "tips")
#write out the output
write.csv(model.ave.rates, file = "model.ave.rates.csv")

# The result of the reconstrution is a matrix with the parameter estimates for each of the tips species averaged
# over all models. Note that for the GeoSSE model there is no "extinction" parameter associated with
# widespread (01) lineages. Also not that one can change the type of model averaging (between tips, nodes,
# and both) when callin the GetModelAveRates function.

head( model.ave.rates )

# Finally, we can plot the use the resulting data matrix to make a plot of the results.

plot.geohisse.states(x = recon.models, rate.param = "net.div", type = "phylogram",
                     show.tip.label = FALSE, legend = FALSE,  legend.kernel = "hist", state.colors = c("light blue", "royal blue", "black"), edge.width = 6)



sim.dat
plot.geohisse.states(x = recon.models, rate.param = "turnover", type = "phylogram",
                     show.tip.label = FALSE, legend = FALSE, state.colors = c("light blue", "blue", "black"), edge.width = 6)
#?plot.geohisse.states

plot.hisse.states.PIES.SMH(recon_list, rate.param="net.div", show.tip.label=TRUE, edge.width.rate=5, node.pie.size=0.075, 
                           tip.pie.size=0.05, fsize=0.5, state.colors=c("green", "purple"), rate.colors=c("blue", "red"), 
                           legend="tips", legend.cex=3) 
?plot.hisse.states()
###References
# Caetano, D.S., B.C. O'Meara, and J.M. Beaulieu. 2018. Hidden state models improve state-dependent
# diversification approaches, including biogeographic models. Evolution, 72:2308-2324.
