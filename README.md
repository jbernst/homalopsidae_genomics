# Phylogenomic Investigation of Homalopsidae
Data and scripts for the phylogenomic investigation of Mud Snakes (Homalopsidae) using SqCL v2 probe set. 
Raw sequence files for target capture can be found on the Sequence Read Archive under BioProject ID PRJNA792597.

If you modify the code used here, please be sure to cite myself and the original authors of the software in the code!
To cite this repository, please use: [![DOI](https://zenodo.org/badge/644001397.svg)](https://zenodo.org/badge/latestdoi/644001397)


## SqCL Loci Processing
* [Phyluce-1.7.1](https://phyluce.readthedocs.io/en/latest/installation.html) was used for data processing. Check here for scripts (raw data on SR Database):  
  >[https://github.com/jbernst/homalopsidae_genomics/tree/main/phyluce-1.7.1](https://github.com/jbernst/homalopsidae_genomics/tree/main/phyluce-1.7.1)
* Raw data can be found on the Sequence Reach Archive under BioProject ID PRJNA792597.
  >TO INSERT LINK
* Specimen, Barcode, and Adapter Information can be found here:
  >[https://github.com/jbernst/homalopsidae_genomics/blob/main/AppendixS1_Specimen-Lists.xlsx](https://github.com/jbernst/homalopsidae_genomics/blob/main/AppendixS1_Specimen-Lists.xlsx). 

## Phylogenetic Analyses
* [IQTREE v1](http://www.iqtree.org/) and [RAxML-NG](https://github.com/amkozlov/raxml-ng) were used to reconstruct phylogenies. RAxML-NG was used simply for constraint tree, and is not included here. However, the
  IQTREE analyses can be found here: 
  >[https://github.com/jbernst/homalopsidae_genomics/tree/main/IQTREE](https://github.com/jbernst/homalopsidae_genomics/tree/main/IQTREE).  
* [ASTRAL-III](https://github.com/smirarab/ASTRAL) was used to generate the species tree. The input gene trees and code can be found here:  
  >[https://github.com/jbernst/homalopsidae_genomics/tree/main/ASTRAL-III](https://github.com/jbernst/homalopsidae_genomics/tree/main/ASTRAL-III)
* [treePL](https://github.com/blackrim/treePL) was used time-calibrate the ASTRAL species tree. Time calibration points and code can be found at:
  >[https://github.com/jbernst/homalopsidae_genomics/tree/main/treePL](https://github.com/jbernst/homalopsidae_genomics/tree/main/treePL)  
* Secondary node calibrations for treePL were extracted from [Burbrink et al. (2020)](https://doi.org/10.1093/sysbio/syz062) using custom R-code: 
   [https://github.com/jbernst/homalopsidae_genomics/tree/main/treePL/node-divergence-extraction](https://github.com/jbernst/homalopsidae_genomics/tree/main/treePL/node-divergence-extraction)

## Ancestral Range Estimation
* [BioGeoBEARS](https://github.com/nmatzke/BioGeoBEARS) in R was used to perform ancestral range estimation on the time-calibrated trees from treePL. The [BioGeoBEARS Wiki](http://phylo.wikidot.com/biogeobears)
  is extrememly helpful in running these analyses. The code and input files for BioGeoBEARS for the fresh (genomic) and fresh+formalin (cyt-b) runs can be found here:
  >Fresh Samples - [https://github.com/jbernst/homalopsidae_genomics/tree/main/BioGeoBEARS](https://github.com/jbernst/homalopsidae_genomics/tree/main/BioGeoBEARS)  
  >Formalin+Fresh Samples - [https://github.com/jbernst/homalopsidae_genomics/tree/main/BioGeoBEARS_FreshFormalin_CYTB](https://github.com/jbernst/homalopsidae_genomics/tree/main/BioGeoBEARS_FreshFormalin_CYTB)

## Ancestral State Reconstruction
* [APE](https://cran.r-project.org/web/packages/ape/index.html) was used for ancestral state reconstructions on homalopsids. Associated R code and input files are found here:
  >[https://github.com/jbernst/homalopsidae_genomics/tree/main/APE_Ancestral-State-Reconstructions](https://github.com/jbernst/homalopsidae_genomics/tree/main/APE_Ancestral-State-Reconstructions) 

## Hidden Geographic State Speciation and Extinction (GeoHiSSE)
* [GeoHiSSE](https://doi.org/10.1093/sysbio/syw022) was run using the [GeoHiSSE vignette](https://cran.r-project.org/web/packages/hisse/vignettes/GeoHiSSE-vignette.pdf). 
* The modified code and input time-calibrated trees are found at:
  >[https://github.com/jbernst/homalopsidae_genomics/tree/main/GeoHiSSE](https://github.com/jbernst/homalopsidae_genomics/tree/main/GeoHiSSE)

#TO EDIT BELOW

## Source Publication

[Link to Paper]()

DOI to paper: [https://ssbbulletin.org/index.php/bssb/article/view/9393/7948](https://ssbbulletin.org/index.php/bssb/article/view/9393/7948)

![Image](https://live.staticflickr.com/65535/53121066876_e0cb8194fd_h.jpg)
