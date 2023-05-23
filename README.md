# Phylogenomic Investigation of Homalopsidae
Data and scripts for the phylogenomic investigation of Mud Snakes (Homalopsidae) using SqCL v2 probe set. 
Raw sequence files for target capture can be found on the Sequence Read Archive under BioProject ID PRJNA792597.




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

## Ancestral Range Estimation
* [BioGeoBEARS](https://github.com/nmatzke/BioGeoBEARS) in R was used to perform ancestral range estimation on the time-calibrated trees from treePL. The [BioGeoBEARS Wiki](http://phylo.wikidot.com/biogeobears)
  is extrememly helpful in running these analyses. The code and input files for BioGeoBEARS can be found here:
  >

## Ancestral State Reconstruction
* [APE](https://cran.r-project.org/web/packages/ape/index.html) was used for ancestral state reconstructions on homalopsids. Associated R code and input files are found here:
  >[https://github.com/jbernst/homalopsidae_genomics/tree/main/APE_Ancestral-State-Reconstructions](https://github.com/jbernst/homalopsidae_genomics/tree/main/APE_Ancestral-State-Reconstructions) 

## Hidden Geographic State Speciation and Extinction (GeoHiSSE)
* [GeoHiSSE](https://doi.org/10.1093/sysbio/syw022) was run using the [GeoHiSSE vignette](https://cran.r-project.org/web/packages/hisse/vignettes/GeoHiSSE-vignette.pdf). 
* The modified code and input time-calibrated trees are found at:
  >[https://github.com/jbernst/homalopsidae_genomics/tree/main/GeoHiSSE](https://github.com/jbernst/homalopsidae_genomics/tree/main/GeoHiSSE)

#TO EDIT BELOW

## Source Publication
Bernstein, J. M., Murphy, J. C., Voris, H. K., Brown, R. M., Ruane, S. 2021. Phylogenetics of Mud Snakes (Squamata: Serpentes: 
Homalopsidae): A Paradox of Both Undescribed Diversity and Taxonomic Inflation. Molecular Phylogenetics and Evolution 160: 107109.

[Link to Paper](https://static1.squarespace.com/static/633a1ad2337f6700f6fcf3de/t/6344ec2cba833d7984fec29e/1665461295127/Bernstein-et-al_Homalopsidae_MPE2021.pdf)

DOI to paper: [https://doi.org/10.1016/j.ympev.2021.107109](https://doi.org/10.1016/j.ympev.2021.107109)

![Image](https://ars.els-cdn.com/content/image/1-s2.0-S1055790321000427-ga1_lrg.jpg)
