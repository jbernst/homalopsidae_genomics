#The below scripts are not to be run consecutively. These are just a few instances of running IQTREE and preparing for ASTRAL-III analysis.

#IQTREE script for running IQTREE v1 on the concatenated alignments
iqtree -s [input.phy] -bb 1000 -alrt 1000 -m GTR

#Loop for iterating through all alignments given individual gene alignment files in a directory
for i in *.phylip \
do \
iqtree -s $i -bb 1000 -alrt 1000 \
done

# Concatenate all Newick trees for treePL 
cat *.treefile >> ALL-gene=trees.contree
