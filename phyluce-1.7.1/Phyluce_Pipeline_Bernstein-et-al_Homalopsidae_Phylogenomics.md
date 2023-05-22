# Phyluce -1.7.1 scripts for processing the SqCL V2 Probe Set from Singhal et al. 2017)                
Code for Phyluce modified from https://phyluce.readthedocs.io/en/latest/index.html (Owner: Brant Faircloth)  
Modification Date: August, 2022  
Dependencies: Anaconda/Miniconda; Phyluce-1.7.1  
#################################################################################################  

I recommend checking the Phyluce ReadtheDocs page for information on configuration files and their formats.
Before we start, it is a good idea to check how many reads your files have and the file sizes. This can be done using 
the sequecning file that was sent to you from your sequencing provider. I found that anything larger than 10-14 million 
reads would likely time out on the HPC cluster (>72 hour run time) during the Assembly step. 

The structure of the workspace I will be working in will be:  
PROJECT_DIRECTORY  
&nbsp;&nbsp;&nbsp;&nbsp;|__RawData/  
&nbsp;&nbsp;&nbsp;&nbsp;|__Assembly_SPADES/

## Subsetting Samples (Optional)

In your RawData directory, you can use seqtk to subset your R1 and R2 files to have less reads (but, less data).
Here, we are obtaining 3.5 million reads from EACH read file (a total of 7000000) in a directory with your R1 and R2 files.

```
READS=3500000
for dir in /path/to/your/clearn/data/dir/from/illumiprocesser/*;
do
        RAND=$RANDOM;
        echo $RAND;
        for file in $dir/split-adapter-quality-trimmed/*-READ[1-2]*;
        do
                echo $file;
                seqtk sample -s $RAND $file $READS | gzip > $file:t:r:r.SUBSAMPLE.fastq.gz
        done;
done
```

# Count the Read Data

Here, we can check how many reads we have per read per sample

```
for i in *_R1_*.fastq.gz; do echo $i; gunzip -c $i | wc -l | awk '{print $1/4}'; done
```

# Clean the Read Data - Trimming with Illumiprocessor

Next, we must remove the adapter sequences and barcodes from the untrimmed data.  Different Illumina instruments use different workflows. 
The Forward Strand Workflow requires the i5 in the forward orientation and the Reverse Complement Workflow requires the reverse complemented i5. 
The NovaSeq v1 chemistry (pre-November 2020) used the Forward Strand Workflow. The v1.5 chemistry was released in November of 2020, and so the 
NovaSeq was switched to the Reverse Complement Workflow. Some sequencing facilities still require that the indexes are provided in the Forward 
Strand orientation because they do a screen on the MiSeq first, which uses the Forward Strand Workflow. So, the indexes that are provided are in 
the Forward Strand Workflow orientation. [Link to workflows here](https://support-docs.illumina.com/SHARE/IndexedSeq/indexed-sequencing.pdf). 

Go to the directory containing your config files and data

```
cd RawData
```


Run illumiprocessor. Remember to change the config file name and number of cores available on your computer.
You will need your illumiprocessor configuration file for this step.

```
illumiprocessor \
    --input raw-fastq/ \
    --output clean-fastq \
    --config illumiprocessor_homalopsidae.conf \
    --cores 16 \
	--r1-pattern _L005_R1_001 \
	--r2-pattern _L005_R2_001
```	

# Quality Control

We can find out how our trimming improved the quality of our reads (though, FASTQC will be a better program to use).

```
# move to the directory holding our cleaned reads
cd clean-fastq/

# run this script against all directories of reads

for i in *;
do
    phyluce_assembly_get_fastq_lengths --input $i/split-adapter-quality-trimmed/ --csv;
done
```

# Data Assembly with SPADES

Make to already have a configuration file for the SPADES assembly step. Here, we will compile our raw data into contigs.
If you are running this on an HPC cluster, you will likely need to do this in batches of 3-5 samples. A small trick is to
take a sample that runs successfully on its own, and then take all other samples and divide their read counts by the test sample'sample
number of reads. Then, multiple this number by how many minutes it took to assemble the test sample. 

So, for example:
Test Sample = 10,000,000 reads. Assembly time = 250 minutes.
Sample 1 = 9,000,000 reads
Sample 2 = 10,000,000 reads
Sample 3 = 11,000,000 reads

If you wanted to see how long it would take to assemble samples 1, 2, and 3 together (and not exceed the 72-hour limit run time found
on many HPC clusters), you would do:
	(9,000,000 + 10,000,000 + 11,000,000) / 10,000,000 = 30,000,000
	(30,000,000 / 10,000,000) \* 250 = 750 minutes (~12.5 hours).

This is a very rough calculation, but could be helpful for assembling your data. I like to copy my clean-fastq directory to the
AssemblySpades directory, so my Assembly steps are separate from everything else (many files and directories will be generated).

```
cp clean-fastq/ ../AssemblySpades
cd ../Assembly_Spades
```

Run the assembly with SPADES. Change cores and memory to reflect your HPC. The below code chunk was run for several batches.
I call all of my output directories 'spades-assembles_{Sample_IDs}.'
```
phyluce_assembly_assemblo_spades \
    --conf HM28_assembly.conf \
    --output spades-assemblies_HM8-11 \
    --cores 28 \
    --memory 182
```

You can copy all of the contigs into a new directory in the parent directory so they are easy to access.

```
mkdir ../homalopsidae_contigs
```

```
for i in * \
do \
cp $i/contigs/*.fasta ../homalopsidae_contigs \
done 
```

Now the PROJECT_DIRECTORY will look like this:  
PROJECT_DIRECTORY  
&nbsp;&nbsp;&nbsp;&nbsp;|__RawData/  
&nbsp;&nbsp;&nbsp;&nbsp;|__Assembly_SPADES/  
&nbsp;&nbsp;&nbsp;&nbsp;|__homalopsidae_contigs/

Let's perform a quality check (QC) on the assemblies.

```
cd ../

for i in spades-assemblies/homalopsidae_contigs/*.fasta;
do
    phyluce_assembly_get_fasta_lengths --input $i --csv;
done
```

The output will show, in order:  
samples,contigs,total bp,mean length,95 CI length,min length,max length,median legnth,contigs >1kb

# Finding the Targeted SqCL Loci

Because we are not using the Tetrapod 5k probe set, we must use the SqCL V2 probe set fast (modifeid to remove 'GN' from gene names,
which was causing downstream syntax issues). The origine probe set is found on [Sonal Singhal's Github](https://github.com/singhal/SqCL).
We will use some regex to capture the appropriate names and loci.

```
phyluce_assembly_match_contigs_to_probes \
  --regex '(\w+-\w+)\_\w+' \
  --contigs /scratch/jmb689/PROCESSED_PHYLUCE_SQCL_homalopsidae/contigs \
  --probes SqCL_mod_noGN.fasta \
  --output SqCL-search-results
```

Note: Subsampled reads will have less loci compared to ones in which all the reads were used!

# Extracting SqCL Loci

At this step, you should be in the PROJECT_DIRECTORY, where we have a taxon set configuration file.

Create an output directory for the taxon set

```
mkdir -p taxon-sets/all
```

Create the initial list of loci for subsequent extraction.

```
phyluce_assembly_get_match_counts \
    --locus-db SqCL-search-results/probe.matches.sqlite \
    --taxon-list-config taxon-set.conf \
    --taxon-group 'all' \
    --incomplete-matrix \
    --output taxon-sets/all/all-taxa-incomplete.conf
```

Now, extract FASTA data that correspond to the loci in the all-taxa-incomplete.conf that was created.

```
cd taxon-sets/all/

mkdir log

phyluce_assembly_get_fastas_from_match_counts \
    --contigs ../../contigs \
    --locus-db ../../SqCL-search-results/probe.matches.sqlite \
    --match-count-output fresh-taxa-incomplete.conf \
    --output fresh-taxa-incomplete.fasta \
    --incomplete-matrix fresh-taxa-incomplete.incomplete \
    --log-path log
```

# Explode the Monolithic FASTA file

This step serves to obtain individual statistics on the UCE assemblies by creating files of each locus AND/OR
for each samples. I do both of these as different analyses will require different input data (if you want individual
gene trees, say, for species tree analysis in ASTRAL-III, get individual loci; for concatenated, get files by sample).

Do it by taxon first.

```
phyluce_assembly_explode_get_fastas_file \
    --input all-taxa-incomplete3.fasta \
    --output exploded-fastas-by-taxon \
    --by-taxon
```

And now explode by locus. The --by-locus flag is the default, so just delete --by-taxon to do this by locus. 
Don't forget to change the output directory name!

```
phyluce_assembly_explode_get_fastas_file \
    --input all-taxa-incomplete3.fasta \
    --output exploded-fastas-by-locus \
```

Get summary stats on the FASTAS

```
for i in exploded-fastas-by-taxon/*.fasta;
do
    phyluce_assembly_get_fasta_lengths --input $i --csv;
done 
```

# Aligning the SqCL Loci

First we need to trim the alignments. The Phyluce docs mention Edge-Trimming is probably enough for taxa with divergences <30-50 mya,
but since we are unsure of the divergences here (and homalopsids and outgroups may push back further than 50 mya, we do internal trimming.

Output FASTA formatted alignments.
```
phyluce_align_seqcap_align \
    --input all-taxa-incomplete3.fasta \
    --output mafft-nexus-internal-trimmed \
    --taxa 169 \
    --aligner mafft \
    --cores 25 \
    --incomplete-matrix \
    --output-format fasta \
    --no-trim \
    --log-path log
```

Then use Gblocks to do the internal trimming.

```
phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed \
    --alignments mafft-nexus-internal-trimmed \
    --output mafft-nexus-internal-trimmed-gblocks \
    --cores 12 \
    --log log
```

Obtain summary stats based on the trimmed alignments.

```
phyluce_align_get_align_summary_data \
    --alignments mafft-nexus-internal-trimmed-gblocks \
    --cores 28 \
    --log-path log
```

# Alignment Cleaning

Next, we will remove the locus names from each alignment, which doesn't look good and contains characters that can cause
issues in subsequent analyses. 

It was this step that the 'GN' in gene names for the nuclear protein coding genes in the SqCL probeset were not being removed. So I took
care of that in the probe set FASTA file.
```
phyluce_align_remove_locus_name_from_files \
    --alignments mafft-nexus-internal-trimmed-gblocks \
    --output mafft-nexus-internal-trimmed-gblocks-clean \
    --cores 25 \
    --log-path log
```

# Final Data Matrices

Create 75% and 95% data matrices. These are called 'completeness' matrices, where 75% completeness 
means that in a study of 100 taxa (total), all alignments will contain at least 75 of these 100 taxa. 
A 95% matrix means that in a study of 100 taxa, all alignments will contain 95 of these 100 taxa.

The 75p matrix:
```
phyluce_align_get_only_loci_with_min_taxa \
    --alignments mafft-nexus-internal-trimmed-gblocks-clean \
    --taxa 169 \
    --percent 0.75 \
    --output mafft-nexus-internal-trimmed-gblocks-clean-75p \
    --cores 15 \
    --log-path log
```

The 95p matrix:
```
phyluce_align_get_only_loci_with_min_taxa \
    --alignments mafft-nexus-internal-trimmed-gblocks-clean \
    --taxa 169 \
    --percent 0.95 \
    --output mafft-nexus-internal-trimmed-gblocks-clean-75p \
    --cores 15 \
    --log-path log
```

# Data Preparation for Downstream Analysis (Concatenation)

We can concatenate all of the alignment for phylogenetic analyses such as RAxML or IQTREE. 
Below is the script just for the 75p matrix.

```
phyluce_align_concatenate_alignments \
    --alignments mafft-nexus-internal-trimmed-gblocks-clean-75p \
    --output mafft-nexus-internal-trimmed-gblocks-clean-75p-raxml \
    --phylip \
    --log-path log
```


