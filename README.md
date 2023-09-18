## Quick start ##
```
git clone https://github.com/zhangrengang/orthoindex.git
cd orthoindex

# install
conda env create -f OrthoIndex.yaml
conda activate OrthoIndex
python3 setup.py install

# test
cd example_data/
sh example.sh

# example.sh:
# dot plots
# A
soi dotplot -s Populus_trichocarpa-Salix_dunnii.collinearity.gz \
        -g Populus_trichocarpa-Salix_dunnii.gff.gz -c Populus_trichocarpa-Salix_dunnii.ctl  \
        --kaks Populus_trichocarpa-Salix_dunnii.collinearity.ks.gz \
        --xlabel '$Populus\ trichocarpa$' --ylabel '$Salix\ dunnii$' \
        --ks-hist --max-ks 1.5 -o Populus_trichocarpa-Salix_dunnii     \
        --plot-ploidy --gene-axis --number-plots
# B
soi dotplot -s Populus_trichocarpa-Salix_dunnii.orthologs.gz    \
        -g Populus_trichocarpa-Salix_dunnii.gff.gz -c Populus_trichocarpa-Salix_dunnii.ctl  \
        --kaks Populus_trichocarpa-Salix_dunnii.collinearity.ks.gz \
        --xlabel '$Populus\ trichocarpa$' --ylabel '$Salix\ dunnii$'  \
        --ks-hist --max-ks 1.5 -o Populus_trichocarpa-Salix_dunnii.o   \
        --plot-ploidy --gene-axis --number-plots  \
        --homology      # homology input
# C
soi dotplot -s Populus_trichocarpa-Salix_dunnii.collinearity.gz \
        -g Populus_trichocarpa-Salix_dunnii.gff.gz -c Populus_trichocarpa-Salix_dunnii.ctl  \
        --xlabel '$Populus\ trichocarpa$' --ylabel '$Salix\ dunnii$' \
        --ks-hist --max-ks 1 -o Populus_trichocarpa-Salix_dunnii.io    \
        --plot-ploidy --gene-axis --number-plots \
        --ofdir OrthoFinder/OrthoFinder/Results_*/ --of-color   # coloring by Orthology Index
# D
soi dotplot -s Populus_trichocarpa-Salix_dunnii.collinearity.gz \
        -g Populus_trichocarpa-Salix_dunnii.gff.gz -c Populus_trichocarpa-Salix_dunnii.ctl  \
        --kaks Populus_trichocarpa-Salix_dunnii.collinearity.ks.gz \
        --xlabel '$Populus\ trichocarpa$' --ylabel '$Salix\ dunnii$' \
        --ks-hist --max-ks 1.5 -o Populus_trichocarpa-Salix_dunnii.io  \
        --plot-ploidy --gene-axis --number-plots \
        --ofdir OrthoFinder/OrthoFinder/Results_*/ --of-ratio 0.6       # filtering by Orthology Index


# filter orthologous synteny
soi filter -s Populus_trichocarpa-Salix_dunnii.collinearity.gz -o OrthoFinder/OrthoFinder/Results_*/ \
        -c 0.6 > Populus_trichocarpa-Salix_dunnii.collinearity.ortho
# or (alter input format)
soi filter -s Populus_trichocarpa-Salix_dunnii.collinearity.gz -o Populus_trichocarpa-Salix_dunnii.orthologs.gz \
        -c 0.6 > Populus_trichocarpa-Salix_dunnii.collinearity.ortho

```
### Example output dot plots ###
![dotplots](example_data/mege_4dot.png)

## Table of Contents

   * [Introduction](#introduction)
   * [Subcommands](#Subcommands)
      - [filter](#filter)
      - [cluster](#cluster)
      - [outgroup](#outgroup)
      - [phylo](#phylo)
      - [dotplot](#dotplot)
   * [Phylogenomics pipeline](#Phylogenomics-pipeline)
   * [Singularity/Apptainer](#Singularity/Apptainer)

## Introduction ##
Orthology Index (OrthoIndex or OI) incorporates algorithmic advances of two methods (orthology inference and synteny detection), to determine the orthology of a syntenic block. 
It is straightforward, representing the proportion of orthologous gene pairs within a syntenic block. 

## Subcommands ##
```
$ soi -h
usage: soi [-h] {dotplot,filter,cluster,outgroup,phylo,stats} ...

Play with Orthology Index

positional arguments:
  {dotplot,filter,cluster,outgroup,phylo,stats}
                        sub-command help
    dotplot             Generate colored dot plots
    filter              Filter synteny with Orthology Index (standard output)
    cluster             Cluster syntenic orthogroups (SOGs)
    outgroup            Add outgroups for SOGs from synteny
    phylo               Build gene trees from SOGs
    stats               Make statistics of SOGs for phylogeny

optional arguments:
  -h, --help            show this help message and exit
```
#### `filter` ####
The subcommand `filter` filters orthologous blocks with a default minimum index of 0.6:
```
$ soi filter -h
usage: soi filter [-h] -s [FILE [FILE ...]] -o [FOLDER/FILE [FOLDER/FILE ...]] [-c FLOAT] [-upper FLOAT] [-n INT]

optional arguments:
  -h, --help            show this help message and exit
  -s [FILE [FILE ...]], -synteny [FILE [FILE ...]]
                        Collinearity files output from MCscanX, WGDI, or MCscan/JCVI. [required]
  -o [FOLDER/FILE [FOLDER/FILE ...]], -orthology [FOLDER/FILE [FOLDER/FILE ...]]
                        Orthologues output from OrthoFinder (folder), or OrthoMCL (file). [required]
  -c FLOAT, -cutoff FLOAT
                        Cutoff (lower limit) of Orthology Index [default=0.6]
  -upper FLOAT          Upper limit of Orthology Index [default=1]
  -n INT, -min_n INT    Minimum gene number in a block [default=0]
```
Usage examples:
```
# from outputs of WGDI and OrthoFinder
soi filter -s wgdi/*.collinearity -o OrthoFinder/*/Result*/ > collinearity.ortho

# from outputs of MCscanX and OrthoMCL
soi filter -s mcscanx/*.collinearity -o orthologs.txt > collinearity.ortho

# from a list file and increase the cutoff
ls wgdi/*.collinearity > collinearity.list
soi filter -s collinearity.list -o OrthoFinder/*/Result*/ -c 0.7 > collinearity.ortho

# filter a paralogous peak
soi filter -s wgdi/*.collinearity -o OrthoFinder/*/Result*/ -c 0.05 -upper 0.4 > collinearity.para
```
#### `cluster` ####
The subcommand ‘cluster’ groups orthologous syntenic genes into syntenic orthogroups (SOGs), through constructing an orthologous syntenic graph 
and applying the Markov Cluster (MCL) algorithm to perform graph clustering and break weak links. 
```
$ soi cluster -h
usage: soi cluster [-h] -s [FILE [FILE ...]] [-o [FOLDER/FILE [FOLDER/FILE ...]]] [-I FLOAT] [-outgroup [TAXON/FILE [TAXON/FILE ...]]] [-ingroup [TAXON/FILE [TAXON/FILE ...]]]
                   [-prefix OUTPRE]

optional arguments:
  -h, --help            show this help message and exit
  -s [FILE [FILE ...]], -synteny [FILE [FILE ...]]
                        Collinearity files from `filter` sub-command. [required]
  -o [FOLDER/FILE [FOLDER/FILE ...]], -orthology [FOLDER/FILE [FOLDER/FILE ...]]
                        Orthologues output from OrthoFinder (folder), or OrthoMCL (file). This will use Orthology Index as weight for MCL [default=None]
  -I FLOAT, -inflation FLOAT
                        Inflation for MCL (varying this parameter affects granularity) [default=1.5]
  -outgroup [TAXON/FILE [TAXON/FILE ...]]
                        Outgroups to exclude from orthogroups (prior to `-ingroup`) [default=None]
  -ingroup [TAXON/FILE [TAXON/FILE ...]]
                        Ingroups that are only included [default=None]
  -prefix OUTPRE        Output prefix [default=cluster]
```
Usage examples:
```
# all species to include
soi cluster -s collinearity.ortho -prefix cluster

# exclude outgroup species that do not share the INGROUP-specific WGD event
soi cluster -s collinearity.ortho -outgroup XXX YYY
```
The defualt output file is `cluster.mcl`, with the orthogroup format of legacy OrthoMCL.

#### `outgroup` ####
The subcommand ‘outgroup’ retrieves syntenic orthologs from outgroups that lack WGDs shared with ingroups. 
```
$ soi outgroup -h
usage: soi outgroup [-h] -s [FILE [FILE ...]] -og FILE -outgroup [TAXON [TAXON ...]] [-cutoff FLOAT]

optional arguments:
  -h, --help            show this help message and exit
  -s [FILE [FILE ...]], -synteny [FILE [FILE ...]]
                        Collinearity files from `filter` sub-command. [required]
  -og FILE, -orthogroup FILE
                        Orthogroups output from `cluster` sub-command. [required]
  -outgroup [TAXON [TAXON ...]]
                        Outgroups to include to orthogroups [required]
  -cutoff FLOAT         Cutoff (lower limit) of links to outgroup genes [default=0.2]
```
Usage examples:
```
# If outgroups are excluded in the last `cluster` step:
soi outgroup -s collinearity.ortho -og cluster.mcl -outgroup XXX YYY > cluster.mcl.plus
```

#### `phylo` ####
The subcommand ‘phylo’ reconstructs multi-copy or single-copy gene trees, 
by aligning protein sequences with MAFFT v7.481 (Standley and Katoh 2013), 
converting protein alignment to codon alignment with PAL2NAL v14 (Suyama et al. 2006), 
trimming alignments with trimAl v1.2 (Capella-Gutierrez et al. 2009) (parameter: -automated1) 
and reconstructing maximum-likelihood trees with IQ-TREE v2.2.0.3 (Minh et al. 2020). 
```
$ soi phylo -h
usage: soi phylo [-h] -og FILE -pep FILE [-cds FILE] [-both] [-fmt STR] [-root [TAXON [TAXON ...]]] [-pre STR] [-mm FLOAT] [-mc INT] [-sc] [-ss FILE] [-concat] [-p INT] [-tmp FOLDER]
                 [-clean]

optional arguments:
  -h, --help            show this help message and exit
  -og FILE, -orthogroup FILE
                        Orthogroups output from `cluster` or `outgroup` sub-commands. [required]
  -pep FILE             Protein fasta file. [required]
  -cds FILE             CDS fasta file. [default=None]
  -both                 To use both CDS and PEP to build gene trees. [default: only CDS when `-cds` is true]
  -fmt STR              Format of `-orthogroup` input. [default=orthomcl]
  -root [TAXON [TAXON ...]], -outgroup [TAXON [TAXON ...]]
                        Outgroups to root gene trees [default=None]
  -pre STR, -prefix STR
                        Output prefix. [default=sog]
  -mm FLOAT, -max_missing FLOAT
                        To allow maximum ratio of missing species. [default=0.4]
  -mc INT, -max_copies INT
                        To limit a common maximum copy number for every species. [default=6]
  -sc, -singlecopy      Only retrieve singlecopy genes (=`-max_copies 1`). [default=None]
  -ss FILE, -spsd FILE  To limit a specific copy number for each species (format: 'TAXON<tab>NUMBER'). [default=None]
  -concat               To concatenate alignments for tools such as IQTREE (valid when `-singlecopy` is true). [default=None]
  -p INT, -ncpu INT     Number of processors. [default=20]
  -tmp FOLDER, -tmpdir FOLDER
                        Temporary folder. [default=./tmp/]
  -clean                Cleanup temporary folder. [default=None]
```
Usage examples:
```
# output multi-copy gene trees of both protein and CDS(-both); rooted with grape (-root)
soi phylo -og cluster.mcl.plus -pep pep.faa -cds cds.fa -both -root Vitis_vinifera -pre mc-sog -p 80

# output single-copy gene trees (-sc) and concatenated alignments (-concat) of both protein and CDS (-both); rooted with grape (-root)
soi phylo -og cluster.mcl.plus -pep pep.faa -cds cds.fa -both -root Vitis_vinifera -pre sc-sog -sc -concat -p 80

# output multi-copy gene trees of protein, allowing up to 50% taxa missing
soi phylo -og cluster.mcl.plus -pep pep.faa -mm 0.5
```

#### `dotplot` ####
The subcommand `dotplot` enables visualization and evaluation of synteny, 
with colored by the Orthology Index or Ks values.

```
$ soi dotplot -h
usage: soi dotplot [-h] -s FILE [FILE ...] -g FILE -c FILE [-o STR] [--format FORMAT] [--homology] [--cluster] [--diagonal] [--gene-axis] [--number-plots] [--min-block INT]
                   [--min-same-block INT] [--xlabel XLABEL] [--ylabel YLABEL] [--figsize NUM] [--fontsize NUM] [--dotsize NUM] [--ofdir FOLDER/FILE [FOLDER/FILE ...]] [--of-ratio FLOAT]
                   [--of-color] [--kaks FILE] [--ks-hist] [--max-ks Ks] [--ks-cmap Ks [Ks ...]] [--ks-step Ks] [--use-median] [--method STR] [--lower-ks Ks] [--upper-ks Ks]
                   [--plot-ploidy] [--window_size INT] [--window_step INT] [--min_block INT] [--max_distance INT] [--max_ploidy INT] [--min_overlap FLOAT] [--color COLOR]
                   [--edgecolor COLOR]

optional arguments:
  -h, --help            show this help message and exit
  -s FILE [FILE ...]    syntenic block file (*.collinearity, output of MCSCANX/WGDI)
  -g FILE               gene annotation gff file (*.gff, one of MCSCANX/WGDI input)
  -c FILE               chromosomes config file (*.ctl, same format as MCSCANX dotplotter)
  -o STR                output file prefix. [default: the same as `-c`]
  --format FORMAT       output figure format [default=['pdf', 'png']]
  --homology            `-s` is in homology format (gene1<tab>gene2). [default=False]
  --cluster             cluster chromosomes. [default=False]
  --diagonal            try to put blocks onto the diagonal. [default=False]
  --gene-axis           use gene as axis instead of base pair. [default=False]
  --number-plots        number subplots with (a-d). [default=False]
  --min-block INT       min gene number in a block. [default=None]
  --min-same-block INT  min gene number in a block on the same chromosome. [default=25]

Art settings:
  art settings for plots

  --xlabel XLABEL       x label for dot plot. [default=None]
  --ylabel YLABEL       y label for dot plot. [default=None]
  --figsize NUM         figure size [default=18]
  --fontsize NUM        font size [default=10]
  --dotsize NUM         dot size [default=0.8]

Orthology Index filter/color:
  filtering or coloring blocks by Orthology Index (prior to Ks color)

  --ofdir FOLDER/FILE [FOLDER/FILE ...]
                        OrthoFinder output folder/ OrthoMCL output pair file. [default=None]
  --of-ratio FLOAT      Orthology Index cutoff [default=0]
  --of-color            coloring dots by Orthology Index [default=None]

Ks plot:
  options to plot with Ks

  --kaks FILE           kaks output from KaKs_Calculator/WGDI. [default=None]
  --ks-hist             plot histogram or not [default=None]
  --max-ks Ks           max Ks (x limit) [default=1]
  --ks-cmap Ks [Ks ...]
                        color map for Ks. [default=None]
  --ks-step Ks          Ks step of histogram [default=0.02]
  --use-median          use median Ks for a block. [default=False]
  --method STR          Ks calculation method [default=NG86]
  --lower-ks Ks         lower limit of median Ks. [default=None]
  --upper-ks Ks         upper limit of median Ks. [default=None]

ploidy plot:
  options to plot relative ploidy (synteny depth)

  --plot-ploidy         plot relative ploidy. [default=False]
  --window_size INT     window_size. [default=50]
  --window_step INT     window_step. [default=10]
  --min_block INT       min genes for a block. [default=10]
  --max_distance INT    max distance. [default=20]
  --max_ploidy INT      x limit. [default=10]
  --min_overlap FLOAT   min overlap. [default=0.4]
  --color COLOR         bar fill color. [default=None]
  --edgecolor COLOR     bar edge color. [default=None]
```

Usage examples: see [Quick Start](#Quick-Start).

### Phylogenomics pipeline ###

See [evolution_example](https://github.com/zhangrengang/evolution_example/) for a pipeline of phylogenomics analyses based on Orthology Index.
