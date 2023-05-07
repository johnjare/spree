# Spree
A tool for automatically generating species trees

## Author
Jared Johnson

## Synopsis

Spree is a tool that helps you classify microbial isolates using whole genome assemblies. It does this by building a neighbor-joining Mash tree with up to five representatives of each species of a specified genus. Representative genomes are automatically downloaded from the NCBI database and are preferentially selected based on multiple quality metrics: database source (RefSeq > GenBank), assembly level (Complete Genome > Chromsome > Scaffold > Contig), unusual genome lengths (± 3 stdevs), and estimated depth of coverage (≥ 30X). The primary results are a phylogenetic tree and the top 10 closest related genomes based on Mash distance. The name Spree is a portmanteau of "species" and "tree".

## Quick Start
```
spree -g Streptococcus -o strep-tree -t 8 strep-1.fasta
```

# Installation
## Docker
```
docker pull johnjare/spree:latest
```
## Source
The scripts below are meant as a template. It is likley you will want to change the path that the executables are moved to. There is no need to 
re-install the packages you already have.
```
% Install required R packages
R -e "install.packages(c('tidyverse', 'ggnewscale', 'phangorn','BiocManager','RColorBrewer','progress'), repos='https://cran.rstudio.com/')"
R -e "BiocManager::install('treeio'); BiocManager::install('ggtree')"

% Install Mash
wget https://github.com/marbl/Mash/releases/download/v2.3/mash-Linux64-v2.3.tar && tar -xvf mash-Linux64-v2.3.tar && rm mash-Linux64-v2.3.tar && mv 
mash-Linux64-v2.3/mash bin/

% Install NCBI Datasets & Dataformat
wget https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets && wget 
https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/dataformat && chmod +x datasets dataformat && mv datasets bin/ && mv dataformat 
bin/

% Install Spree
git clone https://github.com/johnjare/spree.git --branch v1.0 && mv spree/spree bin/
```

# Usage
## Inputs
```
usage: ./spree [-h] -g GENUS [-p PREFIX] [-o OUTPUT] [-t THREADS] [-l LIMIT]
               [--force]
               [input ...]

positional arguments:
  input       input assembly file(s) (.FASTA)

options:
  -h, --help  show this help message and exit
  -g GENUS    NCBI Taxon (e.g., Streptococcus or 1301)
  -p PREFIX   Prefix for run-specific files (e.g., cool_project1) (default
              timestamp)
  -o OUTPUT   output directory (default 'spree')
  -t THREADS  number of threads to use (default 1)
  -l LIMIT    total download size limit (GB) (default 1)
  --force     force download **WARNING** only use this if you are confident
              you have enough storage.
```

## Outputs
### Species-specific outputs
These files can be re-used if the same output directory is supplied. This was designed to save time and hard drive space.
- A tab-separated file containing information about all publicly available genomes for the selected genus in the NCBI database (e.g., `Streptococcus-ncbi.tsv`)
- A directory containing up to five representative genomes of each species of the specified genus (e.g., `Streptococcus/`). These genomes are selected using the tab-separated file described above.
- A tab-separated file containing information about the selected genomes (e.g., `Streptococcus-select-refs.tsv`). This includes all of the information in the original tab-separated file plus some additional metrics used during the selection process.

### Run-specific outputs
These files are specific to the run and will be overwritten if the same prefix is used.
- A newick tree file of the neighbor-joining Mash tree generated from your input assembly/assemblies and the representative assemblies selected from NCBI (e.g., 
`cool_project1-mash-ava.tree`)
- A tree image made from the newick file described above (e.g., `cool_project1-mash-tree.jpg`). Input samples are displayed in bold. References are displayed in light grey and contain information about the species and the genome quality.
- A tab-separated file containing the metadata used to make the make the image file using ggtree (e.g., `cool_project1-metadata.csv`). Use this file to learn more about the quality of the reference assemblies.
- A tab-separated file containing the top 10 closest genomes for each isolate included in the analysis (e.g., `cool_project1-mash-top-10.tsv`)

### Tree Annotations
Each sample in the tree has a colored node and number (displayed in brackets) that correspond to one of the species listed in the legend. The total number of unique assemblies of each species available in NCBI is also displayed in parentheses next to the species names in the lengend. See below for further explanation:
#### Legend Example
<img src="https://github.com/johnjare/spree/blob/main/examples/legend_example.png" width="300">

#### Sample Example
<img src="https://github.com/johnjare/spree/blob/main/examples/node_example.png" width="300">

# Examples
## Blastomyces
### Command
```
spree -g Blastomyces -o blasto -p blast-example -t 8 -l 2 example_genomes/Blastomyces_percursus_GCA_018296075.1.fasta
```
### NJ Mash Tree
In this example, we can see that `Blastomyces_percursus_GCA_018296075.1.fasta` was correctly classified as *Blastomyces_percursus*; however, can also see one of the limitations of relying on species classifications from NCBI. Assembly GCA_003206215.1 was reported as *Blastomyces emzantsi* in NCBI but was likely misclassified and is actually *Blastomyces percursus*. This demonstrates why you should always be cautions when interpreting these trees.

<img src="https://github.com/johnjare/spree/blob/main/examples/blasto-example-mash-tree.jpg" width="1000">

### Top 10 Closest Genomes
|Sample                                      |Reference       |Species                            | Mash Distance| % Matching Kmers|
|:-------------------------------------------|:---------------|:----------------------------------|-------------:|----------------:|
|Blastomyces_percursus_GCA_018296075.1.fasta |GCA_003206225.1 |[5] Blastomyces percursus (n=14)   |     0.0066634|             76.9|
|Blastomyces_percursus_GCA_018296075.1.fasta |GCA_018296065.1 |[5] Blastomyces percursus (n=14)   |     0.0070168|             75.9|
|Blastomyces_percursus_GCA_018296075.1.fasta |GCA_003206215.1 |[2] Blastomyces emzantsi (n=8)     |     0.0083078|             72.4|
|Blastomyces_percursus_GCA_018296075.1.fasta |GCA_003206295.1 |[5] Blastomyces percursus (n=14)   |     0.0085767|             71.7|
|Blastomyces_percursus_GCA_018296075.1.fasta |GCA_001883805.1 |[5] Blastomyces percursus (n=14)   |     0.0091258|             70.3|
|Blastomyces_percursus_GCA_018296075.1.fasta |GCA_003206275.1 |[5] Blastomyces percursus (n=14)   |     0.0092455|             70.0|
|Blastomyces_percursus_GCA_018296075.1.fasta |GCA_000003525.2 |[1] Blastomyces dermatitidis (n=5) |     0.1037400|              6.0|
|Blastomyces_percursus_GCA_018296075.1.fasta |GCA_000151595.1 |[1] Blastomyces dermatitidis (n=5) |     0.1037400|              6.0|
|Blastomyces_percursus_GCA_018296075.1.fasta |GCF_000003525.1 |[1] Blastomyces dermatitidis (n=5) |     0.1037400|              6.0|
|Blastomyces_percursus_GCA_018296075.1.fasta |GCA_000166155.1 |[1] Blastomyces dermatitidis (n=5) |     0.1052640|              5.8|

## Escherichia
### Command
```
spree -g Escherichia -o strep -t 8 -l 2 Streptococcus_thermophilus_ATCC_19258.fasta Streptococcus_pyogenes_ATCC_12344.fasta
```
### NJ Mash Tree
In this example, we can see the result of running several input assemblies at a time. While these assemblies were different species, they still classified correctly in each case. This plot also demonstrates a weakeness. We can see that synthetic *E. coli* were included in the plot. While this does not affect our interpretation, it is a reminder that we are at the mercy of anyone who submits to NCBI.

<img src="https://github.com/johnjare/spree/blob/main/examples/esch-example-mash-tree.jpg" width="1000">

### Top 10 Closest Genomes
|Sample                                  |Reference       |Species                            | Mash Distance| % Matching Kmers|
|:---------------------------------------|:---------------|:----------------------------------|-------------:|----------------:|
|Escherichia_coli_ATCC_11775.fasta       |GCF_000833145.1 |[2] Escherichia coli (n=212702)    |     0.0289561|             37.4|
|Escherichia_coli_ATCC_11775.fasta       |GCF_000988355.1 |[2] Escherichia coli (n=212702)    |     0.0296126|             36.7|
|Escherichia_coli_ATCC_11775.fasta       |GCA_000474035.1 |[7] synthetic Escherichia (n=3)    |     0.0298031|             36.5|
|Escherichia_coli_ATCC_11775.fasta       |GCA_000826905.1 |[7] synthetic Escherichia (n=3)    |     0.0298031|             36.5|
|Escherichia_coli_ATCC_11775.fasta       |GCA_000826925.1 |[7] synthetic Escherichia (n=3)    |     0.0298031|             36.5|
|Escherichia_coli_ATCC_11775.fasta       |GCF_000833635.2 |[2] Escherichia coli (n=212702)    |     0.0298988|             36.4|
|Escherichia_coli_ATCC_11775.fasta       |GCF_000987875.1 |[2] Escherichia coli (n=212702)    |     0.0319897|             34.3|
|Escherichia_coli_ATCC_11775.fasta       |GCF_000967155.2 |[2] Escherichia coli (n=212702)    |     0.0326168|             33.7|
|Escherichia_coli_ATCC_11775.fasta       |GCA_029876145.1 |[5] Escherichia ruysiae (n=14)     |     0.0653942|             14.5|
|Escherichia_coli_ATCC_11775.fasta       |GCA_019839465.1 |[5] Escherichia ruysiae (n=14)     |     0.0659724|             14.3|
|Escherichia_fergusonii_ATCC_35469.fasta |GCF_020097475.1 |[3] Escherichia fergusonii (n=212) |     0.0000716|             99.7|
|Escherichia_fergusonii_ATCC_35469.fasta |GCF_008064895.1 |[3] Escherichia fergusonii (n=212) |     0.0114865|             64.7|
|Escherichia_fergusonii_ATCC_35469.fasta |GCF_008064875.1 |[3] Escherichia fergusonii (n=212) |     0.0117564|             64.1|
|Escherichia_fergusonii_ATCC_35469.fasta |GCF_008064915.1 |[3] Escherichia fergusonii (n=212) |     0.0136002|             60.2|
|Escherichia_fergusonii_ATCC_35469.fasta |GCF_003944565.2 |[3] Escherichia fergusonii (n=212) |     0.0154006|             56.7|
|Escherichia_fergusonii_ATCC_35469.fasta |GCF_000833635.2 |[2] Escherichia coli (n=212702)    |     0.0716227|             12.5|
|Escherichia_fergusonii_ATCC_35469.fasta |GCA_000474035.1 |[7] synthetic Escherichia (n=3)    |     0.0719629|             12.4|
|Escherichia_fergusonii_ATCC_35469.fasta |GCA_000826905.1 |[7] synthetic Escherichia (n=3)    |     0.0719629|             12.4|
|Escherichia_fergusonii_ATCC_35469.fasta |GCA_000826925.1 |[7] synthetic Escherichia (n=3)    |     0.0719629|             12.4|
|Escherichia_fergusonii_ATCC_35469.fasta |GCF_000988355.1 |[2] Escherichia coli (n=212702)    |     0.0719629|             12.4|

# Limitations
- This method requires that you know the genus of the organism prior to use. This is often the case in public health when using MALDI-TOF to classify organisms prior to whole genome sequencing.
- Species classifications are based on what is reported to NCBI, which is sometimes unreliable.
- Garbage in equals garbage out. If your assembly or any one of the reference assemblies from NCBI are low quality (incomplete, contaminated, etc.,), then their placement in the tree will be greatly impacted.

