# Spree
A tool for automatically generating species trees

## Author
Jared Johnson

## Synopsis

Spree is a tool that helps you classify microbial isolates using whole genome assemblies. It does this by building a UPGMA Mash tree with the target isolate(s) and up to up to five representatives of each species of a specified genus. Representative genomes are automatically downloaded from the NCBI database and are preferentially selected based on several quality metrics: database source (RefSeq > GenBank), assembly level (complete > chromsome > scaffold > contig), genome length, and depth of coverage. The primary results are a phylogenetic tree and a list of the top 10 closest related genomes based on Mash distance.

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
These files can be re-used when running a new isolate or batch of isolates of the same species in the same output directory. This is designed to save 
time and hard-drive space.
- A tab-separated file containing all genomes of a specified genus in NCBI (e.g., `${outdir}/Streptococcus-ncbi.tsv`)
- A directory containing up to five representative genomes of each species of a specified genus available in NCBI (e.g., `${outdir}/Streptococcus/`)
- A tab-separated file containing metadata for all of the genomes located in the representative genomes directory (e.g., 
`${outdir}/Streptococcus-select-refs.tsv`)

### Run-specific outputs
The files are specific to the run and will be overwritten if the same prefix is used.
- A newick tree file of the UPGMA Mash tree generated from your input assemblies and the representative assemblies from NCBI (e.g., 
`cool_project1-mash-ava.tree`)
- A tree image made from the newick file described above (e.g., `cool_project1-mash-tree.jpg`)
- A comma-separated file containing the metadata used to make the make the image file using ggtree (e.g., `cool_project1-metadata.csv`)
- A tab-separated file containing the top 10 closest genomes for each isolate included in the analysis (e.g., `cool_project1-mash-top-10.tsv`)

# Examples
## Blastomyces
### Command
```
spree -g Blastomyces -o blasto -p blast-example -t 8 -l 2 example_genomes/Blastomyces_percursus_GCA_018296075.1.fasta
```
### UPGMA Mash Tree
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

## Streptococcus
### Command
```
spree -g Streptococcus -o strep -t 8 -l 2 Streptococcus_thermophilus_ATCC_19258.fasta Streptococcus_pyogenes_ATCC_12344.fasta
```
### UPGMA Mash Tree
<img src="https://github.com/johnjare/spree/blob/main/examples/strep-example-mash-tree.jpg" width="1000">

### Top 10 Closest Genomes
|Sample                                      |Reference       |Species                                                   | Mash Distance| % Matching Kmers|
|:-------------------------------------------|:---------------|:---------------------------------------------------------|-------------:|----------------:|
|Streptococcus_thermophilus_ATCC_19258.fasta |GCF_001008015.1 |[119] Streptococcus thermophilus (n=473)                  |     0.0136002|             60.2|
|Streptococcus_thermophilus_ATCC_19258.fasta |GCF_001280285.1 |[119] Streptococcus thermophilus (n=473)                  |     0.0136002|             60.2|
|Streptococcus_thermophilus_ATCC_19258.fasta |GCF_000971665.1 |[119] Streptococcus thermophilus (n=473)                  |     0.0137986|             59.8|
|Streptococcus_thermophilus_ATCC_19258.fasta |GCF_001514435.1 |[119] Streptococcus thermophilus (n=473)                  |     0.0150285|             57.4|
|Streptococcus_thermophilus_ATCC_19258.fasta |GCF_903886475.1 |[119] Streptococcus thermophilus (n=473)                  |     0.0151342|             57.2|
|Streptococcus_thermophilus_ATCC_19258.fasta |GCF_003438185.1 |[129] Streptococcus vestibularis (n=36)                   |     0.0651085|             14.6|
|Streptococcus_thermophilus_ATCC_19258.fasta |GCF_900636445.1 |[129] Streptococcus vestibularis (n=36)                   |     0.0653942|             14.5|
|Streptococcus_thermophilus_ATCC_19258.fasta |GCF_023110355.1 |[129] Streptococcus vestibularis (n=36)                   |     0.0662649|             14.2|
|Streptococcus_thermophilus_ATCC_19258.fasta |GCF_026781275.1 |[129] Streptococcus vestibularis (n=36)                   |     0.0674584|             13.8|
|Streptococcus_thermophilus_ATCC_19258.fasta |GCF_023110145.1 |[129] Streptococcus vestibularis (n=36)                   |     0.0696435|             13.1|
|Streptococcus_pyogenes_ATCC_12344.fasta     |GCF_000743015.1 |[105] Streptococcus pyogenes (n=2912)                     |     0.0114865|             64.7|
|Streptococcus_pyogenes_ATCC_12344.fasta     |GCF_900475035.1 |[105] Streptococcus pyogenes (n=2912)                     |     0.0115313|             64.6|
|Streptococcus_pyogenes_ATCC_12344.fasta     |GCF_000756485.1 |[105] Streptococcus pyogenes (n=2912)                     |     0.0118471|             63.9|
|Streptococcus_pyogenes_ATCC_12344.fasta     |GCF_000167435.2 |[105] Streptococcus pyogenes (n=2912)                     |     0.0119840|             63.6|
|Streptococcus_pyogenes_ATCC_12344.fasta     |GCF_000772185.1 |[105] Streptococcus pyogenes (n=2912)                     |     0.0137986|             59.8|
|Streptococcus_pyogenes_ATCC_12344.fasta     |GCF_012844385.1 |[33] Streptococcus dysgalactiae subsp. equisimilis (n=75) |     0.0751681|             11.5|
|Streptococcus_pyogenes_ATCC_12344.fasta     |GCF_012844405.1 |[33] Streptococcus dysgalactiae subsp. equisimilis (n=75) |     0.0755413|             11.4|
|Streptococcus_pyogenes_ATCC_12344.fasta     |GCF_013004125.1 |[33] Streptococcus dysgalactiae subsp. equisimilis (n=75) |     0.0755413|             11.4|
|Streptococcus_pyogenes_ATCC_12344.fasta     |GCF_008693725.1 |[31] Streptococcus dysgalactiae (n=71)                    |     0.0807479|             10.1|
|Streptococcus_pyogenes_ATCC_12344.fasta     |GCF_016128095.1 |[31] Streptococcus dysgalactiae (n=71)                    |     0.0816138|              9.9|
