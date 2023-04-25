# libraries
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(treeio))
suppressPackageStartupMessages(library(ggtree))

###---START MESSAGE---###
cat('Thank you for using spree!', sep = "\n")

args <- commandArgs(trailingOnly = TRUE)
taxon <- args[1]
input <- args[2]

###---CHECK INPUT---###
if (!file.exists(input)) {
  # if file doesn't exist, print error message and exit program
  cat("Error: input file not found.\n")
  q("no")
}
input.name <- input %>%
  gsub(pattern = ".*/", replacement = "") %>%
  gsub(pattern = "\\.*$", replacement = "")

###---GET ASSEMBLIES SUMMARY FOR TAXON---###
# run NCBI-Datasets
system(paste0("docker run staphb/ncbi-datasets datasets summary genome taxon ",
             taxon,
" --as-json-lines > ncbi.json && docker run -v $PWD:/data/ staphb/ncbi-datasets dataformat tsv genome --inputfile /data/ncbi.json --fields accession,assmstats-genome-coverage,organism-name,assmstats-total-sequence-len,assminfo-status,assminfo-level,organism-tax-id,source_database > ncbi.tsv && rm ncbi.json"))

# load result
## check for output
# check if file exists
if (!file.exists('ncbi.tsv')) {
  # if file doesn't exist, print error message and exit program
  cat("Error: ncbi.tsv not found.\n")
  q("no")
} else {
  df.ncbi <- read_tsv('ncbi.tsv', show_col_types = FALSE)
}

###---CLEAN NCBI SUMMARY---###
# select only current assemblies
df.ncbi <- df.ncbi %>%
  subset(`Assembly Status` == "current")

# extract species names
extract_species <- function(row) {
  words <- strsplit(row, " ")[[1]]
  if (grepl(pattern = "subsp", row)){
    return(paste(words[1:4], collapse = " "))
  }
  if (length(words) >= 2) {
    return(paste(words[1:2], collapse = " "))
  } else {
    return(NA)
  }
}
df.ncbi$Species <- unlist(lapply(df.ncbi$`Organism Name`, FUN=extract_species))

# remove unwanted species
df.ncbi <- df.ncbi[grep(x=df.ncbi$Species, pattern= "(uncultured| sp.| spp.)",invert=T),]

# total counts to species names
df.ncbi <- df.ncbi %>%
  group_by(Species) %>%
  count() %>%
  merge(df.ncbi, by = "Species") %>%
  mutate(Species = paste0(Species," (n=",n,")"))

# add index to each species name and accession name
n_species <- unique(df.ncbi$Species)
index.table <- data.frame(Species = n_species, index = seq(from=1, to=length(n_species))) %>%
  mutate(index_species = paste0("[",index,"] ",Species))
df.ncbi <- df.ncbi %>%
  merge(index.table, by="Species") %>% 
  mutate(Species = factor(index_species, levels = unique(index_species)),
         name = paste0("[",index,"] ",`Assembly Accession`))

# clean up coverage values
df.ncbi$`Assembly Stats Genome Coverage` <- as.numeric(gsub(pattern = "x", replacement = "", df.ncbi$`Assembly Stats Genome Coverage`))

###---ASSIGN QUALITY STATUS---###
# coverage status
df.ncbi <- df.ncbi %>%
  mutate(Coverage_Status = case_when(`Assembly Stats Genome Coverage` > 30 ~ TRUE,
                                     TRUE ~ FALSE))
# database status
df.ncbi <- df.ncbi %>%
  mutate(Database_Status = case_when(`Source Database` == 'SOURCE_DATABASE_REFSEQ' ~ TRUE,
                                     TRUE ~ FALSE))
# complete status
df.ncbi <- df.ncbi %>%
  mutate(Complete_Status = case_when(`Assembly Level` == 'Complete Genome' ~ TRUE,
                                     TRUE ~ FALSE))
# outlier status
df.ncbi <- df.ncbi %>% 
  group_by(Species) %>%
  mutate(gl_mean = mean(`Assembly Stats Total Sequence Length`),
         gl_sd = sd(`Assembly Stats Total Sequence Length`)) %>%
  mutate(Outlier_Status = case_when(abs(gl_mean - `Assembly Stats Total Sequence Length`) <= gl_sd*3 ~ TRUE,
                                    TRUE ~ FALSE)) %>%
  select(-gl_mean, -gl_sd) %>%
  ungroup()

# create summary
df.ncbi$Status_Counts <- df.ncbi$`Coverage_Status`+df.ncbi$`Database_Status`+df.ncbi$`Complete_Status`+df.ncbi$`Outlier_Status`
df.ncbi <- df.ncbi %>%
  group_by(`Assembly Accession`) %>%
  mutate(Overall_Status = paste(rep('*', Status_Counts), collapse = '')) %>%
  select(-Status_Counts)

###---SELECT REFEREANCES---###
select_ref <- function(species){
  # subset by group
  subdf <- df.ncbi[df.ncbi$Species == species,]
  # select up to 5 assemblies with four stars
  n_refs <- 0
  if(nrow(subdf[subdf$Overall_Status == '****',]) > 0){
    assemblies <- subdf[subdf$Overall_Status == '****',]$`Assembly Accession`
    n_total <- length(assemblies)
    if(n_total > 5){
      refs_4 <- assemblies[1:5]
    }else(refs_4 <- assemblies[1:n_total])
    n_refs <- length(refs_4)
  }else(refs_4 <- NA)
  # if more assemblies are needed, select assemblies with three stars up to 5
  if(nrow(subdf[subdf$Overall_Status == '***',]) > 0){
    assemblies <- subdf[subdf$Overall_Status == '***',]$`Assembly Accession`
    n_total <- length(assemblies)
    needed <- 5-n_refs
    if(n_total > needed){
      refs_3 <- assemblies[1:needed]
    }else(refs_3 <- assemblies[1:n_total])
    n_refs <- n_refs+length(refs_3)
  }else(refs_3 <- NA)
  # if more assemblies are still needed, select assemblies with two stars up to 5
  if(nrow(subdf[subdf$Overall_Status == '**',])){
    assemblies <- subdf[subdf$Overall_Status == '**',]$`Assembly Accession`
    n_total <- length(assemblies)
    needed <- 5-n_refs
    if(n_total > needed){
      refs_2 <- assemblies[1:needed]
    }else(refs_2 <- assemblies[1:n_total])
    n_refs <- n_refs+length(refs_2)
  }else(refs_2 <- NA)
  # if more assemblies are still needed, select assemblies with one stars up to 5
  if(nrow(subdf[subdf$Overall_Status == '*',]) > 0){
    assemblies <- subdf[subdf$Overall_Status == '*',]$`Assembly Accession`
    n_total <- length(assemblies)
    needed <- 5-n_refs
    if(n_total > needed){
      refs_1 <- assemblies[1:needed]
    }else(refs_1 <- assemblies[1:n_total])
    n_refs <- n_refs+length(refs_1)
  }else(refs_1 <- NA)
  ## combine all selections and extract dataframe
  refs <- c(refs_4,refs_3,refs_2,refs_1)
  refs <- refs[!is.na(refs)]
  if(length(refs) > 0){
    result <- subdf[subdf$`Assembly Accession` %in% refs,]
  }
  return(result)
}

df.select <- do.call(rbind, lapply(unique(df.ncbi$Species), select_ref))
write.table(x = df.select, 'select-refs.tsv', row.names = F, quote = F, sep = '\t')

###---DOWNLOAD REFERENCES---###
# remove existing 'refs' directory and make new
if(dir.exists('refs')){
  unlink('refs')
}
dir.create('refs')
# function for downloading reference genomes
download_refs <- function(a){
  system(paste("docker run -v $PWD:/data/ staphb/ncbi-datasets datasets download genome accession ",
               a,
               " && unzip ncbi_dataset.zip && cp ncbi_dataset/data/*/*.fna refs && rm -r ncbi_dataset* README.md"))
}

#lapply(df.select$`Assembly Accession`, FUN=download_refs)

###---RUN MASH---###
system(paste0("docker run -v $PWD:/data/ staphb/mash mash sketch -p 4 -o sketch.msh refs/* ",input," && docker run -v $PWD:/data/ staphb/mash mash dist -p 4 sketch.msh sketch.msh > mash-ava.tsv"))

###---MAKE TREE---###
# prepare metadata
input.name <- str_split(input, pattern = "/") %>% unlist()
input.meta <- data.frame("Assembly Accession" = input, 
                         name=input.name[max(length(input.name))], 
                         Species = NA, 
                         Overall_Status = NA, face = "bold") %>%
  rename("Assembly Accession"=Assembly.Accession)
meta <- df.select %>%
  select(`Assembly Accession`, name, Species, Overall_Status) %>%
  mutate(name = paste(name,Overall_Status, sep=" "),
         face="plain") %>%
  rbind(input.meta)

write.csv(meta, file="meta.csv")
# load mash file
mash.all <- read_tsv("mash-ava.tsv", col_names = F, show_col_types = FALSE)[,1:3]

# clean up names
clean_names <- function(name){
  name <- gsub(pattern = ".*/", replacement = "", name) %>%
    gsub(pattern = "\\.*$", replacement = "")
  if(grepl(pattern = "^(GCA|GCF)_.*", name)){
    name <- substr(start = 1, stop = 15, name)
  }
  return(name)
}
mash.all$X1 <- lapply(mash.all$X1, FUN=clean_names) %>% unlist()
mash.all$X2 <- lapply(mash.all$X2, FUN=clean_names) %>% unlist()

# get top hit info
mash.top.10 <- mash.all %>%
  subset(X1 == input.name) %>%
  subset(X1 != X2) %>% 
  rename("Assembly Accession"=X2) %>%
  rename("Mash Distance"=X3) %>%
  merge(df.ncbi[,c("Assembly Accession","Species")]) %>%
  select(`Assembly Accession`, Species, `Mash Distance`)
mash.top.10 <- mash.top.10[order(mash.top.10$`Mash Distance`),]
mash.top.10[1:10,] %>%
  knitr::kable(row.names = F)
# convert to matrix
mash.all <- mash.all %>% unique() %>% spread(key = "X2", value = "X3") %>% column_to_rownames(var="X1") %>% as.matrix()

# create distance matrix
dist <- dist(mash.all, method = "euclidean")
# create tree and save
upgma <- phangorn::upgma(dist)
write.tree(phy = upgma, file = "mash-ava.tree")

# make tree image
## load tree
mash.tree <- read.tree("mash-ava.tree")

## initial plot
p_mash.tree <- ggtree(mash.tree)
## get sizing info
size_tree <- function(plot){
  # extract plot data
  data <- plot$data
  # determine limits
  ## max sample name length
  max_name <- nchar(data$label) %>% max(na.rm = T)
  if(max_name > 100){
    max_name <- 30
  }
  ## max x-coordinate
  max_x <- max(data$x)
  return(list(max_name,max_x))
}
maxs <- size_tree(p_mash.tree) %>% unlist()
name_size=as.numeric(maxs[1])/15
x_max=as.numeric(maxs[2])*1.2
## re-plot tree
p_mash.tree <- p_mash.tree %<+% meta+ 
  geom_tiplab(aes(label=name, fontface=face),size=name_size, hjust = -0.2)+
  geom_tippoint(aes(color=Species))+
  xlim(0,as.numeric(x_max))+
  ggtitle("UPGMA Mash Tree")
## save image
n_iso <- p_mash.tree$data %>%
  drop_na() %>%
  nrow()
# set image dimensions
wdth <- n_iso/5
if(wdth < 6){
  wdth <- 6
}
hght=n_iso/5
if(hght < 6){
  hght <- 6
}

ggsave(plot = p_mash.tree, filename = "mash-tree.jpg", width = wdth, height = hght, dpi = 300, limitsize = F)