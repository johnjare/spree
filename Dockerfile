# Use a base image with Ubuntu
FROM bioconductor/bioconductor_docker

# required R packages
RUN R -e "install.packages(c('tidyverse', 'ggnewscale', 'phangorn','BiocManager','RColorBrewer','progress'), repos='https://cran.rstudio.com/')"
RUN R -e "BiocManager::install('treeio'); BiocManager::install('ggtree')"

# Install Mash
RUN wget https://github.com/marbl/Mash/releases/download/v2.3/mash-Linux64-v2.3.tar && tar -xvf mash-Linux64-v2.3.tar && rm mash-Linux64-v2.3.tar && mv mash-Linux64-v2.3/mash bin/

# Install NCBI datasets
RUN wget https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets && wget https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/dataformat && chmod +x datasets dataformat && mv datasets bin/ && mv dataformat bin/

# Clone the spree repository
RUN git clone https://github.com/johnjare/spree.git && mv spree/spree bin/
