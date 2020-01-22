FROM jupyter/datascience-notebook:72cca2a7f3ea

USER root

RUN apt-get update \
    && apt-get install -y --no-install-recommends zip \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

USER jovyan

RUN pip2 install --upgrade pip
RUN pip3 install --user --upgrade pip

# Install Python 2 & 3 Additional Packages for Notebooks
ADD requirements.txt requirements.txt
RUN pip2 install -r requirements.txt 
RUN /opt/conda/bin/python -m pip install -r requirements.txt

# Install R Packages used in Treehouse
RUN conda config --system --add channels r
RUN conda config --system --add channels bioconda
RUN conda config --system --add channels conda-forge
RUN conda install --yes -c r r-base=3.3.2 r-data.table=1.10.0 r-xml=3.98_1.5 
RUN conda clean -tipsy

RUN R -e "install.packages(c('openxlsx'), repos='http://cran.rstudio.com/')"

RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("biomaRt")'
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("pathview")'

# New scripts for the isoform report
RUN export TAR="/bin/tar" && R -e "require(devtools); install_version('vctrs', version = '0.2.0', repos = 'http://cran.us.r-project.org'); install_version('tidyverse', version = '1.3.0', repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages(c('treemapify', 'BiocManager'), repos='https://cran.cnr.berkeley.edu/');BiocManager::install('Sushi')"

USER root
RUN mkdir /work; chmod ugo+rwx /work
RUN chmod -R o+w /home/jovyan/.ipython
USER jovyan

ADD . /app

WORKDIR /outputs
ENTRYPOINT ["tini", "--", "python", "/app/setup.py"]
