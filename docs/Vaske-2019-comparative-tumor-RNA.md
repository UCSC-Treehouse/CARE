
# Vaske 2019 Comparative Tumor RNA-Seq

[Resources Hub](https://treehousegenomics.ucsc.edu/p/vaske-2019-comparative-tumor-RNA/)

The docker used to generate outlier results for the Vaske 2019 comparative Tumor RNA-Seq publication is available for download
on [Docker hub](https://cloud.docker.com/u/ucsctreehouse/repository/docker/ucsctreehouse/care):

    ucsctreehouse/care:0.8.1.0

# Usage

Download the background compendium and reference files:
[vaske-2019-comparative-tumor-RNA-CARE-2019-09-14.tgz](https://xena.treehouse.gi.ucsc.edu/download/CARE/vaske-2019-comparative-tumor/vaske-2019-comparative-tumor-RNA-CARE-2019-10-24.tgz)

MD5sum: 401c63e6a99a9b8a6d2d3cce5464bd7c

Size: 4.3G when uncompressed

Place your rsem_genes.results  and diagnosed_disease.txt files in the input/ directory structure:

    inputs/
    └── SAMPLEID
        ├── diagnosed_disease.txt
        └── expression
            └── RSEM
                └── rsem_genes.results


rsem_genes.results :  Generated from the [RNASeq expression pipeline](https://github.com/UCSC-Treehouse/pipelines).
diagnosed_disease.txt : Plain text containing the diagnosis of your sample

Next, you must have an email address registered with the Broad GSEA/MSigDB as the docker uses this site to retrieve overexpressed pathways.

To register (free): http://software.broadinstitute.org/gsea/register.jsp

Then, to process all samples in input/, depositing the results into output/:

    make run EMAIL=user@example.edu

(replacing the email address with your registered address)


## Reproducibility notes

The provided background compendium (v5) contains expression data that is not fully compatible with the [latest version of CARE](https://github.com/UCSC-Treehouse/CARE).
Specifically: CARE works in a HUGO gene name space, while the [Toil RNA-Seq pipeline](https://github.com/UCSC-Treehouse/pipelines) that generates the RNA expression levels
names the genes using Ensembl IDs. This results in a conversion step where the code combines the expression values of multiple Ensembl IDs that map to the same HUGO name.
For this release (v5), the resulting value was the mean of the inputs; however, we have since corrected the algorithm to sum the inputs instead. We have confirmed that the outliers reported in Vaske 2019 are not affected by this change, but we suggest using the latest CARE code for any new samples to be analyzed.

In addition, when using this v5 compendium and pipeline, the results involving MSigDB drug-gene database or  Broad "compute overlaps" API may differ slightly from published due to updates to the gene sets since time of original analysis.
In addition, when performing overlap analysis only the 2025 genes with highest expression in a gene set will be used. This is due to restrictions on the Broad overlap analysis API which were not present at time of original analysis. 
The current version of CARE does not use this API and so is not subject to these restriction.
