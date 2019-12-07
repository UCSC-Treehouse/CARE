<table><tr><td><img src="/docs/CARE-asklepion-v2.png"></td>
 <td><h1> Treehouse Comparative Analysis of RNA Expression (CARE)</h1><br>
 <a href="http://treehousegenomics.ucsc.edu">UCSC Treehouse Childhood Cancer Initiative</a>
</td></tr></table>

## Documentation on the docker used in the [Vaske 2019 comparative tumor RNA-Seq analysis publication](https://treehousegenomics.ucsc.edu/p/vaske-2019-comparative-tumor-RNA/) is available [here](docs/Vaske-2019-comparative-tumor-RNA.md).

Treehouse Comparative Analysis of RNA Expression (CARE) is a dockerized pipeline of Jupyter notebooks that performs outlier analysis on an RNA-Seq expression file generated from a FASTQ by the UCSC
[UCSC Treehouse Toil RNA-Seq expression pipeline](https://github.com/UCSC-Treehouse/pipelines), vs a background compendium of samples. It generates two HTML documents summarizing the results that you can populate with your sample's clinical information.

## Requirements
```
Linux
Git
Make
Docker
6.5 Gb of disk space, approximately, for background compendium and reference files
200 Mb disk space per focus sample for combined input+output
```

## Usage
Clone:
```
git clone https://github.com/UCSC-Treehouse/CARE.git
cd CARE
```

Checkout the [latest release](https://github.com/UCSC-Treehouse/CARE/releases):

```git checkout 0.16.0.0```

Use Make to run the docker:
```make run```

This will:
 - Pull the docker
 - Download the background compendium and reference files if not already present. (May take several minutes; compendium download is ~1.5GB).
 - Download a test sample
 - Run CARE on the test sample

### Manual downloads
If you prefer, you can manually download the input tgz files from the following URLs:
- https://xena.treehouse.gi.ucsc.edu/download/CARE/TumorCompendium_v10_PolyA.tgz
- https://xena.treehouse.gi.ucsc.edu/download/CARE/TreehouseReferences-2019-10-23.tgz
- https://xena.treehouse.gi.ucsc.edu/download/CARE/TreehouseExampleFocusSample-TH03_0013_S02-2019-10-23.tgz
- https://xena.treehouse.gi.ucsc.edu/download/CARE/TreehouseExampleManifest-2019-10.24.tgz

MD5sums:
```
cfafeb5ff7b93591d0a80894061b5f3e  TumorCompendium_v10_PolyA.tgz
d616135b19a31f8c7a0cab4932d2fd01  TreehouseReferences-2019-10-23.tgz
4bc8001b744d5ea6a8d9a94db06a7df9  TreehouseExampleFocusSample-TH03_0013_S02-2019-10-23.tgz
03bb977bb9610a19c84e98a788a96d41  TreehouseExampleManifest-2019-10.24.tgz
```

To use, place them in the "resources" folder as the original tgz archive without extracting.

## Documentation
Docs are available [here](/docs).

## Contact
You can open an issue with your question or get in touch via the [Treehouse contact form](https://treehousegenomics.soe.ucsc.edu/contact/).
