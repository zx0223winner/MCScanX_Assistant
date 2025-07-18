[![Snakemake](https://img.shields.io/badge/snakemake-≥5.7.0-brightgreen.svg)](https://snakemake.readthedocs.io)
[![run with conda ](http://img.shields.io/badge/run%20with-conda%20-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)



# MCScanX_Assistant

## Overview
This repository stored the custom codes we used for the re-analysis of the protocol from Wang et al. Nature Protocols https://doi.org/10.1038/s41596-024-00968-2 (2024). 

It is intended to provide researchers with access to reproduce our work, facilitating future researchers to easily prepare the input data and install of MCScanX tool.


## Usage

Refer to [Usage](./docs/Usage.md) documents for details.

> [!NOTE]
> If you are new to Snakmake, please refer to [this page](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) on how to set-up SnakeMake. Make sure to test the example data below before running the workflow on your data.

```
# Test if you have successfully installed the SnakeMake
mamba activate snakemake
snakemake --help
```

> [!NOTE]
> Begin with a `config.yaml` file as below (detailed the input files requested).
> 
> The six species of NCBI assemblies are used in this analysis, including *B. carinata*, *A. suecica*, *A. arenosa*, *T. arvense*, *A. thaliana*, *B. oleracea*. 


## Arguments
**config.yaml**
```config.yaml

species_name:
  - Athaliana

ncbi_genomes:
    Athaliana:
        ncbi_assembly: "data/ncbi_download/GCF_000001735.4.zip"
        assembly_id: "GCF_000001735.4"
        feature_table: "data/ncbi_download/GCF_000001735.4_TAIR10.1_feature_table.txt.gz"
        species: Athaliana
```

> [!NOTE]
> Optional: To download extra ncbi assembly 'XX.zip' from NCBI, users can substitue the ncbi_assembly id (e.g., GCF_000001735.4) with yours in the command below:

```
curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000001735.4/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=GCF_000001735.4.zip"

```
## Running

Now, you can run the pipeline using the following commands:

```
# Download the package
git clone https://github.com/zx0223winner/MCScanX_Assistant.git

# enter the working directory
cd MCScanX_Assistant
```
> [!NOTE]
>Due to the size of sample files (we have prepared users with the standard input files of NCBI genome assemblies for the six species of NCBI assemblies are used in this analysis, including *B. carinata*, *A. suecica*, *A. arenosa*, *T. arvense*, *A. thaliana*, *B. oleracea*.  ).
>
> please download the test data - `MCScanX_Assistant_data.tar.gz` through the Google drive [link](https://drive.google.com/file/d/13KlaGXuVQysIAXoMjXtt2lHxVb-WNiP7/view?usp=sharing)
>
> (optional) please download the test result - `MCScanX_Assistant_results.tar.gz` through the Google drive [link](https://drive.google.com/file/d/1RorENiC0NPZhclForl9uU1pGLyHhgqVZ/view?usp=sharing). This file includes complete running results for users to check.

```
# Then decompress the file MCScanX_Assistant_data.tar.gz under the MCScanX_Assistant directory,
# This will bring you a data folder with test files ready 
tar -xvzf MCScanX_Assistant_results.tar.gz

# Then you can give a dry run by the following command.
snakemake --use-conda --cores all -s workflow/Snakefile_Input_preparing -n

# If everthing is OK, then you can test the pipeline by running one after another:
snakemake --use-conda --cores all -s workflow/Snakefile_Input_preparing
snakemake --use-conda --cores all -s workflow/Snakefile_Ks_distribution_plot
snakemake --use-conda --cores all -s workflow/Snakefile_MCScanX_6species
```



## License
The code within this repository is licensed under the [MIT License](./LICENSE). Please refer to the license file for more information on the terms and conditions of using and contributing to this project.

## Ciation
If you used the codes in this respository, please cite:

Xi Zhang, David Roy Smith, John M Archibald, Key step missing for input file preparation in MCScanX , Nature Protocols, (in submission).
