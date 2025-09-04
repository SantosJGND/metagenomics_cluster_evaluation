## Mapping Cluster

Metagenomic analysis creates a wealth of data that can be challenging to interpret. One of the main challenges in this field is to discriminate between Cross-Hits - Metagenomic classifiers often return multiple, related, hits. This makes it difficult to determine which hit is the most relevant for a given sample. This project aims to benchmark effectiveness of clustering in resolving these Cross-Hits.

## Project Structure

The project is structured to facilitate the evaluation of clustering methods on metagenomic data. The main components include:

- **simulation**: Simulating reads from reference assemblies to create a controlled dataset for testing.
- **classification/**: Classifying metagenomic reads using Centrifuge and Kraken2 classifiers.
- **reference_management/**: Managing reference assemblies, including fetching representative assemblies from NCBI, local db management, and checking assembly existence.
- **mappings**: Mapping classified reads to reference assemblies to evaluate the accuracy of classifications.
- **clustering/**: Clustering classified reads and evaluating the clusters.
- **reference_management/**: Contains scripts for managing reference assemblies and their taxonomic classifications.

### Deployments

Directories containing deployment scripts and parameters for running

- [**deployment_benchmark/**](deployment_benchmark): Simulate reads, map classified reads to reference assemblies, cluster mapped reads.

- [**deployment_map_cluster/**](deployment_map_cluster): Map reads onto table of provided references, cluster the results.

## Requirements

- Nextflow

- A conda environment with the the following packages:
  - `biopython`
  - `mess` # simulation only

This environment must be activated before running the nextflow workflow.

- a python environment with the following packages:

  - `pandas`
  - `numpy`
  - `Biopython`

- Metagenomic classifiers:
  - Centrifuge
  - Kraken2
    and their respective databases.

# Running the Project

Begin by setting up the parameters for the project in the `deployment/params.json` file. This file contains paths to the necessary scripts and binaries, as well as parameters for the classifiers and clustering methods.

activate the conda environment and run the Nextflow workflow:

```bash
nextflow run main.nf -profile conda --params_file deployment/params.json
```

This will execute the workflow, simulating reads, classifying them, mapping them to references, clustering the results, and evaluating the clusters.

# Reference extraction

References used for mapping arep stored locally, in an local directory specified in the `params.json` file as the `assembly_store` directory. reference identification is performed using the `Entrez` tool suite, as implemented in the `Bio.Entrez` python package.

### Databases

The software will first check for availability for a specified taxid in the `nucleotide` reference database. If not found, it will resort to the `assembly` database.

## Checking for software behaviour

The reference management module will download references on the fly during workflow deployment. It might be useful
to first ensure that the references will be available, or if the references to be downloaded are appropriate.

For this, run the `check`submodule of the `references_management/main.py` script on the prospective taxid table:
flow, simulating re
```bash
python references_management/main.py check --input_table path/to/taxid_table.tsv --assessment /path/to/assessment.tsv
```

The output `assessment.tsv` file will contain information on the availability of references for each taxid in the input table.
