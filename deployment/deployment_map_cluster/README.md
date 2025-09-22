## Cluster Evaluation of Metagenomic Results

# Running the Project

Begin by setting up the parameters for the project in the `deployment/params.json` file. This file contains paths to the necessary scripts and binaries, as well as parameters for the classifiers and clustering methods.

activate the conda environment and run the Nextflow workflow:

```bash
nextflow run main.nf -profile conda --params_file deployment/params.json
```

This will execute the workflow, simulating reads, classifying them, mapping them to references, clustering the results, and evaluating the clusters.

# Nextflow workflow

```mermaid
flowchart TB
    subgraph " "
    subgraph params
    v0["reads"]
    v1["input_table"]
    end
    v5([ExtractReferenceSequences])
    v12([MapMinimap2Paired])
    v14([FilterBamMsamtools])
    v16([sortBam])
    v18([SamtoolsCoverage])
    v24([MergeCoverageStatistics])
    v26([ClusterMappedReads])
    v28([MatchCladeReportWithReferenceSequences])
    v1 --> v5
    v0 --> v12
    v5 --> v12
    v12 --> v14
    v14 --> v16
    v16 --> v18
    v18 --> v24
    v12 --> v26
    v5 --> v28
    v24 --> v28
    v26 --> v28
    end
```
