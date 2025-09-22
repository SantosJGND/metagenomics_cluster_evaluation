## Benchmarking Cluster Evaluation of Metagenomic Results

# Running the Project

Begin by setting up the parameters for the project in the `deployment/params.json` file. This file contains paths to the necessary scripts and binaries, as well as parameters for the classifiers and clustering methods.

activate the conda environment and run the Nextflow workflow:

```bash
nextflow run main.nf -profile conda --params_file deployment/params.json
```

This will execute the workflow, mapping them to references, clustering the results, and evaluating the clusters.

# Nextflow workflow

```mermaid

flowchart TB
    subgraph " "
    v0["Taxid Table"]
    end
    v2([ExtractFastaSequences])
    v3([FormatToMess])
    v5([SimulateReadsMess])
    v7([CentrifugeClassificationPaired])
    v8([Kraken2ClassificationPaired])
    v9([MergeClassificationResults])
    v10([ExtractReferenceSequences])
    v14([MapMinimap2Paired])
    v23([FilterBamMsamtools])
    v15([sortBam])
    v16([SamtoolsCoverage])
    v21([MergeCoverageStatistics])
    v22([ClusterMappedReads])
    subgraph "Input"
    end
    v24([MatchCladeReportWithReferenceSequences])
    v25([EvaluateClusteringResults])
    v1(( ))
    v6(( ))
    v0 --> v1
    v1 --> v2
    v2 --> v3
    v1 --> v3
    v3 --> v5
    v5 --> v6
    v6 --> v7
    v7 --> v9
    v6 --> v8
    v8 --> v9
    v9 --> v10
    v9 --> v24
    v10 --> v24
    v10 --> v6
    v10 --> v14
    v6 --> v14
    v14 --> v23
    v23 --> v15
    v23 --> v22
    v15 --> v16
    v16 --> v21
    v21 --> v24
    v22 --> v24
    v24 --> v25
```
