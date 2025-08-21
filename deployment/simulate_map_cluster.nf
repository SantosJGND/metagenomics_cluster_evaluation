workflow {

    params.input_table = params.input_table ?: ""
    if (params.input_table == "") {
        error("Input table path is not provided. Please set the 'input_table' parameter.")
    }
    input_table_ch = Channel.fromPath(params.input_table)
        .ifEmpty { error("Cannot find the input table: ${params.input_table}") }

    // Get reference sequences for each taxonomic ID in the input table
    matched_table_ch = ExtractFastaSequences(input_table_ch)
    input_table_ch = FormatToMess(input_table_ch, matched_table_ch.input_table_with_sequences)

    // Simulate reads based on the input table
    reads_ch = SimulateReadsMess("illumina", input_table_ch).map { table_id, fastq_files ->
        // Unpack the list of FASTQ files into separate variables
        def (fastq1, fastq2) = fastq_files
        tuple(table_id, fastq1, fastq2)
    }
    // Process reads using classifiers
    centrifuge_classification_ch = CentrifugeClassificationPaired(input_table_ch, reads_ch)
    kraken2_classification_ch = Kraken2ClassificationPaired(input_table_ch, reads_ch)
    merge_classification_results_ch = MergeClassificationResults(centrifuge_classification_ch, kraken2_classification_ch)

    // Extract reference sequences from the classification results
    reference_sequences_ch = ExtractReferenceSequences(input_table_ch, merge_classification_results_ch)

    // Map reads to reference sequences using minimap2
    flattened_reference_sequences_ch = reference_sequences_ch.reference_sequences.flatMap { ref_list ->
        ref_list
    }
    combined_ch = reads_ch.combine(flattened_reference_sequences_ch)

    mapped_reads_ch = MapMinimap2Paired(combined_ch, params.minimap2_illumina_params)

    // Extract mapping statistics from BAM files
    sorted_reads_ch = sortBam(mapped_reads_ch)
    coverage_ch = SamtoolsCoverage(sorted_reads_ch)

    // collect all mapping files, provide directory for clustering
    coverage_ch = coverage_ch.map { file, _filename, _refname, group -> tuple(group, file) }
    coverage_ch = coverage_ch.groupTuple()
    // collect all mapping files, provide directory for clustering
    mapping_files_info = mapped_reads_ch.map { file, group, _refname -> tuple(group, file) }
    mapping_files_info = mapping_files_info.groupTuple()

    // Merge mapping statistics
    merged_coverage_ch = MergeCoverageStatistics(coverage_ch)

    // Cluster mapped reads across alignment files
    clustering_ch = ClusterMappedReads(input_table_ch, mapping_files_info)

    MatchCladeReportWithReferenceSequences(
        input_table_ch,
        clustering_ch.clade_report,
        reference_sequences_ch.matched_assemblies,
        merged_coverage_ch.merged_coverage_statistics,
        merge_classification_results_ch,
    )
}

/*
* Match FormatSequences file to mess input table
*/
process FormatToMess {
    tag "FormatToMess ${input_table.baseName}"
    publishDir "${params.output_dir}/${input_table.baseName}/input", mode: 'symlink'

    input:
    path input_table
    path matched_table

    output:
    path "${input_table.baseName}.tsv", emit: formatted_input_table

    script:
    """
    #!/usr/bin/env python3
    import os
    import pandas as pd
    matched_table = pd.read_csv("${matched_table}", sep="\\t")
    if 'path' not in matched_table.columns:
        if 'assembly_file' not in matched_table.columns:
            raise ValueError("The matched table must contain a 'path' or 'assembly_file' column.")
        matched_table.rename(columns = {'assembly_file': 'path'}, inplace=True)
    matched_table['fasta'] = matched_table['path'].apply(lambda x: os.path.basename(x))
    matched_table.to_csv("${input_table.baseName}.tsv", sep="\\t", index=False)
    """
}

/*
* extract fasta sequences corresponding to the taxonomic IDs in the input_table
*/
process ExtractFastaSequences {
    tag "ExtractFastaSequences ${input_table.baseName}"
    publishDir "${params.output_dir}/${input_table.baseName}/input", mode: 'symlink'

    input:
    path input_table

    output:
    path "reference_sequences/matched_assemblies.tsv", emit: input_table_with_sequences

    script:
    """
    ${params.python_bin} ${params.references_extract_script} \
    --classification_output_path ${input_table} \
    --assembly_store "${params.assembly_store}" \
    --mapping_references_dir "reference_sequences" 
    """
}

/*
* merge coverage statistics from different BAM files
*/
process MergeCoverageStatistics {
    tag "MergeCoverageStatistics ${query_id}"

    publishDir "${params.output_dir}/${query_id}/coverage", mode: 'symlink'

    input:
    tuple val(query_id), path(coverage_files)

    output:
    path "merged_coverage_statistics.tsv", emit: merged_coverage_statistics

    script:
    def coverage_files_string = coverage_files.collect { it[0] }.join(',')
    """
    #!/usr/bin/env python3
    import pandas as pd
    import glob
    
    
    cov_files = "${coverage_files_string}".split(',')
    coverage_data = []
    for file in cov_files:
        df = pd.read_csv(file, sep="\\t")
        df['file'] = file.split('/')[-1]
        coverage_data.append(df)
    merged_df = pd.concat(coverage_data, ignore_index=True)

    merged_df.reset_index(inplace=True)
    merged_df.to_csv("merged_coverage_statistics.tsv", sep="\\t", index=False)
    """
}


/*
* Use samtools to extract BAM file coverage statistics
*/
process SamtoolsCoverage {
    tag "SamtoolsCoverage ${bamfile.baseName}"

    publishDir "${params.output_dir}/${query_id}/coverage", mode: 'symlink'

    input:
    tuple path(bamfile), path(bamindex), val(query_id), val(reference_id)

    output:
    tuple path("${bamfile.baseName}.coverage.txt"), val(bamfile.baseName), val(reference_id), val(query_id)

    script:
    """
    samtools coverage -o ${bamfile.baseName}.coverage.txt ${bamfile}
    """
}


/*
* sort and index bam file, maintain tuple file, query_id, reference_id in channel
*/
process sortBam {
    tag "sortMapping"
    publishDir "${params.output_dir}/${query_id}/sorted_reads", mode: 'symlink'

    input:
    tuple path(bamfile), val(query_id), val(reference_id)

    output:
    tuple path("${query_id}_${reference_id}.sorted.bam"), path("${query_id}_${reference_id}.sorted.bam.bai"), val(query_id), val(reference_id)

    script:
    """
    samtools sort ${bamfile} > ${query_id}_${reference_id}.sorted.bam
    samtools index ${query_id}_${reference_id}.sorted.bam
    samtools addreplacerg -r "ID:${query_id}" -r "SM:${query_id}" -o named.bam ${query_id}_${reference_id}.sorted.bam
    mv named.bam ${query_id}_${reference_id}.sorted.bam
    """
}


/*
* match clustering clade_report.txt and reference_sequences/matched_assemblies.tsv
*/
process MatchCladeReportWithReferenceSequences {
    tag "MatchCladeReportWithReferenceSequences ${input_table.baseName}"
    publishDir "${params.output_dir}/${input_table.baseName}/clustering", mode: 'symlink'

    input:
    path input_table
    path clade_report
    path matched_assemblies
    path coverage_report
    path merge_classification_results

    output:
    path "clade_report_with_references.tsv", emit: clade_report_with_references

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    clade_report = pd.read_csv("${clade_report}", sep="\\t", header=None, names=["clade", "nuniq", "freq", "files"])
    clade_report['clade']
    clade_report['files'] = clade_report['files'].str.split(',')
    clade_report = clade_report.explode('files')

    matched_assemblies = pd.read_csv("${matched_assemblies}", sep="\\t")
    matched_assemblies['filename'] = matched_assemblies['assembly_file'].str.split('/').str[-1]

    coverage_report = pd.read_csv("${coverage_report}", sep="\\t")
    merged_classification_results = pd.read_csv("${merge_classification_results}", sep="\\t")

    def find_assembly_mapping(row):
        accession = row['assembly_accession']
        match = clade_report[clade_report['files'].str.contains(accession, na=False)]
        if match.empty:
            row['clade'] = 'unmapped'
            row['nuniq'] = 0
            row['freq'] = 0
        else:
            row['clade'] = match['clade'].values[0]
            row['nuniq'] = match['nuniq'].values[0]
            row['freq'] = match['freq'].values[0]
        
        return row
    
    def find_assembly_coverage(row):
        accession = row['assembly_accession']
        match = coverage_report[coverage_report['file'].str.contains(accession, na=False)]
        if match.empty:
            row['coverage'] = 0
        else:
            row['coverage'] = match['coverage'].values[0]
        
        return row
    
    def find_assembly_classification(row):
        taxid = row['taxid']
        match = merged_classification_results[merged_classification_results['taxid'] == taxid]
        if match.empty:
            row['classifier'] = 'unclassified'
        else:
            row['classifier'] = match['classification'].values[0]
        
        return row
    
    clade_report_with_references = matched_assemblies.apply(find_assembly_mapping, axis=1)
    clade_report_with_references = clade_report_with_references.apply(find_assembly_classification, axis=1)
    clade_report_with_references = clade_report_with_references.apply(find_assembly_coverage, axis=1)
    clade_report_with_references = clade_report_with_references[['description', 'taxid', 'assembly_accession', \
            'coverage', 'clade', 'nuniq', 'freq', 'classifier']]

    clade_report_with_references.to_csv("clade_report_with_references.tsv", sep="\\t", index=False)
    """
}

/*
* Merge mapping statistics from different BAM files
*/
process MergeMappingStatistics {
    tag "MergeMappingStatistics ${input_table.baseName}"
    publishDir "${params.output_dir}/${input_table.baseName}/mapping_stats", mode: 'symlink'

    input:
    val input_table
    path flagstat_files

    output:
    path "merged_mapping_stats.tsv"

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import glob

    flagstat_files = glob.glob("${flagstat_files}/*.txt")
    dataframes = []
    for file in flagstat_files:
        df = pd.read_csv(file, sep="\\t", header=None, names=["stat", "value"])
        df['file'] = file.split('/')[-1]
        dataframes.append(df)
    merged_df = pd.concat(dataframes, ignore_index=True)
    merged_df = merged_df.pivot(index='file', columns='stat', values='value')
    merged_df.reset_index(inplace=True)
    merged_df.to_csv("merged_mapping_stats.tsv", sep="\\t", index=False)
    """
}

/*
* Extract mapping statistics from BAM files
*/
process ExtractMappingStatistics {
    publishDir "${params.output_dir}/${query_id}/mapping_stats", mode: 'symlink'

    input:
    tuple path(bam_file), val(query_id), val(reference_name)

    output:
    path "${bam_file.baseName}.txt"

    script:
    """
    samtools flagstat -O tsv ${bam_file} > ${bam_file.baseName}.txt
    """
}

/*
* cluster mapped reads accross alignment files
*/
process ClusterMappedReads {
    tag "ClusterMappedReads ${input_table.baseName}"

    publishDir "${params.output_dir}/${input_table.baseName}", mode: 'symlink'

    input:
    path input_table
    tuple val(query_id), path(mapped_reads)

    output:
    path "clustering/clade_report.txt", emit: clade_report
    path "clustering/sample_report.txt", emit: sample_report

    script:
    def mapped_reads_string = mapped_reads.collect { it[0] }.join(',')
    """
    ${params.map_to_matrix_bin} \
    --files ${mapped_reads_string} \
    -o clustering \
    ${params.map_to_matrix_params} 
    """
}


/*
* Map to a reference using minimap2
*/
process MapMinimap2Paired {
    tag "MapMinimap2Paired ${fastq1} ${fastq2} ${reference.baseName}"

    publishDir "${params.output_dir}/${query_id}/mapped_reads", mode: 'symlink'

    input:
    tuple val(query_id), path(fastq1), path(fastq2), path(reference)
    val minimap2_params

    output:
    tuple path("${fastq1.baseName}_${reference.baseName}.bam"), val(query_id), val(reference.baseName)

    script:
    """
    mkdir -p ${query_id}
    minimap2 ${minimap2_params} -ax sr ${reference} ${fastq1} ${fastq2} | samtools view -bS -F 4 - > ${fastq1.baseName}_${reference.baseName}.bam
    """
}

/*
* Extract reference sequences from classifier output results to global reference database
*/
process ExtractReferenceSequences {
    tag "ExtractReferenceSequences ${input_table.baseName}"

    publishDir "${params.output_dir}/${input_table.baseName}", mode: 'symlink'

    input:
    path input_table
    path classifier_output

    output:
    path "reference_sequences/*gz", emit: reference_sequences
    path "reference_sequences/matched_assemblies.tsv", emit: matched_assemblies

    script:
    """
    ${params.python_bin} ${params.references_extract_script} \
    --classification_output_path ${classifier_output} \
    --assembly_store "${params.assembly_store}" \
    --mapping_references_dir "reference_sequences" 
    """
}


/*
* Merge classification results from different classifiers
*/
process MergeClassificationResults {
    tag "MergeClassificationResults ${query_id}"
    publishDir "${params.output_dir}/${query_id}/classification", mode: 'symlink'

    input:
    path centrifuge_output
    path centrifuge_output_processed
    val query_id
    path kraken2_output
    path kraken2_output_processed
    val query_id_2

    output:
    path "${query_id}_merged_classification.tsv"

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    centrifuge_file = "${centrifuge_output_processed}"
    kraken2_file = "${kraken2_output_processed}"  
    query_id = "${query_id}"
    centrifuge_df = pd.read_csv("${centrifuge_output_processed}", sep="\\t").rename(columns={"taxID": "taxid", "name": "description"})
    kraken2_df = pd.read_csv("${kraken2_output_processed}", sep="\\t").rename(columns={"taxID": "taxid", "name": "description"})
    kraken2_df['classifier'] = 'kraken2'
    centrifuge_df['classifier'] = 'centrifuge'
    merged_df = pd.merge(centrifuge_df, kraken2_df, on=["description", "taxid"], how="outer", suffixes=("_centrifuge", "_kraken2"))
    merged_df = merged_df.drop_duplicates(subset=["taxid"])
    def classify(row):
        if pd.notna(row['classifier_kraken2']):
            if pd.notna(row['classifier_centrifuge']):
                    return 'kraken2'
            else:
                return 'kraken2/centrifuge'
        elif pd.notna(row['classifier_centrifuge']):
            return 'centrifuge'
        else:
            return 'unclassified'
    
    merged_df['classification'] = merged_df.apply(classify, axis=1)
    merged_df.to_csv(f"${query_id}_merged_classification.tsv", sep="\\t", index=False)
    """
}

/*
* Process output from classifiers using python script
*/
process ProcessClassifierOutput {
    tag "ProcessClassifierOutput ${query_id}"

    publishDir "${params.output_dir}/${query_id}/classification/${classifier_type}", mode: 'symlink'

    input:
    path classifier_output
    val query_id
    val classifier_type

    output:
    tuple path("${query_id}_${classifier_type}_processed_classifier_output.tsv"), val(query_id)

    script:
    """
    ${params.python_bin} ${params.classifier_process_script} \
    --input ${classifier_output} \
    --output ${query_id}_${classifier_type}_processed_classifier_output.tsv \
    --type ${classifier_type} 
    """
}


/*
* Run Kraken2 classification on paired-end reads.
*/
process Kraken2ClassificationPaired {
    tag "Kraken2Classification ${input_table.baseName}"
    publishDir "${params.output_dir}/${input_table.baseName}/classification/kraken2", mode: 'symlink'

    input:
    path input_table
    tuple val(table_id), path(fastq1), path(fastq2)

    output:
    path "${input_table.baseName}_kraken2_report.txt"
    path "${input_table.baseName}_krk2_processed_classifier_output.tsv"
    val input_table.baseName

    script:
    """
    ${params.kraken2_bin} --db ${params.kraken2_index} \
    --report ${input_table.baseName}_kraken2_report.txt \
    --gzipped-compressed \
    --output ${input_table.baseName}_kraken2_classification.txt \
    ${fastq1} ${fastq2} ${params.kraken2_params}

    ${params.python_bin} ${params.classifier_process_script} \
    --input "${input_table.baseName}_kraken2_report.txt" \
    --output "${input_table.baseName}_krk2_processed_classifier_output.tsv" \
    --type "kraken2" \
    --nuniq_threshold ${params.minimum_uniq_reads}

    """
}

/*
* Classify reads using Centrifuge in paired-end mode
*/
process CentrifugeClassificationPaired {
    tag "CentrifugeClassificationSingle ${input_table.baseName}"

    publishDir "${params.output_dir}/${input_table.baseName}/classification/centrifuge", mode: 'symlink'

    input:
    path input_table
    tuple val(table_id), path(fastq1), path(fastq2)

    output:
    path "${input_table.baseName}_centrifuge_report.txt"
    path "${input_table.baseName}_centrifuge_processed_classifier_output.tsv"
    val "${input_table.baseName}"

    script:
    """
    centrifuge -x ${params.centrifuge_index} -1 ${fastq1} -2 ${fastq2} \
    -S ${input_table.baseName}_centrifuge_classification.tsv \
    --output ${input_table.baseName}_centrifuge_classification.txt \
    --report-file ${input_table.baseName}_centrifuge_report.txt \
    ${params.centrifuge_params}

    ${params.python_bin} ${params.classifier_process_script} \
    --input "${input_table.baseName}_centrifuge_report.txt" \
    --output "${input_table.baseName}_centrifuge_processed_classifier_output.tsv" \
    --type "centrifuge" \
    --nuniq_threshold ${params.minimum_uniq_reads}  

    """
}

/*
* Simulate reads using the mess package and conda environment
*/
process SimulateReadsMess {
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    val technology
    path input_table

    output:
    tuple val("${input_table.baseName}"), path("${input_table.baseName}/fastq/*.fq.gz")

    script:
    """
    mess simulate --input ${input_table} \
    --output "${input_table.baseName}" \
    --threads 3 \
    --tech ${technology} \
    --bam
    """
}
