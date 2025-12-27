params.analysis_id = params.analysis_id ?: UUID.randomUUID().toString()

workflow {

    analysis_id = params.analysis_id


    if (params.reads == "") {
        error("Reads directory path is not provided. Please set the 'reads' parameter.")
    }


    Channel
        .fromPath(["${params.reads}/*1.fq.gz", "${params.reads}/*1.fastq.gz"])
        .map { file -> tuple(file) }
        .ifEmpty { error('Cannot find any paired-end fastq files') }
        .set { reads_ch1 }

    Channel
        .fromPath(["${params.reads}/*2.fq.gz", "${params.reads}/*2.fastq.gz"])
        .map { file -> tuple(file) }
        .ifEmpty { file('NO_FILE') }
        .set { reads_ch2 }

    reads_ch = reads_ch1.combine(reads_ch2)

    qc = QCReadsPrinseqPaired(analysis_id, reads_ch)
    // print
    r1_ch = qc._good_R1
    r2_ch = qc._good_R2

    reads_ch = r1_ch
        .combine(r2_ch.ifEmpty { file('NO_FILE') })
        .map { r1, r2 -> tuple(r1, r2) }

    // Process reads using classifiers
    centrifuge_classification_ch = CentrifugeClassificationPaired(analysis_id, reads_ch)
    kraken2_classification_ch = Kraken2ClassificationPaired(analysis_id, reads_ch)
    diamond_classification_ch = DiamondClassificationPaired(analysis_id, reads_ch)
    krakenunique_classification_ch = KrakenUniqueClassificationPaired(analysis_id, reads_ch)
    merge_classification_results_ch = MergeClassificationResults(
        centrifuge_classification_ch,
        kraken2_classification_ch,
        krakenunique_classification_ch,
        diamond_classification_ch,
    )
    // Extract reference sequences from the classification results
    reference_sequences_ch = ExtractReferenceSequences(analysis_id, merge_classification_results_ch)

    // Check if reference sequences are empty and end workflow gracefully if so
    flattened_reference_sequences_ch = reference_sequences_ch.reference_sequences.flatMap { ref_list -> ref_list }

    combined_ch = reads_ch.combine(flattened_reference_sequences_ch)

    mapped_reads_ch = MapMinimap2Paired(analysis_id, combined_ch, params.minimap2_illumina_params)

    // Extract mapping statistics from BAM files
    filtered_alignments_ch = FilterBamMsamtools(mapped_reads_ch)
    sorted_reads_ch = sortBam(filtered_alignments_ch)
    coverage_ch = SamtoolsCoverage(sorted_reads_ch)

    // collect all mapping files, provide directory for clustering
    coverage_ch = coverage_ch.map { file1, file2, _filename, _refname, group -> tuple(group, file1, file2) }
    coverage_ch = coverage_ch.groupTuple()

    // collect all mapping files, provide directory for clustering
    mapping_files_info = mapped_reads_ch.map { file, group, _refname -> tuple(group, file) }
    mapping_files_info = mapping_files_info.groupTuple()

    // Merge mapping statistics
    merged_coverage_ch = MergeCoverageStatistics(coverage_ch)

    // Cluster mapped reads across alignment files
    clustering_ch = ClusterMappedReads(analysis_id, mapping_files_info)

    clustering_ch.clade_report.ifEmpty {
        log.info("No clustering results generated. Ending workflow.")
        System.exit(0)
    }

    MatchCladeReportWithReferenceSequences(
        analysis_id,
        clustering_ch.clade_report,
        reference_sequences_ch.matched_assemblies,
        merged_coverage_ch.merged_coverage_statistics,
        merge_classification_results_ch,
    )
}


/*
* Check reads_ch, if tuple contains only one file, append a placeholder NO_FILE for the second file
*/
process CheckReads {
    input:
    tuple path(fastq1), path(fastq2)

    output:
    tuple path(fastq1), path(fastq2)

    script:
    """
    if [ ! -f ${fastq2} ]; then
        echo "No second read file found, using placeholder"
        fastq2="NO_FILE"
    fi
    """
}


/*
* Quality control of single end reads using prinseq++
*/
process QCReadsPrinseqSingle {
    input:
    val analysis_id
    tuple path(fastq1)

    output:
    tuple path("${analysis_id}_good.fastq.gz")

    script:
    """
    prinseq++ ${params.prinseq_params} -fastq ${fastq1} -out_good ${analysis_id}_good.fastq -out_bad ${analysis_id}_bad.fastq
    bgzip ${analysis_id}_good.fastq && bgzip ${analysis_id}_bad.fastq
    """
}

/*
* Quality control of paired reads using prinseq++
*/
process QCReadsPrinseqPaired {
    debug true

    input:
    val analysis_id
    tuple path(fastq1), path(fastq2)

    output:
    path "${analysis_id}_good_R1.fastq.gz", emit: _good_R1
    path "${analysis_id}_good_R2.fastq.gz", emit: _good_R2, optional: true

    script:
    """
    if [ -f ${fastq2} ]; then
        prinseq++ ${params.prinseq_params} -fastq ${fastq1} -fastq2 ${fastq2} -out_good ${analysis_id}_good_R1.fastq -out_bad ${analysis_id}_bad_R1.fastq \
        -out_good2 ${analysis_id}_good_R2.fastq -out_bad2 ${analysis_id}_bad_R2.fastq
        bgzip ${analysis_id}_good_R1.fastq && bgzip ${analysis_id}_good_R2.fastq
        bgzip ${analysis_id}_bad_R1.fastq && bgzip ${analysis_id}_bad_R2.fastq
    else
        prinseq++ ${params.prinseq_params} -fastq ${fastq1} -out_good ${analysis_id}_good_R1.fastq -out_bad ${analysis_id}_bad_R1.fastq
        bgzip ${analysis_id}_good_R1.fastq && bgzip ${analysis_id}_bad_R1.fastq
    fi
    """
}


/*
* merge coverage statistics from different BAM files
*/
process MergeCoverageStatistics {
    tag "MergeCoverageStatistics ${analysis_id}"

    input:
    tuple val(analysis_id), path(stats_files), path(coverage_files)

    output:
    path "merged_coverage_statistics.tsv", emit: merged_coverage_statistics

    script:
    def coverage_files_string = coverage_files.collect { it[0] }.join(',')
    def stats_files_string = stats_files.collect { it[0] }.join(',')
    """
    #!/usr/bin/env python3
    import pandas as pd
    import glob
    
    
    cov_files = "${coverage_files_string}".split(',')
    stats_files = "${stats_files_string}".split(',')
    coverage_data = []
    for ix, file in enumerate(cov_files):
        stats_filename = stats_files[ix]
        stats_df = pd.read_csv(stats_filename, sep="\\t", header=None, names=["stat", "value", "comment"])
        df = pd.read_csv(file, sep="\\t")
        df['file'] = file.split('/')[-1]
        df['error_rate'] = stats_df.loc[stats_df['stat'] == 'error rate:', 'value'].values[0]
        coverage_data.append(df)
    merged_df = pd.concat(coverage_data, ignore_index=True)
    

    merged_df.reset_index(inplace=True)
    merged_df.to_csv("merged_coverage_statistics.tsv", sep="\\t", index=False)
    """
}

/*
* Use samtools to extract BAM file statistics
*/
process SamtoolsStats {
    tag "SamtoolsStats ${bamfile.baseName}"

    input:
    tuple path(bamfile), path(bamindex), val(query_id), val(reference_id)

    output:
    tuple path("${bamfile.baseName}.stats.txt"), val(bamfile.baseName), val(reference_id), val(query_id)

    script:
    """
    samtools stats ${bamfile} > ${bamfile.baseName}.stats.txt
    """
}


/*
* Use samtools to extract BAM file coverage statistics
*/
process SamtoolsCoverage {

    tag "SamtoolsCoverage ${bamfile.baseName}"

    input:
    tuple path(bamfile), path(bamindex), val(query_id), val(reference_id)

    output:
    tuple path("${bamfile.baseName}.stats.txt"), path("${bamfile.baseName}.coverage.txt"), val(bamfile.baseName), val(reference_id), val(query_id)

    script:
    """
    samtools coverage -o ${bamfile.baseName}.coverage.txt ${bamfile}
    samtools stats ${bamfile} | grep ^SN | cut -f 2- > ${bamfile.baseName}.stats.txt
    """
}


/*
* sort and index bam file, maintain tuple file, query_id, reference_id in channel
*/
process sortBam {
    tag "sortMapping"

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
* Filter bam file using msamtools
*/
process FilterBamMsamtools {
    tag "FilterBamMsamtools ${bamfile.baseName}"

    input:
    tuple path(bamfile), val(query_id), val(reference_id)

    output:
    tuple path("${query_id}_${reference_id}.filtered.bam"), val(query_id), val(reference_id)

    script:
    """
    msamtools filter -b ${params.msamtools_params} ${bamfile} > ${query_id}_${reference_id}.filtered.bam
    """
}


/*
* match clustering clade_report.txt and reference_sequences/matched_assemblies.tsv
*/
process MatchCladeReportWithReferenceSequences {
    tag "MatchCladeReportWithReferenceSequences ${analysis_id}"
    publishDir "${params.output_dir}/${analysis_id}/output", mode: 'copy'

    input:
    val analysis_id
    path clade_report
    path matched_assemblies
    path coverage_report
    path merge_classification_results

    output:
    path "clade_report_with_references.tsv", emit: clade_report_with_references
    path clade_report
    path matched_assemblies
    path coverage_report

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    clade_report = pd.read_csv("${clade_report}", sep="\\t", header=None,  names=["clade", "nuniq", "freq", "min_pair_dist", "nfiles", "files"])
    clade_report['clade']
    clade_report['files'] = clade_report['files'].str.split(',')
    clade_report = clade_report.explode('files')

    matched_assemblies = pd.read_csv("${matched_assemblies}", sep="\\t")
    matched_assemblies['filename'] = matched_assemblies['assembly_file'].str.split('/').str[-1]

    coverage_report = pd.read_csv("${coverage_report}", sep="\\t")
    merged_classification_results = pd.read_csv("${merge_classification_results}", sep="\\t")

    def find_assembly_mapping(row):
        accession = row['assembly_accession']
        if accession is None or pd.isna(accession):
            row['clade'] = 'unmapped'
            row['nuniq'] = 0
            row['freq'] = 0
            return row
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
        if accession is None or pd.isna(accession):
            row['coverage'] = 0
            return row
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
    tag "ClusterMappedReads ${analysis_id}"
    publishDir "${params.output_dir}/${analysis_id}", mode: 'copy'

    input:
    val analysis_id
    tuple val(query_id), path(mapped_reads)

    output:
    path "clustering/clade_report.tsv", emit: clade_report, optional: true
    path "clustering/sample_report.tsv", emit: sample_report, optional: true
    path "clustering/distance_matrix.tsv", emit: distance_matrix, optional: true
    path "clustering/all_node_statistics.tsv", emit: all_node_statistics, optional: true
    path "clustering/nj_tree_edges.txt", emit: nj_tree_edges, optional: true

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

    input:
    val analysis_id
    tuple path(fastq1), path(fastq2), path(reference)
    val minimap2_params

    output:
    tuple path("${fastq1.baseName}_${reference.baseName}.bam"), val(analysis_id), val(reference.baseName)

    script:
    """
    if [ -f ${fastq2} ]; then
        mkdir -p ${analysis_id}
        minimap2 ${minimap2_params} -ax sr ${reference} ${fastq1} ${fastq2} | samtools view -bS -F 4 - > ${fastq1.baseName}_${reference.baseName}.bam
    else
        minimap2 ${minimap2_params} -ax sr ${reference} ${fastq1} | samtools view -bS -F 4 - > ${fastq1.baseName}_${reference.baseName}.bam
    fi
    """
}

/*
* Extract reference sequences from classifier output results to global reference database
*/
process ExtractReferenceSequences {
    tag "ExtractReferenceSequences ${analysis_id}"

    input:
    val analysis_id
    path classifier_output

    output:
    path "reference_sequences/*gz", emit: reference_sequences, optional: true
    path "reference_sequences/matched_assemblies.tsv", emit: matched_assemblies

    script:
    """
    ${params.python_bin} ${params.references_extract_script} retrieve \
    --input_table ${classifier_output} \
    --assembly_store "${params.assembly_store}" \
    --mapping_references_dir "reference_sequences" \
    --include_term "complete" \
    --exclude_term "plasmid"
    """
}


/*
* Merge classification results from different classifiers
*/
process MergeClassificationResults {
    tag "MergeClassificationResults ${query_id}"
    publishDir "${params.output_dir}/${query_id}/classification", mode: 'copy'

    input:
    path centrifuge_output
    path centrifuge_output_processed
    val query_id
    path kraken2_output
    path kraken2_output_processed
    val query_id_2
    path krakenunique_output
    path krakenunique_output_processed
    val query_id_3
    path diamond_output
    path diamond_output_processed
    val query_id_4

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
    krakenunique_df = pd.read_csv("${krakenunique_output_processed}", sep="\\t").rename(columns={"taxID": "taxid", "name": "description"})
    diamond_df = pd.read_csv("${diamond_output_processed}", sep="\\t").rename(columns={"taxID": "taxid", "name": "description"})
    all_taxids = set(centrifuge_df['taxid']).union(set(kraken2_df['taxid'])).union(set(krakenunique_df['taxid'])).union(set(diamond_df['taxid']))

    merged_df = pd.DataFrame({'taxid': list(all_taxids)})
    def get_classification(row):
        taxid = row['taxid']
        classifiers = []
        uniq_reads = 0
        description = None
        if taxid in kraken2_df['taxid'].values:
            classifiers.append('kraken2')
            uniq_reads += kraken2_df[kraken2_df['taxid'] == taxid]['uniq_reads'].sum()
            description = kraken2_df[kraken2_df['taxid'] == taxid]['description'].values[0]
        if taxid in centrifuge_df['taxid'].values:
            classifiers.append('centrifuge')
            uniq_reads += centrifuge_df[centrifuge_df['taxid'] == taxid]['uniq_reads'].sum()
            if description is None:
                description = centrifuge_df[centrifuge_df['taxid'] == taxid]['description'].values[0]
        if taxid in krakenunique_df['taxid'].values:
            classifiers.append('krakenunique')
            uniq_reads += krakenunique_df[krakenunique_df['taxid'] == taxid]['uniq_reads'].sum()
            if description is None:
                description = krakenunique_df[krakenunique_df['taxid'] == taxid]['description'].values[0]
        if taxid in diamond_df['taxid'].values:
            classifiers.append('diamond')
            uniq_reads += diamond_df[diamond_df['taxid'] == taxid]['uniq_reads'].sum()
            if description is None:
                description = diamond_df[diamond_df['taxid'] == taxid]['description'].values[0]
        classifiers = '/'.join(classifiers) if classifiers else 'unclassified'
        row['classification'] = classifiers
        row['total_uniq_reads'] = uniq_reads
        row['description'] = description
        return row

    merged_df = merged_df.apply(get_classification, axis=1)
    merged_df.to_csv(f"${query_id}_merged_classification.tsv", sep="\\t", index=False)
    """
}

/*
 * Run Kraken2 classification on paired-end reads.
 */
process Kraken2ClassificationPaired {
    tag "Kraken2Classification ${analysis_id}"
    publishDir "${params.output_dir}/${analysis_id}/classification/kraken2", mode: 'copy'

    input:
    val analysis_id
    tuple path(fastq1), path(fastq2)

    output:
    path "${analysis_id}_kraken2_report.txt"
    path "${analysis_id}_krk2_processed_classifier_output.tsv"
    val analysis_id

    script:
    """

    if [ -f ${fastq2} ]; then
        ${params.kraken2_bin} --db ${params.kraken2_index} \
        --report ${analysis_id}_kraken2_report.txt \
        --output ${analysis_id}_kraken2_classification.txt \
        ${fastq1} ${fastq2} ${params.kraken2_params}
    else
        ${params.kraken2_bin} --db ${params.kraken2_index} \
        --report ${analysis_id}_kraken2_report.txt \
        --output ${analysis_id}_kraken2_classification.txt \
        ${fastq1} ${params.kraken2_params}
    fi

    ${params.python_bin} ${params.classifier_process_script} \
    --input "${analysis_id}_kraken2_report.txt" \
    --output "${analysis_id}_krk2_processed_classifier_output.tsv" \
    --type "kraken2" \
    --nuniq_threshold ${params.minimum_uniq_reads}

    """
}


/*
* Classifyt reads using Diamond in paired-end mode
*/
process DiamondClassificationPaired {
    tag "DiamondClassificationSingle ${analysis_id}"
    publishDir "${params.output_dir}/${analysis_id}/classification/diamond", mode: 'copy'

    input:
    val analysis_id
    tuple path(fastq1), path(fastq2)

    output:
    path "${analysis_id}_diamond_classification.tsv"
    path "${analysis_id}_diamond_processed_classifier_output.tsv"
    val "${analysis_id}"

    script:
    """

    if [ -f ${fastq2} ]; then 
        ${params.diamond_bin} blastx \
        --db ${params.diamond_index} \
        --query ${fastq1} \
        --query ${fastq2} \
        --out ${analysis_id}_diamond_classification.tsv \
        --outfmt 6 \
        ${params.diamond_params}
    else
        ${params.diamond_bin} blastx \
        --db ${params.diamond_index} \
        --query ${fastq1} \
        --out ${analysis_id}_diamond_classification.tsv \
        --outfmt 6 \
        ${params.diamond_params}
    fi

    ${params.python_bin} ${params.classifier_process_script} \
    --input "${analysis_id}_diamond_classification.tsv" \
    --output "${analysis_id}_diamond_processed_classifier_output.tsv" \
    --type "diamond" \
    --nuniq_threshold ${params.minimum_uniq_reads}
    """
}

/*
* Classify reads using KrakenUnique in paired-end mode
*/
process KrakenUniqueClassificationPaired {
    tag "KrakenUniqueClassificationSingle ${analysis_id}"
    publishDir "${params.output_dir}/${analysis_id}/classification/krakenunique", mode: 'copy'

    input:
    val analysis_id
    tuple path(fastq1), path(fastq2)

    output:
    path "${analysis_id}_krakenunique_classification.txt"
    path "${analysis_id}_krakenunique_processed_classifier_output.tsv"
    val "${analysis_id}"

    script:
    """
    if [ -f ${fastq2} ]; then
        ${params.krakenunique_bin} --db ${params.krakenunique_index} \
        --report ${analysis_id}_krakenunique_report.txt \
        --output ${analysis_id}_krakenunique_classification.txt \
        ${fastq1} ${fastq2} ${params.krakenunique_params}
    else
        ${params.krakenunique_bin} --db ${params.krakenunique_index} \
        --report ${analysis_id}_krakenunique_report.txt \
        --output ${analysis_id}_krakenunique_classification.txt \
        ${fastq1} ${params.krakenunique_params}
    fi

    ${params.python_bin} ${params.classifier_process_script} \
    --input "${analysis_id}_krakenunique_classification.txt" \
    --output "${analysis_id}_krakenunique_processed_classifier_output.tsv" \
    --type "kuniq" \
    --nuniq_threshold ${params.minimum_uniq_reads}

    """
}



/*
* Classify reads using Centrifuge in paired-end mode
*/
process CentrifugeClassificationPaired {
    tag "CentrifugeClassificationSingle ${analysis_id}"
    publishDir "${params.output_dir}/${analysis_id}/classification/centrifuge", mode: 'copy'

    input:
    val analysis_id
    tuple path(fastq1), path(fastq2)

    output:
    path "${analysis_id}_centrifuge_report.txt"
    path "${analysis_id}_centrifuge_processed_classifier_output.tsv"
    val "${analysis_id}"

    script:
    """
    if [ -f ${fastq1} ] && [ -f ${fastq2} ]; then
        centrifuge -x ${params.centrifuge_index} -1 ${fastq1} -2 ${fastq2} \
        -S ${analysis_id}_centrifuge_classification.tsv \
        --output ${analysis_id}_centrifuge_classification.txt \
        --report-file ${analysis_id}_centrifuge_report.txt \
        ${params.centrifuge_params}
    elif [ -f ${fastq1} ]; then
        centrifuge -x ${params.centrifuge_index} -U ${fastq1} \
        -S ${analysis_id}_centrifuge_classification.tsv \
        --output ${analysis_id}_centrifuge_classification.txt \
        --report-file ${analysis_id}_centrifuge_report.txt \
        ${params.centrifuge_params}
    fi

    ${params.python_bin} ${params.classifier_process_script} \
    --input "${analysis_id}_centrifuge_report.txt" \
    --output "${analysis_id}_centrifuge_processed_classifier_output.tsv" \
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
