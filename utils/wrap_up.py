#!/usr/bin/env python3
import pandas as pd

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="Match clade report with reference sequences")
    parser.add_argument("output_dir", type=str, help="workflow output directory with clade report and matched assemblies")
    parser.add_argument("--clade_report", type=str, default="clade_report.tsv", help="clade report file name")
    parser.add_argument("--matched_assemblies", type=str, default="matched_assemblies.tsv", help="matched assemblies file name")
    parser.add_argument("--coverage_report", type=str, default="merged_coverage_statistics.tsv", help="coverage report file name")
    args = parser.parse_args()
    return args

def main():

    args = parse_args()
    clade_report = f"{args.output_dir}/{args.clade_report}"
    matched_assemblies = f"{args.output_dir}/{args.matched_assemblies}"
    coverage_report = f"{args.output_dir}/{args.coverage_report}"

    clade_report = pd.read_csv(clade_report, sep="\t", header=None, names=["clade", "nuniq", "freq", "files"])
    clade_report['files'] = clade_report['files'].str.split(',')
    clade_report = clade_report.explode('files')

    matched_assemblies = pd.read_csv(matched_assemblies, sep="\t")
    matched_assemblies['filename'] = matched_assemblies['assembly_file'].str.split('/').str[-1]

    coverage_report = pd.read_csv(coverage_report, sep="\t")

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

    clade_report_with_references = matched_assemblies.apply(find_assembly_mapping, axis=1)
    clade_report_with_references = clade_report_with_references.apply(find_assembly_coverage, axis=1)
    clade_report_with_references = clade_report_with_references[['description', 'taxid', 'assembly_accession', \
            'coverage', 'clade', 'nuniq', 'freq']]

    clade_report_with_references.to_csv(f"{args.output_dir}/clade_report_with_references.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()