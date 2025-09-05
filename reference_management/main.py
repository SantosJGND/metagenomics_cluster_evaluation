from utils.ncbi_tools import Passport, NCBITools
from utils.reference_utils import AssemblyStore
import os
import pandas as pd

import argparse

def get_args():
    """
    Define the argument parser with subcommands.
    """
    parser = argparse.ArgumentParser(description="Manage taxid-to-assembly retrieval and reference setup.")
    subparsers = parser.add_subparsers(dest="command", required=True, help="Subcommands: retrieve or check")

    # Subcommand: retrieve
    retrieve_parser = subparsers.add_parser("retrieve", help="Retrieve assemblies based on the input table.")
    retrieve_parser.add_argument(
        "--input_table",
        type=str,
        required=True,
        help="Path to the classification output file."
    )
    retrieve_parser.add_argument(
        "--assembly_store",
        type=str,
        default="assemblies",
        help="Directory to store downloaded assemblies."
    )
    retrieve_parser.add_argument(
        "--mapping_references_dir",
        type=str,
        default="references_to_map",
        help="Directory to store mapping references."
    )

    # Subcommand: check
    check_parser = subparsers.add_parser("check", help="Check if mapping ids can be retrieved.")
    check_parser.add_argument(
        "--input_table",
        type=str,
        required=True,
        help="Path to the classification output file."
    )

    check_parser.add_argument(
        "--assessment",
        type=str,
        default="assembly_assessment.tsv",
        help="Path to the assessment file to check assemblies."
    )

    return parser.parse_args()


def retrieve_assemblies(args):
    """
    Retrieve assemblies based on the input table and store them in the specified directory.
    """
    classification_output_path = args.input_table
    assembly_store = args.assembly_store
    mapping_references_dir = args.mapping_references_dir

    assembly_store = AssemblyStore(assembly_store)
    df = assembly_store.match_taxid_to_assembly(classification_output_path)
    assembly_store.setup_mapping_references(df, mapping_references_dir=mapping_references_dir)
    df.to_csv(
        os.path.join(mapping_references_dir, "matched_assemblies.tsv"),
        index=False,
        sep='\t'
    )

def check_assemblies_exist(args):
    """
    Check if the assemblies can be retrieved from the input table.
    """

    input_table = args.input_table
    df = pd.read_csv(input_table, sep='\t')
    taxid_col = False
    if 'taxid' in df.columns:
        df.rename(columns={'taxid': 'taxid'}, inplace=True)
        taxid_col = True
    elif 'TaxID' in df.columns:
        df.rename(columns={'TaxID': 'taxid'}, inplace=True)
        taxid_col = True
    elif 'taxon' in df.columns:
        df.rename(columns={'taxon': 'taxid'}, inplace=True)
        taxid_col = True


    accid_col = False
    if 'assembly_accession' in df.columns:
        df.rename(columns={'assembly_accession': 'accid'}, inplace=True)
        accid_col = True
    elif "accession" in df.columns:
        df.rename(columns={'accession': 'accid'}, inplace=True)
        accid_col = True    
    elif "accID" in df.columns:
        df.rename(columns={'accID': 'accid'}, inplace=True)
        accid_col = True
    elif "accid" in df.columns:
        df.rename(columns={'accid': 'accid'}, inplace=True)
        accid_col = True
    

    
    if taxid_col is False and accid_col is False:
        raise ValueError("The classification output file must contain a taxonomic ID column [taxid, taxID or taxon] or an accession column [assembly_accession, accession, accID or accid].")

    def check_assembly_exists(row, taxid_col= False, accid_col=False):
        """
        Check if the assembly for the given taxid exists.
        """
        taxid = None
        accid = None

        if taxid_col:
            taxid = str(int(row["taxid"]))
        if accid_col:
            accid = str(row["accid"])

        passport = Passport(taxid = taxid, accession = accid)
        print(f"Processing passport: {passport}")
        ncbi_tools = NCBITools()
        reference_data = ncbi_tools.query_sequence_databases(passport)

        print("retrieved reference data:", reference_data)

        row['assembly_accession'] = reference_data.accession
        row['description'] = reference_data.description
        row['nucleotide_id'] = reference_data.nucleotide_id
        row['assembly_id'] = reference_data.assembly_id

        return row

    df = df.apply(lambda row: check_assembly_exists(row, taxid_col=taxid_col, accid_col=accid_col), axis=1)
    df.to_csv(args.assessment, index=False, sep='\t')

def main():

    args = get_args()

    if args.command == "retrieve":
        retrieve_assemblies(args)
    elif args.command == "check":
        check_assemblies_exist(args)
    


if __name__ == "__main__":
    main()