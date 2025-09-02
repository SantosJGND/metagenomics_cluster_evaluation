from utils.ncbi_tools import get_representative_assembly, retrieve_sequence_databases, query_sequence_databases
import os
from typing import Optional, Tuple
import pandas as pd


def dl_file(url: str, dest: str) -> None:
    """
    Download a file from a given URL to a specified destination.
    """
    import requests
    response = requests.get(url)
    if response.status_code == 200:
        with open(dest, 'wb') as f:
            f.write(response.content)
        print(f"Downloaded {url} to {dest}")
        return True
    else:
        print(f"Failed to download {url}: {response.status_code}")
        return False 

def dl_file_wget(url: str, dest: str) -> None:
    """
    Download a file using wget.
    """
    import os
    command = f"wget -O {dest} {url}"
    os.system(command)
    print(f"Downloaded {url} to {dest}")

    if not os.path.exists(dest):
        print(f"Failed to download {url} to {dest}")
        return False
    return True


class AssemblyStore:
    """
    Class to manage assembly storage and retrieval.
    """
    def __init__(self, store_path: str):
        self.store_path = store_path
        os.makedirs(self.store_path, exist_ok=True)

    def get_assembly_path(self, taxid: str) -> str:
        return os.path.join(self.store_path, taxid)


    def retrieve_local_assembly(self, taxid) -> Optional[Tuple[str, str]]:
        """
        Check if the assembly for the given taxid exists in the assembly store.
        """
        taxid_subdir = f"{self.store_path}/{taxid}"
        if os.path.exists(taxid_subdir) is False:
            print(f"No assembly directory found for taxid {taxid}")
            return None 
        taxid_files = os.listdir(taxid_subdir)
        if not taxid_files:
            print(f"No files found in assembly directory for taxid {taxid}")
            return None
        
        # Assuming the first file is the assembly file
        assembly_file = [x for x in taxid_files if x.endswith('.gz')]

        accid = "_".join(taxid_files[0].split('_')[:2]) if assembly_file else None

        return accid, os.path.join(taxid_subdir, assembly_file[0]) if assembly_file else None


    def retrieve_assembly(self, taxid: str) -> Optional[Tuple[str, str]]:
        """
        Retrieve the assembly for the given taxid, either from local storage or NCBI.
        """
        # First, check if the assembly is available locally
        local_assembly = self.retrieve_local_assembly(taxid)
        if local_assembly:
            print(f"Using local assembly for taxid {taxid}: {local_assembly[0]}")
            return local_assembly
        
        # If not found locally, fetch from NCBI
        print(f"Fetching assembly for taxid {taxid} from NCBI...")
        reference_data = query_sequence_databases(taxid)
        assembly_dir = os.path.join(self.store_path, taxid)
        os.makedirs(assembly_dir, exist_ok=True)
        
        output_path = os.path.join(assembly_dir, f"{taxid}_sequence.fasta.gz")

        success_dl, assembly_file_path = retrieve_sequence_databases(reference_data, output_path, gzipped=True)


        if not success_dl:
            print(f"Failed to download assembly for taxid {taxid}")
            return None
        
        return reference_data.accession, assembly_file_path

    def match_taxid_to_assembly(
        self,
        classification_output_path: str,
    ) -> pd.DataFrame:
        """
        Match taxids from the classification output to their respective assemblies.
        """
        if not os.path.exists(classification_output_path):
            raise FileNotFoundError(f"Classification output file not found: {classification_output_path}")

        df = pd.read_csv(classification_output_path, sep='\t')
        taxid_col = None
        if 'taxid' in df.columns:
            taxid_col = 'taxid'
        elif 'TaxID' in df.columns:
            taxid_col = 'TaxID'
        elif 'taxon' in df.columns:
            taxid_col = 'taxon'
        else:
            raise ValueError("The classification output file must contain a taxonomic ID column [taxid, taxID or taxon].")

        df['assembly_accession'] = None
        df['assembly_file'] = None

        for index, row in df.iterrows():
            taxid = str(int(row[taxid_col]))
            print(f"Processing taxid {taxid}...")
            accession, assembly_file = self.retrieve_assembly(taxid)
            if accession and assembly_file:
                df.at[index, 'assembly_accession'] = accession
                df.at[index, 'assembly_file'] = assembly_file

        return df

    @staticmethod
    def setup_mapping_references(
        classification_output_path: pd.DataFrame,
        mapping_references_dir: str = "references_to_map"):

        if not os.path.exists(mapping_references_dir):
            os.makedirs(mapping_references_dir)
        
        if "assembly_accession" not in classification_output_path.columns or \
        "assembly_file" not in classification_output_path.columns:
            raise ValueError("The DataFrame must contain 'assembly_accession' and 'assembly_file' columns.")
        
        for index, row in classification_output_path.iterrows():
            accession = row['assembly_accession']
            assembly_file = row['assembly_file']
            
            if pd.isna(accession) or pd.isna(assembly_file):
                print(f"Skipping taxid {row['taxid']} due to missing assembly data.")
                continue
            
            dest_path = os.path.join(mapping_references_dir, f"{accession}.fna.gz")
            if not os.path.exists(dest_path):
                print(f"Copying {assembly_file} to {dest_path}")
                os.system(f"cp {assembly_file} {dest_path}")
            else:
                print(f"File {dest_path} already exists, skipping copy.")
        
        print(f"Mapping references setup complete in {mapping_references_dir}.")



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
    taxid_col = None
    if 'taxid' in df.columns:
        taxid_col = 'taxid'
    elif 'TaxID' in df.columns:
        taxid_col = 'TaxID'
    elif 'taxon' in df.columns:
        taxid_col = 'taxon'
    else:
        raise ValueError("The classification output file must contain a taxonomic ID column [taxid, taxID or taxon].")

    def check_assembly_exists(row, taxid_col):
        """
        Check if the assembly for the given taxid exists.
        """
        taxid = str(int(row[taxid_col]))

        reference_data = query_sequence_databases(taxid)

        row['assembly_accession'] = reference_data.accession
        row['description'] = reference_data.description
        row['nucleotide_id'] = reference_data.nucleotide_id
        row['assembly_id'] = reference_data.assembly_id

        return row
        
    df = df.apply(lambda row: check_assembly_exists(row, taxid_col), axis=1)
    df.to_csv(args.assessment, index=False, sep='\t')

def main():

    args = get_args()

    if args.command == "retrieve":
        retrieve_assemblies(args)
    elif args.command == "check":
        check_assemblies_exist(args)
    


if __name__ == "__main__":
    main()