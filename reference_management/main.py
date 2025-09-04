from utils.ncbi_tools import Passport, retrieve_sequence_databases, query_sequence_databases, LocalAssembly
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


    def retrieve_local_assembly(self, taxid) -> Optional[LocalAssembly]:
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

        if not assembly_file:
            print(f"No assembly file found for taxid {taxid}")
            return None

        accid = "_".join(taxid_files[0].split('_')[:2]) if assembly_file else None

        return LocalAssembly(taxid=taxid, accession=accid, file_path=os.path.join(taxid_subdir, assembly_file[0])) if assembly_file else None



    def retrieve_assembly(self, passport: Passport) -> Optional[LocalAssembly]:
        """
        Retrieve the assembly for the given taxid, either from local storage or NCBI.
        """
        # First, check if the assembly is available locally
        local_assembly = self.retrieve_local_assembly(passport.taxid)
        if local_assembly:
            print(f"Using local assembly for taxid {passport.taxid}: {local_assembly.file_path}")
            return local_assembly
        
        # If not found locally, fetch from NCBI
        print(f"Fetching assembly for taxid {passport.taxid} from NCBI...")

        reference_data = query_sequence_databases(passport)
        assembly_dir = os.path.join(self.store_path, passport.prefix)
        os.makedirs(assembly_dir, exist_ok=True)
        
        assembly_file_path = os.path.join(assembly_dir, f"{passport.prefix}_sequence.fasta.gz")
        success_dl = retrieve_sequence_databases(reference_data, assembly_file_path, gzipped=True)


        if not success_dl:
            print(f"Failed to download assembly for passport {passport.taxid}")
            return None
        
        return LocalAssembly(taxid=passport.taxid, accession=reference_data.accession, file_path=assembly_file_path)

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
            passport = Passport(taxid = taxid, accession = None)
            local_assembly = self.retrieve_assembly(passport)
            if local_assembly:
                df.at[index, 'assembly_accession'] = local_assembly.accession
                df.at[index, 'assembly_file'] = local_assembly.file_path

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

        reference_data = query_sequence_databases(passport)

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