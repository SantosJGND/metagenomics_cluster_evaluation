
import logging
import os
import pandas as pd
from typing import Optional
from utils.ncbi_tools import Passport, NCBITools, LocalAssembly

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
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.DEBUG)
        handler = logging.StreamHandler()
        handler.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        handler.setFormatter(formatter)
        self.logger.addHandler(handler)
        self.logger.propagate = False

        self.ncbi = NCBITools()

    def get_assembly_path(self, taxid: str) -> str:
        return os.path.join(self.store_path, taxid)


    def retrieve_local_assembly(self, passport: Passport) -> Optional[LocalAssembly]:
        """
        Check if the assembly for the given taxid exists in the assembly store.
        """
        taxid_subdir = f"{self.store_path}/{passport.taxid}"

        # Assuming the first file is the assembly file
        assembly_file = os.path.join(taxid_subdir, f"{passport.prefix}_sequence.fasta.gz")
        print(passport.prefix)
        print(f"Checking for local assembly at {assembly_file}")

        if not os.path.exists(assembly_file):
            self.logger.warning(f"No assembly file found for taxid {passport.taxid} and accession {passport.accession}")
            return None
        accid = passport.accession

        return LocalAssembly(taxid=passport.taxid, accession=accid, file_path=assembly_file) if assembly_file else None


    def retrieve_assembly(self, passport: Passport) -> Optional[LocalAssembly]:
        """
        Retrieve the assembly for the given taxid, either from local storage or NCBI.
        """
        # First, check if the assembly is available locally
        local_assembly = self.retrieve_local_assembly(passport)

        print("### local assembly ###")
        print("local_assembly", local_assembly)
        if local_assembly:
            self.logger.info(f"Using local assembly for taxid {passport.taxid}: {local_assembly.file_path}")
            return local_assembly
        
        # If not found locally, fetch from NCBI
        self.logger.info(f"Fetching assembly for taxid {passport.taxid} from NCBI...")

        reference_data = self.ncbi.query_sequence_databases(passport)
        assembly_dir = os.path.join(self.store_path, str(passport.taxid))
        os.makedirs(assembly_dir, exist_ok=True)

        assembly_file_path = os.path.join(assembly_dir, f"{reference_data.prefix}_sequence.fasta.gz")
        success_dl = self.ncbi.retrieve_sequence_databases(reference_data, assembly_file_path, gzipped=True)

        if not success_dl:
            self.logger.error(f"Failed to download assembly for passport {passport.taxid}")
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
        if taxid_col is False:
            raise ValueError("The classification output file must contain a taxonomic ID column [taxid, taxID or taxon].")

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

        for index, row in df.iterrows():
            taxid = str(int(row['taxid'])) if taxid_col == True else None
            accession = str(row['accid']) if accid_col == True else None
            if taxid is None and accession is None:
                self.logger.warning(f"Skipping row {index} due to missing taxid and accession.")
                continue
            self.logger.info(f"Processing taxid {taxid}...")
            passport = Passport(taxid = taxid, accession = accession)
            local_assembly = self.retrieve_assembly(passport)

            if local_assembly:
                df.at[index, 'assembly_accession'] = local_assembly.accession
                df.at[index, 'assembly_file'] = local_assembly.file_path

        return df

    
    def setup_mapping_references(
        self,
        classification_output_path: pd.DataFrame,
        mapping_references_dir: str = "references_to_map"):

        if not os.path.exists(mapping_references_dir):
            os.makedirs(mapping_references_dir)
        
        if "assembly_accession" not in classification_output_path.columns or \
        "assembly_file" not in classification_output_path.columns:
            raise ValueError("The DataFrame must contain 'assembly_accession' and 'assembly_file' columns.")
        
        for _, row in classification_output_path.iterrows():
            accession = row['assembly_accession']
            assembly_file = row['assembly_file']

            if pd.isna(accession) or pd.isna(assembly_file):
                self.logger.warning(f"Skipping taxid {row['taxid']} due to missing assembly data.")
                continue
            
            dest_filename = f"{accession}.fna.gz"
            if "taxid" in row and not pd.isna(row["taxid"]):
                dest_filename = f"{row['taxid']}_{dest_filename}"
            dest_path = os.path.join(mapping_references_dir, dest_filename)
            
            if not os.path.exists(dest_path):
                self.logger.info(f"Copying {assembly_file} to {dest_path}")
                os.system(f"cp {assembly_file} {dest_path}")
            else:
                self.logger.warning(f"File {dest_path} already exists, skipping copy.")

        self.logger.info(f"Mapping references setup complete in {mapping_references_dir}.")


