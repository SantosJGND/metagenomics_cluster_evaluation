
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


    def retrieve_local_assembly(self, taxid) -> Optional[LocalAssembly]:
        """
        Check if the assembly for the given taxid exists in the assembly store.
        """
        taxid_subdir = f"{self.store_path}/{taxid}"
        if os.path.exists(taxid_subdir) is False:
            self.logger.warning(f"No assembly directory found for taxid {taxid}")
            return None
        taxid_files = os.listdir(taxid_subdir)
        if not taxid_files:
            self.logger.warning(f"No files found in assembly directory for taxid {taxid}")
            return None
        
        # Assuming the first file is the assembly file
        assembly_file = [x for x in taxid_files if x.endswith('.gz')]

        if not assembly_file:
            self.logger.warning(f"No assembly file found for taxid {taxid}")
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
            self.logger.info(f"Using local assembly for taxid {passport.taxid}: {local_assembly.file_path}")
            return local_assembly
        
        # If not found locally, fetch from NCBI
        self.logger.info(f"Fetching assembly for taxid {passport.taxid} from NCBI...")

        reference_data = self.ncbi.query_sequence_databases(passport)
        assembly_dir = os.path.join(self.store_path, passport.prefix)
        os.makedirs(assembly_dir, exist_ok=True)
        
        assembly_file_path = os.path.join(assembly_dir, f"{passport.prefix}_sequence.fasta.gz")
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
            self.logger.info(f"Processing taxid {taxid}...")
            passport = Passport(taxid = taxid, accession = None)
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
            
            dest_path = os.path.join(mapping_references_dir, f"{accession}.fna.gz")
            if not os.path.exists(dest_path):
                self.logger.info(f"Copying {assembly_file} to {dest_path}")
                os.system(f"cp {assembly_file} {dest_path}")
            else:
                self.logger.warning(f"File {dest_path} already exists, skipping copy.")

        self.logger.info(f"Mapping references setup complete in {mapping_references_dir}.")


