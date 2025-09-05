# You need to install Biopython: pip install biopython
import logging
from Bio import Entrez
import os
import subprocess
from dataclasses import dataclass
from typing import Tuple, Optional
# read email from local .env file

import dotenv
dotenv.load_dotenv()

Entrez.email = os.getenv("NCBI_EMAIL", None)
if Entrez.email is None:
    raise ValueError("NCBI_EMAIL environment variable not set. Please set it to your email address.")

@dataclass
class Passport:
    taxid: Optional[str]
    accession: Optional[str]



    def __str__(self):
        return f"TaxID: {self.taxid}, Accession: {self.accession}"

    @property
    def prefix(self):
        if self.accession:
            return f"{self.taxid}_{self.accession}"
        else:
            return f"{self.taxid}"


@dataclass
class LocalAssembly(Passport):

    file_path: str

    def __str__(self):
        return f"TaxID: {self.taxid}, Accession: {self.accession}, File Path: {self.file_path}"

    
@dataclass
class ReferenceData(Passport):

    description: Optional[str] = None
    nucleotide_id: Optional[str] = None
    assembly_id: Optional[str] = None

    def __str__(self):
        return f"TaxID: {self.taxid}, Accession: {self.accession}, Description: {self.description}, Nucleotide ID: {self.nucleotide_id}, Assembly ID: {self.assembly_id}"


def retrieve_reference_sequence_id(accID: str) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    """
    Retrieve NCBI sequence ID for a given accession ID.
    """
    try:
        handle = Entrez.esearch(db="nucleotide", term=accID, retmax=1)
        record = Entrez.read(handle)
        handle.close()
        if not record['IdList']:
            raise ValueError(f"No sequence found for accession {accID}")
        sequence_id = record['IdList'][0]
        
        handle = Entrez.esummary(db="nucleotide", id=sequence_id)
        summary = Entrez.read(handle)
        handle.close()
        if summary is None:
            raise ValueError(f"No summary found for sequence ID {sequence_id}")
        accession = summary[0]['AccessionVersion']
        if accession != accID:
            print(f"Warning: Retrieved accession {accession} does not match requested {accID}")
        description = summary[0].get('Title', 'No description available')
        return accession, description, sequence_id


    except ValueError as e:
        print(e)
        return None, None, None
    except Exception as e:
        print(f"An error occurred: {e}")
        return None, None, None



def get_reference_sequence_url(taxid) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    """
    Retrieve the reference sequence URL for a given taxid from the nucleotide database.
    """
    try:
        # Search for nucleotide sequences for the given taxid
        handle = Entrez.esearch(db="nucleotide", term=f"txid{taxid}[Organism:exp] AND refseq", retmax=1)
        record = Entrez.read(handle)
        handle.close()
        if not record['IdList']:
            raise ValueError(f"No reference sequences found for taxid {taxid}")

        # Get the first sequence ID
        sequence_id = record['IdList'][0]

        # Fetch the summary for the sequence
        handle = Entrez.esummary(db="nucleotide", id=sequence_id)
        summary = Entrez.read(handle)
        handle.close()

        if summary is None:
            raise ValueError(f"No summary found for sequence ID {sequence_id}")

        # Extract the accession number
        docsum = summary[0]
        accession = docsum['AccessionVersion']
        description = docsum.get('Title', 'No description available')

        return accession, description, sequence_id

    except ValueError as e:
        print(e)
        return None, None, None
    except Exception as e:
        print(f"An error occurred: {e}")
        return None, None, None


def get_representative_assembly(taxid) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    try:
        # Search for assemblies for the given taxid
        handle = Entrez.esearch(db="assembly", term=f"txid{taxid}[Organism:exp]", retmax=5)
        record = Entrez.read(handle)
        handle.close()
        if not record['IdList']:
            raise ValueError(f"No representative genomes found for taxid {taxid}")
        # Get summary for the first assembly (could filter for 'representative genome' here)
        assembly_id = record['IdList'][0]

        handle = Entrez.esummary(db="assembly", id=assembly_id, report="full")
        summary = Entrez.read(handle)
        handle.close()
        docsum = summary['DocumentSummarySet']['DocumentSummary'][0]
        accession = docsum['AssemblyAccession']
        description = docsum.get('SpeciesName', 'No description available')

        return accession, description, assembly_id
    except Exception as e:
        print(f"An error occurred while fetching assembly for taxid {taxid}: {e}")
        return None, None, None
    except ValueError as e:
        print(e)
        return None, None, None


def retrieve_reference_sequence(nucleotide_id, output_path, gzipped=True) -> bool:
    """
    Download the reference sequence file given a nucleotide ID.
    """
    try:
        handle = Entrez.efetch(db="nucleotide", id=nucleotide_id, rettype="fasta", retmode="text")
        fasta_data = handle.read()
        handle.close()
        if gzipped:
            import gzip
            with gzip.open(output_path, 'wt') as f:
                f.write(fasta_data)
        else:
            with open(output_path, 'w') as f:
                f.write(fasta_data)
        return True
    except Exception as e:
        print(f"An error occurred while downloading sequence: {e}")
        return False


def retrieve_assembly_sequence(assembly_id, output_path) -> bool:
    """
    Download the assembly sequence file given an assembly ID.
    """
    try:
        handle = Entrez.esummary(db="assembly", id=assembly_id, report="full")
        summary = Entrez.read(handle)
        handle.close()
        docsum = summary['DocumentSummarySet']['DocumentSummary'][0]
        ftp_path = docsum['FtpPath_RefSeq'] or docsum['FtpPath_GenBank']
        if not ftp_path:
            print(f"No FTP path found for assembly ID {assembly_id}")
            return False
        asm_name = ftp_path.split('/')[-1]
        fasta_url = f"{ftp_path}/{asm_name}_genomic.fna.gz"
        # Download the file
        result = subprocess.run(['wget', '-O', output_path, fasta_url], capture_output=True, check= False)
        if result.returncode != 0:
            raise RuntimeError(f"Failed to download file: {result.stderr.decode()}")
        return True

    except Exception as e:
        print(f"An error occurred while downloading assembly: {e}")
        return False



class NCBITools:
    def __init__(self):
        self.logger = logging.getLogger('NCBITools')
        self.logger.setLevel(logging.DEBUG)
        ch = logging.StreamHandler()
        ch.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        ch.setFormatter(formatter)
        self.logger.addHandler(ch)
        self.logger.propagate = False

    def query_sequence_databases(self, passport:Passport) -> ReferenceData:
        """
        use both strategies above
        """

        if passport.accession is not None:

            accession, description, nucleotide_id = retrieve_reference_sequence_id(passport.accession)
            if nucleotide_id is not None:
                
                return ReferenceData(
                    taxid=passport.taxid,
                    accession=passport.accession,
                    description=description,
                    nucleotide_id=nucleotide_id
                )
            else:
                self.logger.warning(f"No nucleotide ID found for accession {passport.accession}, falling back to taxid search.")
            

        accession, description, nucleotide_id = get_reference_sequence_url(passport.taxid)

        if accession is not None and nucleotide_id is not None:
            return ReferenceData(
                taxid=passport.taxid,
                accession=accession,
                description=description,
                nucleotide_id=nucleotide_id
            )
        
        accession, description, assembly_id = get_representative_assembly(passport.taxid)
        
        return ReferenceData(
            taxid=passport.taxid,
            accession=accession,
            description=description,
            assembly_id=assembly_id
        )


    def retrieve_sequence_databases(self, reference_data:ReferenceData, output_path:str, gzipped=True) -> bool:
        """
        use both strategies above
        """

        if reference_data.nucleotide_id is not None:
            success = retrieve_reference_sequence(reference_data.nucleotide_id, output_path, gzipped)
            if not success:
                self.logger.warning(f"Failed to retrieve reference sequence for taxid {reference_data.taxid}")
            return success
        
        if reference_data.assembly_id is not None:
            success = retrieve_assembly_sequence(reference_data.assembly_id, output_path)
            if not success:
                self.logger.warning(f"Failed to retrieve assembly sequence for taxid {reference_data.taxid}")
            return success

        self.logger.warning(f"No sequence data found for taxid {reference_data.taxid}")
        return False
