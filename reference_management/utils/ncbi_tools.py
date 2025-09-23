# You need to install Biopython: pip install biopython
import logging
from Bio import Entrez
import os
import subprocess
from dataclasses import dataclass
from typing import Tuple, Optional
# read email from local .env file
import numpy as np
import dotenv
dotenv.load_dotenv()
import pandas as pd

Entrez.email = os.getenv("NCBI_EMAIL", None)
if Entrez.email is None:
    raise ValueError("NCBI_EMAIL environment variable not set. Please set it to your email address.")

NCBI_TAXONOMY_LEVELS = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

def compare_lineages(lineage1: Optional[str], lineage2: Optional[str]) -> tuple[float, Optional[str]]:
    """
    Compare two lineages. Split by ';' and compare each level.
    add points for each match, normalize by the length of the longer lineage.
    """
    
    if lineage1 is None or lineage2 is None:
        return 0.0, None

    if pd.isna(lineage1) or pd.isna(lineage2):
        return 0.0, None

    levels1 = lineage1.split('; ')
    levels2 = lineage2.split('; ')
    min_length = min(len(levels1), len(levels2))
    max_length = max(len(levels1), len(levels2))
    score = 0
    level= None
    for i in range(min_length):
        if i < min_length and levels1[i] == levels2[i]:
            score += 1

            if i < len(NCBI_TAXONOMY_LEVELS):
                level = NCBI_TAXONOMY_LEVELS[i]
            else:
                level = f"level_{i+1}"
        else:
            break
    if max_length == 0:
        return 0.0, None
    return score / (len(levels1)), level

@dataclass
class Passport:
    taxid: Optional[str]
    accession: Optional[str] = None
    lineage: Optional[str] = None
    description: Optional[str] = None

    def __str__(self):
        return f"TaxID: {self.taxid}, Accession: {self.accession}"

    @property
    def prefix(self):
        if self.accession is not None:
            return f"{self.taxid}_{self.accession}"
        else:
            return f"{self.taxid}"

    def compare_lineage(self, other_lineage: Optional[str]) -> tuple[float, Optional[str]]:
        return compare_lineages(self.lineage, other_lineage)

    
@dataclass
class LocalAssembly(Passport):

    file_path: Optional[str] = None

    def __str__(self):
        return f"TaxID: {self.taxid}, Accession: {self.accession}, File Path: {self.file_path}"

    
@dataclass
class ReferenceData(Passport):

    nucleotide_id: Optional[str] = None
    assembly_id: Optional[str] = None

    def __str__(self):
        return f"TaxID: {self.taxid}, Accession: {self.accession}, Description: {self.description}, Nucleotide ID: {self.nucleotide_id}, Assembly ID: {self.assembly_id}"


def retrieve_passport_taxonomy(taxid: str) -> Optional[str]:
    """
    Retrieve taxonomy information for a given taxid.
    """
    try:
        handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        if not records:
            raise ValueError(f"No taxonomy found for taxid {taxid}")
        lineage = records[0]['Lineage']
        # standardize lineage to taxonomy levels
        return lineage
    except Exception as e:
        print(f"An error occurred while fetching taxonomy for taxid {taxid}: {e}")
        return None


def retrieve_reference_sequence_id(accID: str, include_term = None, exclude_term=None) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    """
    Retrieve NCBI sequence ID for a given accession ID.
    """
    try:
        term = accID
        if exclude_term is not None:
            term += f" NOT {exclude_term}"

        if include_term is not None:
            term += f" AND {include_term}"
        
        print("Searching NCBI with term:", term)

        handle = Entrez.esearch(db="nucleotide", term=term, retmax=20)
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



def get_reference_sequence_url(taxid, include_term=None, exclude_term=None) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    """
    Retrieve the reference sequence URL for a given taxid from the nucleotide database.
    """
    try:
        term = f"txid{taxid}[Organism:exp] AND refseq"
        if include_term is not None:
            term += f" AND {include_term}"
        if exclude_term is not None:
            term += f" NOT {exclude_term}"

        print("Searching NCBI Nucleotide with term:", term)
        # Search for nucleotide sequences for the given taxid
        handle = Entrez.esearch(db="nucleotide", term=term, retmax=1)
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


def get_representative_assembly(taxid, include_term= None, exclude_term = None) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    try:
        # Search for assemblies for the given taxid
        term = f"txid{taxid}[Organism:exp]" 
        if exclude_term is not None:
            term += f" NOT {exclude_term}"
        if include_term is not None:
            term += f" AND {include_term}"
        
        print("Searching NCBI Assembly with term:", term)

        handle = Entrez.esearch(db="assembly", term=term, retmax=5)
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
    
    def retrieve_passport_taxonomy(self, passport:Passport) -> Optional[str]:
        """
        Retrieve taxonomy information for a given taxid.
        """
        try:
            handle = Entrez.efetch(db="taxonomy", id=passport.taxid, retmode="xml")
            records = Entrez.read(handle)
            handle.close()
            if not records:
                raise ValueError(f"No taxonomy found for taxid {passport.taxid}")
            lineage = records[0]['Lineage']
            return lineage
        except Exception as e:
            self.logger.error(f"An error occurred while fetching taxonomy for taxid {passport.taxid}: {e}")
            return None

    def query_sequence_databases(self, passport:Passport, include_term = None, exclude_term = None) -> ReferenceData:
        """
        use both strategies above
        """

        lineage = self.retrieve_passport_taxonomy(passport)
        passport.lineage = lineage

        if passport.accession is not None:

            accession, description, nucleotide_id = retrieve_reference_sequence_id(passport.accession, include_term=include_term, exclude_term=exclude_term)
            if nucleotide_id is not None:
                
                return ReferenceData(
                    taxid=passport.taxid,
                    accession=passport.accession,
                    lineage=passport.lineage,
                    description=description,
                    nucleotide_id=nucleotide_id
                )
            else:
                self.logger.warning(f"No nucleotide ID found for accession {passport.accession}, falling back to taxid search.")
            

        accession, description, nucleotide_id = get_reference_sequence_url(passport.taxid, include_term=include_term, exclude_term=exclude_term)

        if accession is not None and nucleotide_id is not None:
            passport.accession = accession
            return ReferenceData(
                taxid=passport.taxid,
                accession=accession,
                lineage=passport.lineage,
                description=description,
                nucleotide_id=nucleotide_id
            )
        print("No nucleotide sequence found, trying assembly...")
        accession, description, assembly_id = get_representative_assembly(passport.taxid)

        if accession is None and assembly_id is None:
            self.logger.error(f"No reference data found for taxid {passport.taxid}")
            return ReferenceData(
                taxid=passport.taxid,
                accession=None,
                lineage=passport.lineage,
                description=None,
                assembly_id=None
            )
        
        passport.accession = accession
        
        return ReferenceData(
            taxid=passport.taxid,
            accession=accession,
            lineage=passport.lineage,
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
