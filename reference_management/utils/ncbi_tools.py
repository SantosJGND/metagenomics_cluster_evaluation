# You need to install Biopython: pip install biopython
from Bio import Entrez
import os
import subprocess
from dataclasses import dataclass

# read email from local .env file

import dotenv
dotenv.load_dotenv()


Entrez.email = os.getenv("NCBI_EMAIL", None)
if Entrez.email is None:
    raise ValueError("NCBI_EMAIL environment variable not set. Please set it to your email address.")

@dataclass
class ReferenceData:
    taxid: str = "None"
    accession: str = "None"
    description: str = "None"
    nucleotide_id: str = None
    assembly_id: str = None

    def __str__(self):
        return f"Accession: {self.accession}, Description: {self.description}, Nucleotide ID: {self.nucleotide_id}, Assembly ID: {self.assembly_id}"


def get_reference_sequence_url(taxid):
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
        print(f"Found sequence ID: {sequence_id} for taxid {taxid}")

        # Fetch the summary for the sequence
        handle = Entrez.esummary(db="nucleotide", id=sequence_id)
        summary = Entrez.read(handle)
        handle.close()

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


def retrieve_reference_sequence(nucleotide_id, output_path, gzipped=True):
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

def get_representative_assembly(taxid):
    try:
        # Search for assemblies for the given taxid
        handle = Entrez.esearch(db="assembly", term=f"txid{taxid}[Organism:exp]", retmax=5)
        record = Entrez.read(handle)
        handle.close()
        if not record['IdList']:
            raise ValueError(f"No representative genomes found for taxid {taxid}")
        # Get summary for the first assembly (could filter for 'representative genome' here)
        assembly_id = record['IdList'][0]
        print(f"Found assembly ID: {assembly_id} for taxid {taxid}")
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

def retrieve_assembly_sequence(assembly_id, output_path):
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

def query_sequence_databases(taxid:str) -> ReferenceData:
    """
    use both strategies above
    """

    accession, description, nucleotide_id = get_reference_sequence_url(taxid)
    print(f"Accession: {accession}")
    print(f"Description: {description}")
    if accession is not None and nucleotide_id is not None:
        return ReferenceData(
            taxid=taxid,
            accession=accession,
            description=description,
            nucleotide_id=nucleotide_id
        )
    accession, description, assembly_id = get_representative_assembly(taxid)


    return ReferenceData(
        taxid=taxid,
        accession=accession,
        description=description,
        assembly_id=assembly_id
    )

def retrieve_sequence_databases(reference_data:ReferenceData, output_path:str, gzipped=True):
    """
    use both strategies above
    """

    if reference_data.nucleotide_id is not None:
        success = retrieve_reference_sequence(reference_data.nucleotide_id, output_path, gzipped)
        if not success:
            print(f"Failed to retrieve reference sequence for taxid {reference_data.taxid}")
        return success, output_path
    
    if reference_data.assembly_id is not None:
        success = retrieve_assembly_sequence(reference_data.assembly_id, output_path)
        if not success:
            print(f"Failed to retrieve assembly sequence for taxid {reference_data.taxid}")
        return success, output_path
    
    print(f"No sequence data found for taxid {reference_data.taxid}")
    return False, output_path
