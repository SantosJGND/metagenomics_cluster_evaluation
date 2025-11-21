# You need to install Biopython: pip install biopython
import logging
from Bio import Entrez
import os
import subprocess
from dataclasses import dataclass
# read email from local .env file
from typing import List, Optional, Generator
import json
import os
import pandas as pd
import time
import dotenv
from typing import Tuple
dotenv.load_dotenv()

Entrez.email = os.getenv("NCBI_EMAIL", None)
if Entrez.email is None:
    raise ValueError("NCBI_EMAIL environment variable not set. Please set it to your email address.")

NCBI_TAXONOMY_LEVELS = ['acellular root', 'realm', 'domain', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

def estimate_level_pass(level, threshold_level = "species") -> bool:
    if threshold_level not in NCBI_TAXONOMY_LEVELS:
        raise ValueError(f"Invalid threshold level: {threshold_level}")
    if level not in NCBI_TAXONOMY_LEVELS:
        return True
    

    return NCBI_TAXONOMY_LEVELS.index(level) >= NCBI_TAXONOMY_LEVELS.index(threshold_level)
    
def get_lineage(taxid: str) -> Optional[str]:
    passport = Passport(taxid=taxid, accession=None)
    ncbi_tools = NCBITools()
    lineage = ncbi_tools.retrieve_passport_taxonomy(passport)
    return lineage


class NCBITaxonomistWrapper:
    """ Wrapper for ncbi-taxonomist command line tool. """
    def __init__(self, db: Optional[str] = None, temp_dir: Optional[str] = None):
        self.db_path = db
        self.temp_dir = temp_dir
        self.lineages: dict[int, dict[str, dict[str, str]]] = {}
        self.max_fetch = 50 # max taxids to fetch at once

    def retrieve_lineages_cmd_local(self, taxids: List[int]) -> str:
        """
        Construct the command to retrieve taxid lineage information using ncbi-taxonomist and a local database.
        """

        taxids_str = ",".join(str(t) for t in taxids)

        cmd = f"ncbi-taxonomist resolve -t {taxids_str} -db {self.db_path} "
        print("Retrieving taxids from local database:", cmd)

        return cmd

    def retrieve_lineages_cmd_import(self, taxids: List[int]) -> str:
        """
        Construct the command to import taxid lineage information from NCBI using ncbi-taxonomist"""
        taxids_str = ",".join(str(t) for t in taxids)

        cmd = f"ncbi-taxonomist resolve -t {taxids_str} --remote | ncbi-taxonomist import --database {self.db_path} "
        print("Importing taxids from NCBI:", cmd)

        return cmd
    

    def split_taxids(self, taxids: List[int], chunk_size=50) -> Generator[List[int], None, None]:
        """ Split a list of taxids into chunks of specified size. """
        for i in range(0, len(taxids), chunk_size):
            yield taxids[i:i + chunk_size]


    def resolve_lineages(self, taxids: List[int]) -> None:
        """
        Retrieve taxid lineage information using ncbi-taxonomist.
        - First, try to retrieve from local database.
        - If some taxids are missing, try to import them from NCBI.

        """
        taxids = [int(float(t)) for t in taxids]

        cmd = self.retrieve_lineages_cmd_local(taxids)

        stream = os.popen(cmd)
        output = stream.read()
        lineage_dict = self.parse_lineages_output(output)

        missing = set(taxids) - set(lineage_dict.keys())
        print("Missing taxids after local retrieval:", missing)
        if missing:
            print(f"Missing taxids: {missing}. Trying to import them from NCBI...")
            for taxid_chunk in self.split_taxids(list(missing), chunk_size=self.max_fetch):
                cmd = self.retrieve_lineages_cmd_import(list(taxid_chunk))
                stream = os.popen(cmd)
                output = stream.read()
                lineage_dict_update = self.parse_lineages_output(output)
                lineage_dict.update(lineage_dict_update)
                time.sleep(3)

        self.lineages.update(lineage_dict)

        still_missing = set(taxids) - set(self.lineages.keys())
        if still_missing:
            print(f"Still missing taxids: {still_missing}.")
        
        for taxid in still_missing:
            lineage = get_lineage(str(taxid))
            if lineage is not None:
                parser = NCBIlineageParser(lineage)
                self.lineages[taxid] = {
                    level: {'name': parser.get_level(level), 'taxid': None, 'level': i}
                    for i, level in enumerate(NCBI_TAXONOMY_LEVELS) if parser.get_level(level) is not None
                }
        
        print(f"Updated lineage dictionary with {len(self.lineages)} entries.")
        print(f"Still missing taxids after NCBI fetch: {set(taxids) - set(self.lineages.keys())}")

    def update_lineages(self, taxids: List[int]) -> None:
        """
        Update the lineage dictionary with new taxids.
        - First, check which taxids are missing from the current dictionary.
        - Try to retrieve from local database.
        - If some taxids are still missing, try to import them from NCBI."""
        missing = set(taxids) - set(self.lineages.keys())
        if missing:
            print(f"Missing taxids: {missing}. Trying to import them from NCBI...")
            for taxid_chunk in self.split_taxids(list(missing), chunk_size=self.max_fetch):
                cmd = self.retrieve_lineages_cmd_import(list(taxid_chunk))

                stream = os.popen(cmd)
                output = stream.read()
                lineage_dict_update = self.parse_lineages_output(output)
                self.lineages.update(lineage_dict_update)
                time.sleep(3)

        still_missing = set(taxids) - set(self.lineages.keys())
        if still_missing:
            print(f"Still missing taxids: {still_missing}.")
        
        for taxid in still_missing:
            lineage = get_lineage(str(taxid))
            if lineage is not None:
                parser = NCBIlineageParser(lineage)
                self.lineages[taxid] = {
                    level: {'name': parser.get_level(level), 'taxid': None, 'level': i}
                    for i, level in enumerate(NCBI_TAXONOMY_LEVELS) if parser.get_level(level) is not None
                }
        
        print(f"Updated lineage dictionary with {len(self.lineages)} entries.")
        print(f"Still missing taxids after NCBI fetch: {set(taxids) - set(self.lineages.keys())}")

    def parse_lineages_output(self, output: str) -> dict:
        """
        Parse the output from ncbi-taxonomist resolve command.
        - Return a dictionary mapping taxid to its lineage information.
            - Dict structure: {taxid: {rank: {'name': name, 'taxid': taxid, 'level': level}}}
            - Level is the depth in the lineage (0 for root, increasing downwards).
        """
        lineage_dict = {}
        for line in output.splitlines():
            try:
                data = json.loads(line)
                if 'query' in data and 'lineage' in data:
                    #lineage_dict[data['query']] = data['lineage']
                    len_lineage = len(data['lineage'])
                    lineage_dict[int(data['query'])] = {
                        x['rank']: {'name': x['name'], 'taxid': int(x['taxid']), 'level': len_lineage - i}
                        for i, x in enumerate(data['lineage'])
                    }
            except json.JSONDecodeError:
                continue
        return lineage_dict
    
    def get_name(self, taxid: int) -> Optional[str]:
        """
        Get the scientific name for the given taxid.
        Return None if the taxid is not found in the lineage dictionary.
        Return the deepest level name (species if available, else genus, etc.)
        """
        if taxid not in self.lineages:
            return None
        
        lineage = self.lineages.get(taxid, {})
        lineage = sorted(lineage.items(), key=lambda x: x[1]['level'], reverse=True)
        if lineage:
            return lineage[0][1]['name']
        return None

    def get_level(self, taxid: int, level: str) -> Optional[dict]:
        """
        Get the name of the given taxonomic level for the given taxid.
        Return None if the taxid or level is not found in the lineage dictionary.
        """
        if taxid not in self.lineages:
            return None
        return self.lineages.get(taxid, {}).get(level, {}).get('name', None)
    
    def get_lineage(self, taxid: int) -> list:
        """
        Get the full lineage as a list of (rank, name) tuples for the given taxid.
        Return an empty list if the taxid is not found in the lineage dictionary.
        """
        if taxid not in self.lineages:
            return []
        lineage = self.lineages.get(taxid, {})
        lineage = sorted(lineage.items(), key=lambda x: x[1]['level'])
        return [(rank, info['name']) for rank, info in lineage]
    
    def compare_lineages_relative(self, taxid1: int, taxid2: int) -> Tuple[float, Optional[str]]:
        """
        Compare two taxids based on their lineages.
        Return a similarity score between 0 and 1, and the deepest common taxonomic level name.
        Comparison is done relative to taxid1's lineage."""
        lineage1 = None

        if taxid1 not in self.lineages:
            return 0.0, None
        if taxid2 not in self.lineages:
            return 0.0, None

        min_levels = min(
            min([x['level'] for x in self.lineages[taxid1].values()]),
            min([x['level'] for x in self.lineages[taxid2].values()])
        )

        if min_levels == 0:
            lineage1 = self.get_lineage(taxid1)
            lineage2 = self.get_lineage(taxid2)
            tax, level = compare_lineages(
                "; ".join([x[1] for x in lineage1]),
                "; ".join([x[1] for x in lineage2])
            )
            return tax, level

        lineage1 = sorted(self.lineages[taxid1].items(), key=lambda x: x[1]['level'])

        for tax, level in lineage1:

            if self.get_level(taxid2, tax) != level['name']:
                return (level.get('level', 0) - 1) / len(lineage1), tax
        
        return 1.0, lineage1[-1][0] if lineage1 else None
    
    def level_is_atleast(level: str, threshold_level = "species") -> bool:
        if threshold_level not in NCBI_TAXONOMY_LEVELS:
            raise ValueError(f"Invalid threshold level: {threshold_level}")
        if level not in NCBI_TAXONOMY_LEVELS:
            return True

        return NCBI_TAXONOMY_LEVELS.index(level) >= NCBI_TAXONOMY_LEVELS.index(threshold_level)


class NCBIlineageParser:
    def __init__(self, lineage: Optional[str]):
        self.lineage = lineage
        self.levels = {}
        if lineage is not None:
            parts = [part.strip() for part in lineage.split(';')]
            for i, part in enumerate(parts):
                if i < len(NCBI_TAXONOMY_LEVELS):
                    self.levels[NCBI_TAXONOMY_LEVELS[i]] = part
                else:
                    self.levels[f"level_{i+1}"] = part

    def get_level(self, level: str) -> Optional[str]:
        return self.levels.get(level, None)
    
    def get_order(self) -> Optional[str]:
        return self.get_level('order')
    
    def get_family(self) -> Optional[str]:
        return self.get_level('family')

    def get_genus(self) -> Optional[str]:
        return self.get_level('genus')

    def get_species(self) -> Optional[str]:
        return self.get_level('species')

def compare_lineages(lineage1: Optional[str], lineage2: Optional[str]) -> Tuple[float, Optional[str]]:
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
    accession: Optional[str]
    lineage: Optional[str] = None
    description: Optional[str] = None


    def __init__(self, taxid: Optional[str], accession: Optional[str] = None, lineage: Optional[str] = None, description: Optional[str] = None):
        if "." in str(taxid):
            taxid = str(taxid).split(".")[0]
        self.taxid = taxid
        self.accession = accession
        self.lineage = lineage
        self.description = description

    def __str__(self):
        return f"TaxID: {self.taxid}, Accession: {self.accession}"

    @property
    def prefix(self):
        if self.accession:
            return f"{self.taxid}_{self.accession}"
        else:
            return f"{self.taxid}"

    def compare_lineage(self, other_lineage: Optional[str]) -> float:
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

    def query_sequence_databases(self, passport:Passport) -> ReferenceData:
        """
        use both strategies above
        """
        import numpy as np

        if passport.taxid is None or pd.isna(passport.taxid) or np.isnan(float(passport.taxid)):
            self.logger.error("No taxid provided in passport.")
            return ReferenceData(
                taxid=None,
                accession=None,
                lineage=None,
                description=None
            )
        

        lineage = self.retrieve_passport_taxonomy(passport)
        print(f"Lineage for taxid {passport.taxid}: {lineage}")
        passport.lineage = lineage

        if passport.accession is not None:

            accession, description, nucleotide_id = retrieve_reference_sequence_id(passport.accession)
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
            

        accession, description, nucleotide_id = get_reference_sequence_url(passport.taxid)

        if accession is not None and nucleotide_id is not None:
            return ReferenceData(
                taxid=passport.taxid,
                accession=accession,
                lineage=passport.lineage,
                description=description,
                nucleotide_id=nucleotide_id
            )
        
        accession, description, assembly_id = get_representative_assembly(passport.taxid)
        
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
