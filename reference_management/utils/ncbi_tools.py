# You need to install Biopython: pip install biopython
from Bio import Entrez

Entrez.email = "joao.dourado@insa.min-saude.pt"  # Always set your email

def get_representative_assembly(taxid):
    # Search for assemblies for the given taxid
    handle = Entrez.esearch(db="assembly", term=f"txid{taxid}[Organism:exp]", retmax=5)
    record = Entrez.read(handle)
    handle.close()
    if not record['IdList']:
        print(f"No assemblies found for taxid {taxid}")
        return None, None
    # Get summary for the first assembly (could filter for 'representative genome' here)
    assembly_id = record['IdList'][0]
    print(f"Found assembly ID: {assembly_id} for taxid {taxid}")
    handle = Entrez.esummary(db="assembly", id=assembly_id, report="full")
    summary = Entrez.read(handle)
    handle.close()
    docsum = summary['DocumentSummarySet']['DocumentSummary'][0]
    ftp_path = docsum['FtpPath_RefSeq'] or docsum['FtpPath_GenBank']
    accession = docsum['AssemblyAccession']
    if not ftp_path:
        print(f"No FTP path found for taxid {taxid}")
        return accession, None
    # Find the genomic fasta file
    asm_name = ftp_path.split('/')[-1]
    fasta_url = f"{ftp_path}/{asm_name}_genomic.fna.gz"
    return accession, fasta_url
