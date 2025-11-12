import pandas as pd
import networkx as nx
from abc import ABC
import pandas as pd
from Bio import Entrez
import os 
import dotenv
from typing import Optional
import urllib.error

dotenv.load_dotenv()

Entrez.email = os.getenv("NCBI_EMAIL", None)
if Entrez.email is None:
    raise ValueError("NCBI_EMAIL environment variable not set. Please set it to your email address.")

def taxid_to_description(taxid: int) -> Optional[str]:
    """
    Given a taxid, return the corresponding scientific name using NCBI Entrez.
    """

    if taxid >= 10**7:
        return None

    try:
        Entrez.email = os.getenv("NCBI_EMAIL", None)
        if Entrez.email is None:
            raise ValueError("NCBI_EMAIL environment variable not set. Please set it to your email address.")

        handle = Entrez.esummary(db="taxonomy", id=taxid)
        record = Entrez.read(handle)
        handle.close()
        return record[0]["ScientificName"]
    except RuntimeError as e:
        print(f"Error retrieving description for taxid {taxid}: {e}")
        return None
    except urllib.error.HTTPError as e:
        print(f"HTTP Error retrieving description for taxid {taxid}: {e}")
        return None

def protein_accession_to_taxid(accession: str) -> int:
    """
    Given a protein accession, return the corresponding taxid using NCBI Entrez.
    """
    Entrez.email = "your_email@example.com"
    handle = Entrez.esearch(db="protein", term=accession)
    record = Entrez.read(handle)
    handle.close()
    record_handle = Entrez.esummary(db="protein", id=record["IdList"])
    record = Entrez.parse(record_handle)
    taxid = None
    for rec in record:
        taxid = rec["TaxId"]
        taxid = int(taxid)
    
    return taxid

class ClassifierOutputProcesseor(ABC):

    def __init__(self, class_output_path):
        self.output_path  = class_output_path
        self.final_report: pd.DataFrame = pd.DataFrame(columns = ["description", "taxID"])
    
    @classmethod
    def from_file(cls, class_output_path):
        pass
    
    @classmethod
    def process(self):
        pass

    def prep_final_report(self):
        self.final_report = self.final_report.drop_duplicates(subset=["description"])
        self.final_report["description"] = self.final_report["description"].str.replace(" ", "_").str.lower()
        return self
    
    def save(self, output_path):

        self.final_report.to_csv(output_path, sep="\t", index=False)


class KrakenUniqueOutputProcessor(ClassifierOutputProcesseor):
    
    def __init__(self, class_output_path, min_uniq_reads: float = 10):
        super().__init__(class_output_path)
        self.min_uniq_reads = min_uniq_reads
    
    def from_file(self):    
        try:
            krakenunique_report = pd.read_csv(self.output_path, sep="\t")
            krakenunique_report.columns = ["status", "sequence_id", "taxID", "query_length", "map_length"]
            krakenunique_report = krakenunique_report[krakenunique_report["status"] == "C"]
            self.report = krakenunique_report
        except pd.errors.EmptyDataError:
            self.report = pd.DataFrame(columns=["status", "sequence_id", "taxID", "query_length", "map_length"])

        return self

    def process(self):
        summary = self.report.groupby("taxID").size().reset_index(name="Nreads")
        summary = summary[summary["Nreads"] > self.min_uniq_reads]
        summary['description'] = summary['taxID'].apply(taxid_to_description)
        summary = summary.dropna(subset=['description']).sort_values("Nreads", ascending=False).reset_index(drop=True)
        
        self.final_report = summary[["description", "taxID", "Nreads"]].rename(columns={"Nreads": "uniq_reads"})
        return self


class DiamondOutputProcessor(ClassifierOutputProcesseor):

    def __init__(self, class_output_path, min_uniq_reads: int =5):
        super().__init__(class_output_path)
        self.min_uniq_reads = min_uniq_reads
    
    def from_file(self):
        try:
            diamond_report = pd.read_csv(self.output_path, sep="\t", header=None)
            diamond_report.columns = ["query_id", "subject_id", "perc_identity", "alignment_length", "mismatches", "gap_opens", "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score"]

            self.report = diamond_report

        except pd.errors.EmptyDataError:
            self.report = pd.DataFrame(columns=["query_id", "subject_id", "perc_identity", "alignment_length", "mismatches", "gap_opens", "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score"])
        return self
    
    def process(self):
        summary = self.report.groupby("subject_id").size().reset_index(name="Nreads")
        summary = summary[summary["Nreads"] >= self.min_uniq_reads]
        summary['taxID'] = summary['subject_id'].apply(protein_accession_to_taxid)
        summary['description'] = summary['taxID'].apply(taxid_to_description)
        summary = summary.dropna(subset=['description']).sort_values("Nreads", ascending=False).reset_index(drop=True)
        self.final_report = summary[["description", "taxID", "Nreads"]].rename(columns={"Nreads": "uniq_reads"})

        return self


class CentrifugeOutputProcessor(ClassifierOutputProcesseor):

    def __init__(self, class_output_path, nuniq_threshold: int=5):
        super().__init__(class_output_path)
        self.nuniq_threshold = nuniq_threshold
    
    def from_file(self):
        try:
            centrifuge_report = pd.read_csv(self.output_path, sep="\t")
            self.report = centrifuge_report
        except pd.errors.EmptyDataError:
            self.report = pd.DataFrame(columns=["name", "taxID", "numUniqueReads"])
        return self
    
    def process(self):
        self.report.sort_values("numUniqueReads", ascending=False, inplace=True)
        self.report = self.report[self.report["numUniqueReads"] >= self.nuniq_threshold].reset_index(drop=True)
        self.final_report = self.report[["name", "taxID", "numUniqueReads"]].rename(columns={"name": "description", "numUniqueReads": "uniq_reads"})
        return self


def count_prefix_spaces(s):
    return len(s) - len(s.lstrip())

class KrakenOutputProcessor(ClassifierOutputProcesseor):

    def __init__(self, class_output_path, min_uniq_reads: float = 10):
        super().__init__(class_output_path)
        self.min_uniq_reads = min_uniq_reads
    
    def from_file(self):    
        try:
            kraken_report = pd.read_csv(self.output_path, sep="\t")
            kraken_report.columns = ["PercReads", "NumReadsRoot", "Nreads", "RankCode", "taxID", "name"]
            kraken_report["prefix_spaces"] = kraken_report["name"].apply(count_prefix_spaces)
        except pd.errors.EmptyDataError:
            kraken_report = pd.DataFrame(columns=["PercReads", "NumReadsRoot", "Nreads", "RankCode", "taxID", "name"])

        nodes, edges = self.kraken_report_to_tree(kraken_report)
        self.nodes_dict = {node[0]: node for node in nodes}

        self.edges = edges
        self.report = kraken_report
        self.nodes = nodes

        return self

    def get_node_info(self, node):

        return self.nodes_dict[node]


    def process(self):
        leaves_simple = self.get_simplified_leaves(self.nodes, self.edges)
        leaves_simple = {self.get_node_info(parent): [self.get_node_info(leaf) for leaf in leaves] for parent, leaves in leaves_simple.items()}
        leaves_summary = self.summarize_leaves(leaves_simple)
        leaves_summary = leaves_summary[leaves_summary["Nreads"] > self.min_uniq_reads]
        self.final_report = leaves_summary.rename(columns={"name": "description"})
        self.final_report = self.final_report[["description", "taxID", "Nreads"]].rename(columns={"Nreads": "uniq_reads"})
        return self

    @staticmethod
    def kraken_report_to_tree(report):
        taxid_dict = {}
        nodes_list = []
        edges_dict = {}
        edges = []

        for _, row in report.iterrows():
            name = row["name"].strip()
            tax_id = row["taxID"]
            prefix_spaces = count_prefix_spaces(row["name"])
            perc_reads = row["PercReads"]
            nreads = row["Nreads"]
            tax_rank = row["RankCode"]

            nodes_list.append((tax_id, name, prefix_spaces, perc_reads, nreads, tax_rank))
            taxid_dict[tax_id] = (name, prefix_spaces, perc_reads, nreads, tax_rank)

            # find parent node
            parent_node = None
            if prefix_spaces > 0:

                i = len(nodes_list) - 1
                while i > 0 and nodes_list[i][2] >= prefix_spaces:
                    i -= 1
                if i >= 0:
                    parent_node = nodes_list[i]
            
            if parent_node:
                parent_tax_id = parent_node[0]
                if parent_tax_id not in edges_dict:
                    edges_dict[parent_tax_id] = []
                edges_dict[parent_tax_id].append(tax_id)
                edges.append((parent_tax_id, tax_id))
            
        return nodes_list, edges

    @staticmethod
    def get_leaves(nodes, edges):
        G = nx.DiGraph()
        G.add_edges_from(edges)
        leaves = []
        for node in nodes:
            if G.out_degree(node[0]) == 0:
                leaves.append(node[0])
        return leaves

    @staticmethod
    def get_leaf_parents(nodes, edges):
        G = nx.DiGraph()
        G.add_edges_from(edges)
        leaf_parents = {}
        for node in nodes:
            if G.out_degree(node[0]) == 0:
                parent = list(G.predecessors(node[0]))
                if parent:
                    leaf_parents[node] = parent
        return leaf_parents


    def get_simplified_leaves(self, nodes, edges):
        """
        returns dict of {parent: [leaves]}
        """
        G = nx.DiGraph()
        G.add_edges_from(edges)
        ## starting from leaves:
        # - if leaf rank does not contain 'S', keep leaf.
        # - if leaf rank contains 'S', move up the until we find a parent without 'S' in its rank, keep the one before it.
        simplified_leaves = {}
        for node in nodes:
            if G.out_degree(node[0]) == 0:
                leaf = node[0]
                parent = list(G.predecessors(leaf))
                parent_detail =  self.get_node_info(parent[0]) if parent else None
                while parent and "S" in parent_detail[5]:
                    leaf = parent[0]
                    parent = list(G.predecessors(leaf))
                    parent_detail =  self.get_node_info(parent[0]) if parent else None
                if parent:
                    parent = parent[0]
                    if parent not in simplified_leaves:
                        simplified_leaves[parent] = []
                    simplified_leaves[parent].append(leaf)
        return simplified_leaves


    @staticmethod
    def summarize_leaves(leaves_simple) -> pd.DataFrame:
        """
        for each parent:
        if parent is species, return parent, else, return leaf with highest numUniqueReads
        """
        summary = []
        for parent, leaves in leaves_simple.items():
            if "S" in parent[5]:
                summary.append(parent)
            else:
                # find leaf with highest numUniqueReads
                best_leaf = max(leaves, key=lambda x: x[3])
                summary.append(best_leaf)
        return pd.DataFrame(summary, columns=["taxID", "description", "prefix_spaces", "perc_reads", "Nreads", "rank_code"])

