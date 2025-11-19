import pandas as pd

from utils.ncbi_tools import NCBITaxonomistWrapper, compare_lineages


from typing import Optional

def get_lineage(taxid: str) -> Optional[str]:
    passport = Passport(taxid=taxid, accession=None)
    ncbi_tools = NCBITools()
    lineage = ncbi_tools.retrieve_passport_taxonomy(passport)
    return lineage

def find_assembly_mapping(row, stats_matrix):
    accession = row['assembly_accession']
    if accession is None or pd.isna(accession):
        row['clade'] = 'unmapped'
        row['nuniq'] = 0
        row['freq'] = 0
        row['Min_Pairwise_Dist'] = 0
        row['nleaves'] = 0
        return row

    match = stats_matrix[stats_matrix['leaves'].str.contains(accession, na=False) | (stats_matrix['clade'].str.contains(accession, na=False))]

    if match.empty:
        row['clade'] = None
        row['nuniq'] = 0
        row['freq'] = 0
        row['Min_Pairwise_Dist'] = 0
        row['nleaves'] = 0
    else:
        row['clade'] = match['clade'].values[0]
        row['nuniq'] = match['nuniq'].values[0]
        row['freq'] = match['freq'].values[0]
        row['Min_Pairwise_Dist'] = match['Min_Pairwise_Dist'].values[0]
        row['nleaves'] = match['nleaves'].values[0]

    return row

def find_best_match(taxid1, taxid_list, ncbi_wrapper: NCBITaxonomistWrapper) -> tuple[Optional[int], Optional[str], float]:
    best_taxid = None
    best_level = None
    best_score = 0.0
    for taxid2 in taxid_list:
        score, level = ncbi_wrapper.compare_lineages_relative(taxid2, taxid1)
        if level is not None and score > best_score:
            best_level = level
            best_taxid = taxid2
            best_score = score
    return best_taxid, best_level, best_score

def update_df_best_match(row, taxid_list, ncbi_wrapper: NCBITaxonomistWrapper, min_score = 0.7):
    taxid = row['taxid']
    if pd.isna(taxid):
        row['best_match_taxid'] = None
        row['best_match_level'] = -1
        row['best_match_score'] = 0.0
        row['name'] = None

    best_taxid, best_level, best_score = find_best_match(taxid, taxid_list, ncbi_wrapper)

    if best_score < min_score:
        best_taxid = None
        best_level = -1
        best_score = 0.0

    row['best_match_taxid'] = best_taxid
    row['best_match_level'] = best_level
    row['best_match_score'] = best_score
    
    if best_taxid is not None:
        row['name'] = ncbi_wrapper.get_name(best_taxid)
    else:
        row['name'] = None

    return row


def match_leaf(accid, list):
    for item in list:
        if accid in item:
            return item
    return None 