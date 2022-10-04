import json

import pandas as pd
import pprint
import sqlite3

import yaml

con = sqlite3.connect("../../../resources/felix_dump.db")

swissprot_entries_fp = "../../../resources/swiss_entries.json"

with open(swissprot_entries_fp, 'r') as infile:
    swissprot_entries = json.load(infile)

# pprint.pprint(swissprot_entries)

cur = con.cursor()

df = pd.read_sql_query("""
SELECT 
brbic.qacc as qacc, 
parts.alias as part_alias, 
parts_sequences_plus.seq_name as seq_name, 
parts_sequences_plus.type as seq_type,
parts_sequences_plus.date_added as seq_date_added,
parts_sequences_plus.sequence as nt_sequence,
brbic.sacc as sacc
FROM blast_results_best_in_cluster brbic
join main.parts_sequences_plus on brbic.qacc = parts_sequences_plus.id
join parts on parts_sequences_plus.part_id = parts.id
""", con)

# print(df)

blast_res_plus_lod = df.to_dict('records')

# pprint.pprint(blast_res_plus_lod[0])

dod = {}
for i in blast_res_plus_lod:
    if i['qacc'] not in dod:
        dod[i['qacc']] = {
            'id': f"sequence:{i['qacc']}",
            'associated_part': f"IF:{i['part_alias'].replace('IF', '')}",
            'date_added': i['seq_date_added'],
            'nt_sequence': i['nt_sequence'],
            'seq_name': i['seq_name'],
            'seq_type': i['seq_type'],
            'go_term_ids': set(),
            'go_term_labels': set(),
            'match_ec_numbers': set(),
            'match_gene_symbols_etc': set(),
            'match_names': set(),
            'scientific_names': set(),
            'taxon_ids': set(),
            'uniprot_accessions': set(),

        }
    dod[i['qacc']]['uniprot_accessions'].add(i['sacc'])
    print(i['sacc'])
    if i['sacc'] in swissprot_entries:
        # print(swissprot_entries[i['sacc']]['proteinDescription']['recommendedName']['fullName']['value'])
        dod[i['qacc']]['match_names'].add(
            swissprot_entries[i['sacc']]['proteinDescription']['recommendedName']['fullName']['value'])
        if 'organism' in swissprot_entries[i['sacc']]:
            dod[i['qacc']]['scientific_names'].add(swissprot_entries[i['sacc']]['organism']['scientificName'])
            dod[i['qacc']]['taxon_ids'].add(swissprot_entries[i['sacc']]['organism']['taxonId'])
        dbxrs = swissprot_entries[i['sacc']]['uniProtKBCrossReferences']
        for dbxr in dbxrs:
            if dbxr['database'] == 'GO':
                dod[i['qacc']]['go_term_ids'].add(dbxr['id'])
                for prop in dbxr['properties']:
                    if prop['key'] == 'GoTerm':
                        dod[i['qacc']]['go_term_labels'].add(prop['value'])
        if 'ecNumbers' in swissprot_entries[i['sacc']]['proteinDescription']['recommendedName']:
            ec_numbers_obj = swissprot_entries[i['sacc']]['proteinDescription']['recommendedName']['ecNumbers']
            for ec_number in ec_numbers_obj:
                dod[i['qacc']]['match_ec_numbers'].add(ec_number['value'])
        if 'genes' in swissprot_entries[i['sacc']]:
            # print(swissprot_entries[i['sacc']]['genes'])
            genes_obj = swissprot_entries[i['sacc']]['genes']
            for gene in genes_obj:
                if 'geneName' in gene:
                    dod[i['qacc']]['match_gene_symbols_etc'].add(gene['geneName']['value'])
                if 'synonyms' in gene:
                    for synonym in gene['synonyms']:
                        dod[i['qacc']]['match_gene_symbols_etc'].add(synonym['value'])
parts_seqs_lod = []
sets_to_lists = [
    'go_term_ids',
    'go_term_labels',
    'match_ec_numbers',
    'match_gene_symbols_etc',
    'match_names',
    'scientific_names',
    'taxon_ids',
    'uniprot_accessions',
]
for k, v in dod.items():
    for s2l_item in sets_to_lists:
        v[s2l_item] = list(v[s2l_item])
    parts_seqs_lod.append(v)

database = {"sequence_set": parts_seqs_lod}

with open("../../../resources/synbio_database.yaml", 'w') as outfile:
    yaml.dump(database, outfile)
