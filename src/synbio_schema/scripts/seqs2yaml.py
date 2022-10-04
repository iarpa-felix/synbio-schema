import json
import re

import pandas as pd
import pprint
import sqlite3

import yaml

# todo disregarding asserted PI email. using email from the auth_user table.

# could have taken an instantiation approach

pd.set_option('display.max_columns', None)

sqlite_input_fp = "resources/felix_dump.db"

swissprot_entries_fp = "resources/swiss_entries.json"

yaml_out = "resources/synbio_database.yaml"

sets_to_lists = [
    'go_term_ids',
    'go_term_labels',
    'match_ec_numbers',
    'match_gene_symbols_etc',
    'match_names',
    'other_accessions',
    'scientific_names',
    'taxon_ids',
    'uniprot_accessions',
]

con = sqlite3.connect(sqlite_input_fp)

with open(swissprot_entries_fp, 'r') as infile:
    swissprot_entries = json.load(infile)

cur = con.cursor()

df = pd.read_sql_query("""
SELECT 
brbic.qacc as qacc, 
brbic.sacc as sacc,
brbic.blast_db as blast_db,
parts.alias as part_alias, 
parts_sequences_plus.date_added as seq_date_added,
parts_sequences_plus.seq_name as seq_name, 
parts_sequences_plus.sequence as nt_sequence,
parts_sequences_plus.type as seq_type
FROM blast_results_best_in_cluster brbic
join main.parts_sequences_plus on brbic.qacc = parts_sequences_plus.id
join parts on parts_sequences_plus.part_id = parts.id
""", con)

blast_res_plus_lod = df.to_dict('records')

parts_seqs_dod = {}
for i in blast_res_plus_lod:
    if i['qacc'] not in parts_seqs_dod:
        parts_seqs_dod[i['qacc']] = {
            'id': f"sequence:{i['qacc']}",
            'associated_part': f"IF:{i['part_alias'].replace('IF', '')}",
            'date_added': i['seq_date_added'],
            'nt_sequence': i['nt_sequence'],
            'seq_name': i['seq_name'],
            'seq_type': i['seq_type'],
        }

    for s2l_item in sets_to_lists:
        parts_seqs_dod[i['qacc']][s2l_item] = set()

    print(i['sacc'])

    if i['sacc'] in swissprot_entries:
        parts_seqs_dod[i['qacc']]['uniprot_accessions'].add(i['sacc'])
        parts_seqs_dod[i['qacc']]['match_names'].add(
            swissprot_entries[i['sacc']]['proteinDescription']['recommendedName']['fullName']['value'])

        if 'organism' in swissprot_entries[i['sacc']]:
            parts_seqs_dod[i['qacc']]['scientific_names'].add(
                swissprot_entries[i['sacc']]['organism']['scientificName'])
            parts_seqs_dod[i['qacc']]['taxon_ids'].add(swissprot_entries[i['sacc']]['organism']['taxonId'])
        dbxrs = swissprot_entries[i['sacc']]['uniProtKBCrossReferences']
        for dbxr in dbxrs:
            if dbxr['database'] == 'GO':
                parts_seqs_dod[i['qacc']]['go_term_ids'].add(dbxr['id'])
                for prop in dbxr['properties']:
                    if prop['key'] == 'GoTerm':
                        parts_seqs_dod[i['qacc']]['go_term_labels'].add(prop['value'])

        if 'ecNumbers' in swissprot_entries[i['sacc']]['proteinDescription']['recommendedName']:
            ec_numbers_obj = swissprot_entries[i['sacc']]['proteinDescription']['recommendedName']['ecNumbers']
            for ec_number in ec_numbers_obj:
                parts_seqs_dod[i['qacc']]['match_ec_numbers'].add(ec_number['value'])

        if 'genes' in swissprot_entries[i['sacc']]:
            genes_obj = swissprot_entries[i['sacc']]['genes']
            for gene in genes_obj:
                if 'geneName' in gene:
                    parts_seqs_dod[i['qacc']]['match_gene_symbols_etc'].add(gene['geneName']['value'])
                if 'synonyms' in gene:
                    for synonym in gene['synonyms']:
                        parts_seqs_dod[i['qacc']]['match_gene_symbols_etc'].add(synonym['value'])
    else:
        if i['blast_db'] == 'uniprot':
            parts_seqs_dod[i['qacc']]['uniprot_accessions'].add(i['sacc'])
        if i['blast_db'] == 'fpbase':
            parts_seqs_dod[i['qacc']]['other_accessions'].add(i['sacc'])

parts_seqs_lod = []

for k, v in parts_seqs_dod.items():
    for s2l_item in sets_to_lists:
        v[s2l_item] = list(v[s2l_item])
    parts_seqs_lod.append(v)

modifications_query = """
select * from parts 
join modifications on parts.alias = modifications.element_id
left join auth_user on parts.creator_id = auth_user.id
where "type" = 'MODIFICATION'
"""

modifications_frame = pd.read_sql_query(sql=modifications_query, con=con)

modifications_lod = modifications_frame.to_dict('records')

modification_set = []

people_query = """
select * from auth_user
"""

people_frame = pd.read_sql_query(sql=people_query, con=con)

people_frame['fullname_lc'] = people_frame['first_name'].str.lower() + ' ' + people_frame['last_name'].str.lower()

for mod in modifications_lod:
    for_modification_set = ({
        "id": f"IF:{mod['alias'].replace('IF', '')}",
        "bio_safety_level": f"{mod['bio_safety_level']}",
        "creator": f"person:{mod['creator_id']}",
        "el_name_long": f"{mod['element']}",
        "status": f"{mod['status']}",
    })

    optional_mod_attribs = {
        'aa_change': 'aa_change',
        'category': 'category',
        'descriptor': 'descriptor',
        'notes': 'notes',
        'position': 'position',
        'subcategory_size': 'subcategory_size',
        'modification_type': 'modification_type',
        'el_name_short': 'element_name',
    }

    # todo add extra whitespace removal ?

    for omak, omav in optional_mod_attribs.items():
        if mod[omav] and mod[omav] not in ['None', 'NA']:
            for_modification_set[omak] = mod[omav]

    # todo should really try/except
    if mod['size'] and mod['size'] not in ['None', 'NA']:
        for_modification_set['size_bp'] = int(mod['size'].replace('bp', '').strip())

    asserted_pi_lc = mod['principal_investigator'].lower()

    discovered_pi_id = people_frame.loc[people_frame['fullname_lc'] == asserted_pi_lc, 'id'].values

    if len(discovered_pi_id) != 1:
        print('There must be one and only one principal investigator, with a name that matches the auth_users table')
        exit()
    else:
        for_modification_set['principal_investigator'] = f"person:{discovered_pi_id[0]}"

    modification_set.append(for_modification_set)

person_set = []
people_lod = people_frame.to_dict('records')
for person in people_lod:
    for_person_set = {
        "id": f"person:{person['id']}",
        # "date_joined": f"{person['date_joined']}",
        "email": f"{person['email']}",
        "first_name": f"{person['first_name']}",
        "is_staff": bool(person['is_staff']),
        "is_superuser": bool(person['is_superuser']),
        "last_name": f"{person['last_name']}",
        "username": f"{person['username']}",
    }
    person_set.append(for_person_set)

strain_query = """
select * from parts where "type" = 'STRAIN'
"""

strain_frame = pd.read_sql_query(sql=strain_query, con=con)

strain_lod = strain_frame.to_dict('records')
strain_set = []
pattern = re.compile(r"^MS_\d+$")

optional_strain_attribs = {
    'funding_source': 'funding_source',
    'genotype_phenotype': 'genotype_phenotype',
    'intellectual_property': 'intellectual_property',
    'keywords': 'keywords',
    'summary': 'summary',
    'references': 'references',
    'notes': 'notes',
    'status': 'status',
}

for strain in strain_lod:
    if 'alias' in strain and strain['alias'] and pattern.match(strain['alias']):
        for_strain_set = {
            "id": f"MS:{strain['alias'].replace('MS_', '')}",
            "name": f"{strain['name']}",
            "bio_safety_level": f"{strain['bio_safety_level']}",
            "creator": f"person:{strain['creator_id']}",
        }
        for osak, osav in optional_strain_attribs.items():
            if strain[osav] and strain[osav] not in [
                'N/A',
                'NA',
                'None',
            ]:
                for_strain_set[osak] = strain[osav]

        asserted_pi_lc = strain['principal_investigator'].lower()

        discovered_pi_id = people_frame.loc[people_frame['fullname_lc'] == asserted_pi_lc, 'id'].values

        if len(discovered_pi_id) != 1:
            print(
                'There must be one and only one principal investigator, with a name that matches the auth_users table')
            exit()
        else:
            for_strain_set['principal_investigator'] = f"person:{discovered_pi_id[0]}"

        strain_set.append(for_strain_set)
    else:
        print('Skipping strain whose alias does not match the pattern MS_\\d+')
        pprint.pprint(strain)

#  {'alias': 'MS_946',
#   'bio_safety_level': 'Level 1',
#   'creator_id': 11,
#   'funding_source': 'IARPA-FELIX',
#   'genotype_phenotype': None,
#   'host_species_id': nan,
#   'id': 3592,
#   'intellectual_property': None,
#   'keywords': None,
#   'name': 'MS_946',
#   'notes': None,
#   'principal_investigator': 'Susan Celniker',
#   'principal_investigator_email': 'celniker@fruitfly.org',
#   'references': None,
#   'status': 'In Progress',
#   'summary': '',
#   'type': 'STRAIN'}]

database = {
    "modification_set": modification_set,
    "person_set": person_set,
    "sequence_set": parts_seqs_lod,
    "strain_set": strain_set,
}

with open(yaml_out, 'w') as outfile:
    yaml.dump(database, outfile)
