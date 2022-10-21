import json
import logging
import pprint
import re
import sqlite3

import click
import click_log
import pandas as pd
import requests_cache
import yaml

logger = logging.getLogger(__name__)
click_log.basic_config(logger)


@click.command()
@click_log.simple_verbosity_option(logger)
@click.option("--sqlite_input_fp", type=click.Path(exists=True), default="resources/felix_dump.db")
@click.option("--swiss_entries_fp", type=click.Path(exists=True), default="resources/swiss_entries.json")
@click.option("--fpbase_entries_fp", type=click.Path(exists=True), default='resources/fpbase.json')
@click.option("--yaml_out", type=click.Path(), default="resources/synbio_database.yaml")
@click.option("--uniprot_cache_name", default='uniprot_entries_cache')
def cli(sqlite_input_fp: str, swiss_entries_fp: str, yaml_out: str, uniprot_cache_name: str, fpbase_entries_fp: str):
    # todo look for sequence:3609,IF:00429
    # todo disregarding asserted PI email. using email from the auth_user table.
    # could have taken a LinkML instantiation approach, using the model's data classes

    # sqlite_input_fp = "resources/felix_dump.db"
    # swiss_entries_fp = "resources/swiss_entries.json"
    # uniprot_cache_name = 'uniprot_entries_cache'
    # yaml_out = "resources/synbio_database.yaml"

    pd.set_option('display.max_columns', None)

    missing_data_indicators = [
        'N/A',
        'NA',
        'None',
        'none',
    ]

    annotation_sets_to_lists = [
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

    subparts_query = """select
        ppl.alias as parent_alias,
        spl.alias as subpart_alias
    from
        sub_parts sp
    join parts ppl on
        sp.parent_part_id = ppl.id
    join parts spl on
        sp.sub_part_id = spl.id
    where
        ppl."type" = 'STRAIN'
        and spl."type" = 'MODIFICATION';"""

    best_blasts_query = """
    SELECT 
    brbic.qacc as qacc, 
    brbic.sacc as sacc,
    brbic.qcovs as qcovs,
    brbic.pident as pident,
    brbic.blast_db as blast_db,
    parts.alias as part_alias, 
    parts_sequences_plus.date_added as seq_date_added,
    parts_sequences_plus.seq_name as seq_name, 
    parts_sequences_plus.sequence as nt_sequence,
    parts_sequences_plus.type as seq_type
    FROM blast_results_best_in_cluster brbic
    join main.parts_sequences_plus on brbic.qacc = parts_sequences_plus.id
    join parts on parts_sequences_plus.part_id = parts.id
    """

    modifications_query = """
    select
        *
    from
        parts
    join modifications on
        parts.alias = modifications.element_id
    left join auth_user on
        parts.creator_id = auth_user.id
    left join modifications_genes mg on 
        modifications.modifications_genes_id = mg.id
    where
        "type" = 'MODIFICATION';
    """

    required_mod_attribs = {
        'bio_safety_level': 'bio_safety_level',
        'creator': 'creator_id',
        'el_name_long': 'element',
        'status': 'status',
    }

    optional_mod_attribs = {
        'aa_change': 'aa_change',
        'category': 'category',
        'curated_enzyme_name': 'enzyme_name',
        'curated_protein_name': 'protein_name',
        'descriptor': 'descriptor',
        'el_name_short': 'element_name',
        'modification_type': 'modification_type',
        'notes': 'notes',
        'position': 'position',
        'subcategory_size': 'subcategory_size',
    }

    people_query = """
    select * from auth_user
    """

    required_person_attribs = {
        'email': 'email',
        'first_name': 'first_name',
        'last_name': 'last_name',
        'username': 'username',
    }

    strain_query = """
    select * from parts where "type" = 'STRAIN'
    """

    optional_strain_attribs = {
        'funding_source': 'funding_source',
        'genotype_phenotype': 'genotype_phenotype',
        'intellectual_property': 'intellectual_property',
        'keywords': 'keywords',
        'notes': 'notes',
        'references': 'references',
        'status': 'status',
        'summary': 'summary',
    }

    required_strain_attribs = {
        'bio_safety_level': 'bio_safety_level',
        'creator': 'creator_id',
        'name': 'name',
    }

    blast_plus_required_attributes = {
        'date_added': 'seq_date_added',
        'nt_sequence': 'nt_sequence',
        'seq_name': 'seq_name',
        'seq_type': 'seq_type',
    }

    ###

    session = requests_cache.CachedSession(uniprot_cache_name)

    con = sqlite3.connect(sqlite_input_fp)

    with open(swiss_entries_fp, 'r') as infile:
        swiss_entries = json.load(infile)

    best_blasts_frame = pd.read_sql_query(best_blasts_query, con)

    # best_blasts_frame.to_csv("best_blasts_frame.tsv", sep="\t")

    blast_res_plus_lod = best_blasts_frame.to_dict('records')

    for_sequence_set_dod = {}

    with open(fpbase_entries_fp, 'r') as infile:
        fpbase_entries = json.load(infile)

    for blast_plus_record in blast_res_plus_lod:
        if blast_plus_record['qacc'] not in for_sequence_set_dod:
            # print(f"seeding blast_plus_record {blast_plus_record['qacc']}")
            for_sequence_set_dod[blast_plus_record['qacc']] = {
                'id': f"sequence:{blast_plus_record['qacc']}",
                'associated_part': f"IF:{blast_plus_record['part_alias'].replace('IF', '')}",
            }

            for bprak, bprav in blast_plus_required_attributes.items():
                for_sequence_set_dod[blast_plus_record['qacc']][bprak] = blast_plus_record[bprav]

            for s2l_item in annotation_sets_to_lists:
                for_sequence_set_dod[blast_plus_record['qacc']][s2l_item] = set()

            # pprint.pprint(for_sequence_set_dod)
        else:
            # print(f"annotating {blast_plus_record['qacc']}")
            pass

        if blast_plus_record['blast_db'] == "swissprot":
            print(f"swissprot {for_sequence_set_dod[blast_plus_record['qacc']]['associated_part']}")
            for_sequence_set_dod[blast_plus_record['qacc']]['uniprot_accessions'].add(blast_plus_record['sacc'])

            if blast_plus_record['sacc'] in swiss_entries:

                for_sequence_set_dod[blast_plus_record['qacc']]['match_names'].add(
                    swiss_entries[blast_plus_record['sacc']]['proteinDescription']['recommendedName']['fullName'][
                        'value'])

                print(f"has entry {for_sequence_set_dod[blast_plus_record['qacc']]['match_names']}")

                if 'organism' in swiss_entries[blast_plus_record['sacc']]:
                    print('has organism')
                    for_sequence_set_dod[blast_plus_record['qacc']]['scientific_names'].add(
                        swiss_entries[blast_plus_record['sacc']]['organism']['scientificName'])
                    for_sequence_set_dod[blast_plus_record['qacc']]['taxon_ids'].add(
                        swiss_entries[blast_plus_record['sacc']]['organism']['taxonId'])
                else:
                    print('no organism')

                dbxrs = swiss_entries[blast_plus_record['sacc']]['uniProtKBCrossReferences']

                for dbxr in dbxrs:
                    if dbxr['database'] == 'GO':
                        for_sequence_set_dod[blast_plus_record['qacc']]['go_term_ids'].add(dbxr['id'])
                        for prop in dbxr['properties']:
                            if prop['key'] == 'GoTerm':
                                for_sequence_set_dod[blast_plus_record['qacc']]['go_term_labels'].add(prop['value'])

                if 'ecNumbers' in swiss_entries[blast_plus_record['sacc']]['proteinDescription']['recommendedName']:
                    ec_numbers_obj = swiss_entries[blast_plus_record['sacc']]['proteinDescription']['recommendedName'][
                        'ecNumbers']
                    for ec_number in ec_numbers_obj:
                        for_sequence_set_dod[blast_plus_record['qacc']]['match_ec_numbers'].add(ec_number['value'])

                if 'genes' in swiss_entries[blast_plus_record['sacc']]:
                    genes_obj = swiss_entries[blast_plus_record['sacc']]['genes']
                    for gene in genes_obj:
                        if 'geneName' in gene:
                            for_sequence_set_dod[blast_plus_record['qacc']]['match_gene_symbols_etc'].add(
                                gene['geneName']['value'])
                        if 'synonyms' in gene:
                            for synonym in gene['synonyms']:
                                for_sequence_set_dod[blast_plus_record['qacc']]['match_gene_symbols_etc'].add(
                                    synonym['value'])
            else:
                print('no entry')
        else:
            print('not swissprot')
            for_sequence_set_dod[blast_plus_record['qacc']]['other_accessions'].add(blast_plus_record['sacc'])

    sequence_set = []

    for k, v in for_sequence_set_dod.items():
        for s2l_item in annotation_sets_to_lists:
            v[s2l_item] = list(v[s2l_item])
        sequence_set.append(v)

    # DONE with sequence_set

    # BUILD modifications_set

    modifications_frame = pd.read_sql_query(sql=modifications_query, con=con)

    modifications_lod = modifications_frame.to_dict('records')

    # todo
    # /Users/MAM/Documents/gitrepos/semantic-synbio/synbio-schema/src/synbio_schema/scripts/celniker2yaml.py:149: UserWarning: DataFrame columns are not unique, some columns will be omitted.
    #   modifications_lod = modifications_frame.to_dict('records')

    modification_set = []

    people_frame = pd.read_sql_query(sql=people_query, con=con)

    people_frame['fullname_lc'] = people_frame['first_name'].str.lower() + ' ' + people_frame['last_name'].str.lower()

    subparts_frame = pd.read_sql_query(subparts_query, con)
    subparts_frame['parent_alias'] = subparts_frame['parent_alias'].str.replace("MS_", "MS:")
    subparts_frame['subpart_alias'] = subparts_frame['subpart_alias'].str.replace("IF", "IF:")
    subparts_lod = subparts_frame.to_dict('records')
    has_parts_dos = {}
    part_ofs_dos = {}
    for subparts_row in subparts_lod:
        if subparts_row['parent_alias'] not in has_parts_dos:
            has_parts_dos[subparts_row['parent_alias']] = set()
        has_parts_dos[subparts_row['parent_alias']].add(subparts_row['subpart_alias'])
        if subparts_row['subpart_alias'] not in part_ofs_dos:
            part_ofs_dos[subparts_row['subpart_alias']] = set()
        part_ofs_dos[subparts_row['subpart_alias']].add(subparts_row['parent_alias'])

    for mod in modifications_lod:
        for_modification_set = ({
            "id": f"IF:{mod['alias'].replace('IF', '')}",
        })

        for rmak, rmav in required_mod_attribs.items():
            for_modification_set[rmak] = mod[rmav] \
                # .strip()

        # todo add extra whitespace removal ?

        for omak, omav in optional_mod_attribs.items():
            if mod[omav] and mod[omav] not in missing_data_indicators:
                for_modification_set[omak] = mod[omav] \
                    # .strip()

        # todo should really try/except
        if mod['size'] and mod['size'] not in missing_data_indicators:
            for_modification_set['size_bp'] = int(mod['size'].replace('bp', '').strip())

        asserted_pi_lc = mod['principal_investigator'].lower()

        discovered_pi_id = people_frame.loc[people_frame['fullname_lc'] == asserted_pi_lc, 'id'].values

        if len(discovered_pi_id) != 1:
            print(
                'There must be one and only one principal investigator, with a name that matches the auth_users table')
            exit()
        else:
            for_modification_set['principal_investigator'] = f"person:{discovered_pi_id[0]}"

        gene_symbol_set = set()
        if 'gene_name' in mod and mod['gene_name']:
            gene_symbol_set.add(mod['gene_name'])
        if 'standardized_gene_name' in mod and mod['standardized_gene_name']:
            gene_symbol_set.add(mod['standardized_gene_name'])

        if gene_symbol_set:
            for_modification_set['curated_gene_symbols'] = list(gene_symbol_set)

        if for_modification_set['id'] in part_ofs_dos and part_ofs_dos[for_modification_set['id']]:
            for_modification_set['part_ofs'] = list(part_ofs_dos[for_modification_set['id']])

        if 'uniprot_id' in mod and mod['uniprot_id']:
            if mod['uniprot_id'].isalnum():
                swiss_entry = session.get(f"https://www.uniprot.org/uniprot/{mod['uniprot_id']}.json").json()
                #   "messages": [
                #     "The 'accession' value has invalid format. It should be a valid UniProtKB accession"
                #   ]
                if 'entryType' in swiss_entry:
                    for_modification_set['curated_uniprot_accession'] = mod['uniprot_id']
                else:
                    print(f"{mod['uniprot_id']} doesn't seem to be a Uniprot entry")
            else:
                print(f"not propagating claimed 'uniprot_id' of {mod['uniprot_id']} for {for_modification_set['id']}")

        modification_set.append(for_modification_set)

    # DONE with modification_set

    # BUILD person_set

    person_set = []
    people_lod = people_frame.to_dict('records')

    for person in people_lod:
        for_person_set = {
            "id": f"person:{person['id']}",
            "is_staff": bool(person['is_staff']),
            "is_superuser": bool(person['is_superuser']),
            # todo keep and fix date_joined
            # "date_joined": f"{person['date_joined']}",
        }

        # todo distinguish required from optional
        for rpak, rpav in required_person_attribs.items():
            if rpav in person and person[rpav] and person[rpav] not in missing_data_indicators:
                for_person_set[rpak] = person[rpav] \
                    # .strip()

        person_set.append(for_person_set)

    # DONE with person_set

    # BUILD strain_set

    strain_frame = pd.read_sql_query(sql=strain_query, con=con)

    strain_lod = strain_frame.to_dict('records')
    strain_set = []
    pattern = re.compile(r"^MS_\d+$")

    bad_strains = []

    for strain in strain_lod:
        if 'alias' in strain and strain['alias'] and pattern.match(strain['alias']):
            for_strain_set = {
                "id": f"MS:{strain['alias'].replace('MS_', '')}",
            }
            for rsak, rsav in required_strain_attribs.items():
                for_strain_set[rsak] = strain[rsav] \
                    # .strip()

            for osak, osav in optional_strain_attribs.items():
                if strain[osav] and strain[osav] not in missing_data_indicators:
                    for_strain_set[osak] = strain[osav] \
                        # .strip()

            asserted_pi_lc = strain['principal_investigator'].lower()

            discovered_pi_id = people_frame.loc[people_frame['fullname_lc'] == asserted_pi_lc, 'id'].values

            if len(discovered_pi_id) != 1:
                print(
                    'There must be one and only one principal investigator, with a name that matches the auth_users table')
                exit()
            else:
                for_strain_set['principal_investigator'] = f"person:{discovered_pi_id[0]}"

            if for_strain_set['id'] in has_parts_dos and has_parts_dos[for_strain_set['id']]:
                for_strain_set['has_parts'] = list(has_parts_dos[for_strain_set['id']])

            strain_set.append(for_strain_set)
        else:
            # todo dump these to a file instead of printing
            bad_strains.append(strain)

    print('Skipping these strain records, because their alias value do not match the pattern MS_\\d+')
    pprint.pprint(bad_strains)

    # todo
    #   'host_species_id': nan,
    #   'principal_investigator_email': 'celniker@fruitfly.org',

    # DONE with strain_set

    database = {
        "modification_set": modification_set,
        "person_set": person_set,
        "sequence_set": sequence_set,
        "strain_set": strain_set,
    }

    # print("about to write")

    with open(yaml_out, 'w') as outfile:
        yaml.dump(database, outfile)


if __name__ == "__main__":
    cli()
