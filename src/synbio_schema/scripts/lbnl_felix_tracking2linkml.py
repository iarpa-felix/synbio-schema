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
from linkml_runtime import SchemaView
from linkml_runtime.dumpers import yaml_dumper
import validators
import math
import csv

logger = logging.getLogger(__name__)
click_log.basic_config(logger)


@click.command()
@click_log.simple_verbosity_option(logger)
@click.option("--sqlite_input_fp", type=click.Path(exists=True), default="resources/felix_dump.db")
@click.option("--swiss_entries_fp", type=click.Path(exists=True), default="resources/swiss_entries.json")
@click.option("--fpbase_entries_fp", type=click.Path(exists=True), default='resources/fpbase.json')
@click.option("--schema_fp", type=click.Path(exists=True), default='src/synbio_schema/schema/synbio_schema.yaml')
@click.option("--species_id_taxid_curations_fp", type=click.Path(exists=True),
              default='resources/curations/species_id_taxid_curations.tsv')
@click.option("--yaml_out", type=click.Path(), default="resources/synbio_database.yaml")
@click.option("--uniprot_cache_name", default='uniprot_entries_cache')
def cli(sqlite_input_fp: str, swiss_entries_fp: str, yaml_out: str, uniprot_cache_name: str, fpbase_entries_fp: str,
        schema_fp: str, species_id_taxid_curations_fp: str):
    # todo look for sequence:3609,IF:00429
    # todo disregarding asserted PI email. using email from the auth_user table.
    # could have taken a LinkML instantiation approach, using the model's data classes

    # todo further muddying waters by making a view
    synbio_view = SchemaView(schema_fp)

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
        # 'creator': 'creator_id',
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
    select
        "references",
        alias,
        bio_safety_level,
        creator_id,
        funding_source,
        genotype_phenotype,
        intellectual_property,
        keywords,
        name,
        notes,
        principal_investigator,
        group_concat(selection_marker, '|') as selection_marker,
        status,
        summary,
        external_url,
        group_concat(pa.accession || '_' || pa.source || '_' || pa.type, "|") as genome_accessions,
        host_species_id
    from parts p 
    left join selection_markers sm on p.id = sm.part_id
    left join external_urls eu on p.id = eu.part_id
    left join parts_accessions pa on p.id = pa.part_id
    where p."type" = 'STRAIN'
    group by
        "references",
        alias,
        bio_safety_level,
        creator_id,
        funding_source,
        genotype_phenotype,
        intellectual_property,
        keywords,
        name,
        notes,
        principal_investigator,
        status,
        summary,
        external_url,
        host_species_id
    """

    optional_strain_attribs = {
        'funding_source': 'funding_source',
        'genotype_phenotype': 'genotype_phenotype',
        'intellectual_property': 'intellectual_property',
        'keywords': 'keywords',
        'notes': 'notes',
        'references': 'references',
        'selection_markers': 'selection_marker',
        'status': 'status',
        'summary': 'summary',
        'external_urls': 'external_url',
        'genome_accessions': 'genome_accessions',
    }

    required_strain_attribs = {
        'bio_safety_level': 'bio_safety_level',
        'name': 'name',
    }

    blast_plus_required_attributes = {
        'date_added': 'seq_date_added',
        'nt_sequence': 'nt_sequence',
        'seq_name': 'seq_name',
        'seq_type': 'seq_type',
    }

    species_query = """
    select distinct
        s.id as id,
        species,
        strain,
        "comment",
        s."name" as name
    from
        species s
    join parts p on
        p.host_species_id = s.id;
    """

    org_required_attributes = {
        'name': 'species',
    }

    org_optional_attribs = {
        "comment": "comment",
        "species_name": "name",
        "strain_value": "strain",
    }

    ###

    session = requests_cache.CachedSession(uniprot_cache_name)

    con = sqlite3.connect(sqlite_input_fp)

    with open(swiss_entries_fp, 'r') as infile:
        swiss_entries = json.load(infile)

    best_blasts_frame = pd.read_sql_query(best_blasts_query, con)

    blast_res_plus_lod = best_blasts_frame.to_dict('records')

    for_sequence_set_dod = {}

    # with open(fpbase_entries_fp, 'r') as infile:
    #     fpbase_entries = json.load(infile)

    for blast_plus_record in blast_res_plus_lod:
        if blast_plus_record['qacc'] not in for_sequence_set_dod:
            for_sequence_set_dod[blast_plus_record['qacc']] = {
                'id': f"sequence:{blast_plus_record['qacc']}",
                'associated_part': f"IF:{blast_plus_record['part_alias'].replace('IF', '')}",
            }

            for bprak, bprav in blast_plus_required_attributes.items():
                for_sequence_set_dod[blast_plus_record['qacc']][bprak] = blast_plus_record[bprav]

            for s2l_item in annotation_sets_to_lists:
                for_sequence_set_dod[blast_plus_record['qacc']][s2l_item] = set()

        else:
            pass

        if blast_plus_record['blast_db'] == "swissprot":
            # logger.info(f"swissprot {for_sequence_set_dod[blast_plus_record['qacc']]['associated_part']}")
            for_sequence_set_dod[blast_plus_record['qacc']]['uniprot_accessions'].add(blast_plus_record['sacc'])

            if blast_plus_record['sacc'] in swiss_entries:

                for_sequence_set_dod[blast_plus_record['qacc']]['match_names'].add(
                    swiss_entries[blast_plus_record['sacc']]['proteinDescription']['recommendedName']['fullName'][
                        'value'])

                # logger.info(f"has entry {for_sequence_set_dod[blast_plus_record['qacc']]['match_names']}")

                if 'organism' in swiss_entries[blast_plus_record['sacc']]:
                    for_sequence_set_dod[blast_plus_record['qacc']]['scientific_names'].add(
                        swiss_entries[blast_plus_record['sacc']]['organism']['scientificName'])
                    for_sequence_set_dod[blast_plus_record['qacc']]['taxon_ids'].add(
                        swiss_entries[blast_plus_record['sacc']]['organism']['taxonId'])
                else:
                    logger.warning('no organism')

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
                logger.warning('no entry')
        else:
            logger.warning('not swissprot')
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
            "creator": f"person:{mod['creator_id']}",
        })

        for rmak, rmav in required_mod_attribs.items():
            for_modification_set[rmak] = mod[rmav]

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
            logger.error(
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

                if 'entryType' in swiss_entry:
                    for_modification_set['curated_uniprot_accession'] = mod['uniprot_id']
                else:
                    logger.warning(f"{mod['uniprot_id']} doesn't seem to be a Uniprot entry")
            else:
                logger.warning(
                    f"not propagating claimed 'uniprot_id' of {mod['uniprot_id']} for {for_modification_set['id']}")

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

    strain_induced_class = synbio_view.induced_class('Strain')
    for strain in strain_lod:
        if 'alias' in strain and strain['alias'] and pattern.match(strain['alias']) and 'host_species_id' in strain and \
                strain['host_species_id'] and not math.isnan(strain['host_species_id']):
            for_strain_set = {
                "id": f"MS:{strain['alias'].replace('MS_', '')}",
                "creator": f"person:{strain['creator_id']}",
                "host_organism": f"organism:{str(int(strain['host_species_id']))}",
            }
            for rsak, rsav in required_strain_attribs.items():
                for_strain_set[rsak] = strain[rsav]

            for osak, osav in optional_strain_attribs.items():
                if osav in strain and strain[osav] and strain[osav] not in missing_data_indicators:
                    slot_props = strain_induced_class.attributes[osak]
                    # todo but the delimiter might be different for plasmid attributes etc
                    #   maybe make this slot-specific
                    # if slot_props.multivalued:
                    if osak == 'selection_markers':
                        splitted = strain[osav].split("|")
                        for_strain_set[osak] = splitted
                    elif osak == 'external_urls':
                        splitted = re.split(" |;", strain[osav])
                        splitted = [i for i in splitted if i]
                        url_validated = [i.rstrip(".") for i in splitted if validators.url(i)]
                        url_failures = list(set(splitted) - set(url_validated))
                        logger.info(f"problematic external urls: {url_failures}")
                        logger.info(f"external urls with reasonable FORMATS: {url_validated}")
                        for_strain_set[osak] = url_validated
                    elif osak == 'genome_accessions':
                        splitted = strain[osav].split("|")

                        for one_split in splitted:
                            split_components = one_split.split("_")
                            if len(split_components) == 3:
                                if split_components[1] == 'NCBI Genome':
                                    if split_components[2] == 'assembly':

                                        split_accessions = split_components[0].split("-")
                                        if 'genome_accessions' in for_strain_set:
                                            for_strain_set['genome_accessions'].extend(split_accessions)
                                        else:
                                            for_strain_set['genome_accessions'] = split_accessions
                                    elif split_components[2] == 'sample':
                                        if 'biosample_accessions' in for_strain_set:
                                            for_strain_set['biosample_accessions'].append(split_components[0])
                                        else:
                                            for_strain_set['biosample_accessions'] = [split_components[0]]
                                    else:
                                        logger.warning(
                                            f"don't know the meaning of  {split_components[2]} in {split_components}")
                                else:
                                    logger.warning(
                                        f"don't know the meaning of  {split_components[1]} in{split_components}")
                            else:
                                logger.warning(
                                    f"don't know how to process the {len(one_split)} parts of {split_components} as an accession")
                    else:
                        for_strain_set[osak] = strain[osav]

            asserted_pi_lc = strain['principal_investigator'].lower()

            discovered_pi_id = people_frame.loc[people_frame['fullname_lc'] == asserted_pi_lc, 'id'].values

            if len(discovered_pi_id) != 1:
                logger.error(
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

    logger.warning(
        'Skipping these strain records, because they lack a host_species_id or because their alias value does not match the pattern \\MS_\\d+\\')
    logger.warning(pprint.pformat(bad_strains))

    # todo
    #   'host_species_id': nan,
    #   'principal_investigator_email': 'celniker@fruitfly.org',

    # DONE with strain_set

    org_frame = pd.read_sql_query(sql=species_query, con=con)

    species_id_taxid_curations = {}
    with open(species_id_taxid_curations_fp, newline='') as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter="\t")
        for row in reader:
            species_id_taxid_curations[row['id']] = row['ncbi']

    org_lod = org_frame.to_dict('records')
    org_set = []

    for org in org_lod:
        for_org_set = {
            "id": f"organism:{org['id']}",
        }

        for orak, orav in org_required_attributes.items():
            for_org_set[orak] = org[orav]

        for ooak, ooav in org_optional_attribs.items():
            if ooav in org and org[ooav] and org[ooav] not in missing_data_indicators:
                for_org_set[ooak] = org[ooav]

        if str(org['id']) in species_id_taxid_curations:
            for_org_set['strain_agnostic_taxid'] = species_id_taxid_curations[str(org['id'])]

        org_set.append(for_org_set)

    database = {
        "modification_set": modification_set,
        "person_set": person_set,
        "sequence_set": sequence_set,
        "strain_set": strain_set,
        "organism_set": org_set,
    }

    with open(yaml_out, 'w') as outfile:
        yaml.dump(database, outfile)


if __name__ == "__main__":
    cli()
