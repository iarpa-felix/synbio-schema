import pprint
import json
import sqlite3

import pandas as pd
import requests_cache
import yaml

import logging

import click
import click_log

import pandas as pd

# from typing import List

pd.set_option('display.max_columns', None)

logger = logging.getLogger(__name__)
click_log.basic_config(logger)


@click.command()
@click_log.simple_verbosity_option(logger)
@click.option("--sqlite_db_fp", type=click.Path(), default="resources/felix_dump.db")
@click.option("--swiss_entries_dump_fp", type=click.Path(), default="resources/swiss_entries.json")
@click.option("--blast_results_table", default="blast_results")
@click.option("--min_ident", default=50)
@click.option("--requests_cache_name", default="uniprot_entries_cache", help="basename of SQLite file")
def cli(sqlite_db_fp: str, swiss_entries_dump_fp: str, blast_results_table: str, requests_cache_name: str,
        min_ident: float):
    session = requests_cache.CachedSession(requests_cache_name)

    sqlite_db_conn = sqlite3.connect(sqlite_db_fp)
    swiss_accessions_q = f"select distinct sacc from {blast_results_table} where blast_db = 'swissprot' and pident >= {str(min_ident)} order by sacc"
    swiss_accessions_res = pd.read_sql_query(swiss_accessions_q, sqlite_db_conn)
    swiss_accessions = swiss_accessions_res['sacc'].to_list()

    swiss_entries = {}
    annotation_lod = []

    for i in swiss_accessions[0:]:
        logger.info(i)
        swiss_entries[i] = session.get(f"https://www.uniprot.org/uniprot/{i}.json").json()
        if len(swiss_entries[i]) < 15:
            logger.info(f"  only has {len(swiss_entries[i])} entries")
        annotation_lod.append({"accession": i, "annotationScore": swiss_entries[i]['annotationScore']})

    annotation_frame = pd.DataFrame(annotation_lod)
    annotation_frame.to_sql(name="uniprot_annotation_scores", con=sqlite_db_conn, if_exists="replace", index=False)
    sqlite_db_conn.close()

    logger.info("writing to json")

    with open(swiss_entries_dump_fp, 'w') as outfile:
        json.dump(swiss_entries, outfile)

    # print("writing to yaml")
    #
    # with open(swiss_entries_dump_fp.replace('json', 'yaml'), 'w') as outfile:
    #     yaml.dump(swiss_entries, outfile)


if __name__ == "__main__":
    cli()
