import pprint
import json
import sqlite3

import pandas as pd
import requests_cache
import yaml

# from typing import List

blast_results_table = 'blast_results'
pd.set_option('display.max_columns', None)
session = requests_cache.CachedSession('uniprot_entries_cache')
sqlite_db_fp = "resources/felix_dump.db"
swiss_entries_dump_fp = "resources/swiss_entries.json"

sqlite_db_conn = sqlite3.connect(sqlite_db_fp)
swiss_accessions_q = f"select distinct sacc from {blast_results_table} where blast_db = 'swissprot' and qcovs > 90 and pident > 90 order by sacc"
swiss_accessions_res = pd.read_sql_query(swiss_accessions_q, sqlite_db_conn)
swiss_accessions = swiss_accessions_res['sacc'].to_list()

swiss_entries = {}
annotation_lod = []

for i in swiss_accessions[0:]:
    print(i)
    swiss_entries[i] = session.get(f"https://www.uniprot.org/uniprot/{i}.json").json()
    annotation_lod.append({"accession": i, "annotationScore": swiss_entries[i]['annotationScore']})

annotation_frame = pd.DataFrame(annotation_lod)
annotation_frame.to_sql(name="uniprot_annotation_scores", con=sqlite_db_conn, if_exists="replace", index=False)
sqlite_db_conn.close()

# pprint.pprint(swiss_entries)

with open(swiss_entries_dump_fp, 'w') as outfile:
    json.dump(swiss_entries, outfile)

with open(swiss_entries_dump_fp.replace('json', 'yaml'), 'w') as outfile:
    yaml.dump(swiss_entries, outfile)
