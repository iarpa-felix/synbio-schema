from bx.intervals.cluster import ClusterTree
import collections
import pandas as pd
import pprint
import sqlite3

fpbase_annotation_score = 5

con = sqlite3.connect("resources/felix_dump.db")

cur = con.cursor()

df = pd.read_sql_query("""
SELECT * FROM blast_results
left join uniprot_annotation_scores on blast_results.sacc = uniprot_annotation_scores.accession
""", con)

df['index'] = df.index

df['annotationScore'] = df['annotationScore'].fillna(fpbase_annotation_score)

lod = df.to_dict('records')

cluster_trees = collections.defaultdict(lambda: ClusterTree(1, 1))

for d in lod:
    pprint.pprint(f"{d['qacc']} {d['index']}")
    qmin = min(d['qstart'], d['qend'])
    qmax = max(d['qstart'], d['qend'])
    cluster_trees[d['qacc']].insert(qmin, qmax, d['index'])

cluster_lod = []
for qseq, cluster_tree in cluster_trees.items():
    for start, end, indices in cluster_tree.getregions():
        for index in indices:
            cluster_lod.append({'qacc': qseq,
                                'index': index,
                                "cluster_id": f"{qseq}_{start}_{end}"})

cluster_frame = pd.DataFrame(cluster_lod)

blast_results_clustered_scored = df.merge(cluster_frame, on=['qacc', 'index'], how='left')

blast_results_clustered_scored['composite_score'] = blast_results_clustered_scored['qcovs'] * \
                                                    blast_results_clustered_scored['pident'] * \
                                                    blast_results_clustered_scored['annotationScore'] / 5000

blast_results_clustered_scored[f'composite_score_rank'] = blast_results_clustered_scored.groupby('cluster_id')[
    'composite_score'].rank("average", ascending=False)

blast_results_clustered_scored.to_sql(name="blast_results_clustered_scored", con=con, if_exists="replace", index=False)

blast_results_best_in_cluster = blast_results_clustered_scored.groupby('cluster_id')['composite_score_rank'].min()

blast_results_best_in_cluster = blast_results_clustered_scored.merge(blast_results_best_in_cluster,
                                                                     on=['cluster_id', 'composite_score_rank'])

blast_results_best_in_cluster.to_sql(name="blast_results_best_in_cluster", con=con, if_exists="replace", index=False)
