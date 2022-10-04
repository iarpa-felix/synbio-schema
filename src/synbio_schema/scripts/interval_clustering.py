from bx.intervals.cluster import ClusterTree
import collections
import pandas as pd
import pprint
import sqlite3

con = sqlite3.connect("resources/felix_dump.db")

cur = con.cursor()

df = pd.read_sql_query("""
SELECT * FROM blast_results
join uniprot_annotation_scores on blast_results.sacc = uniprot_annotation_scores.accession
""", con)

df['index'] = df.index

lod = df.to_dict('records')

cluster_trees = collections.defaultdict(lambda: ClusterTree(1, 1))

for d in lod:
    pprint.pprint(d)
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

merged_cluster_frame = df.merge(cluster_frame, on=['qacc', 'index'], how='left')

merged_cluster_frame['composite_score'] = merged_cluster_frame['qcovs'] * merged_cluster_frame['pident'] * \
                                          merged_cluster_frame['annotationScore'] / 5000

merged_cluster_frame[f'composite_score_rank'] = merged_cluster_frame.groupby('cluster_id')[
    'composite_score'].rank("average", ascending=False)

merged_cluster_frame.to_sql(name="blast_results_clustered_scored", con=con, if_exists="replace", index=False)

min_rank_by_cluster = merged_cluster_frame.groupby('cluster_id')['composite_score_rank'].min()

min_rank_by_cluster = merged_cluster_frame.merge(min_rank_by_cluster, on=['cluster_id', 'composite_score_rank'])

# print(min_rank_by_cluster)

min_rank_by_cluster.to_sql(name="blast_results_best_in_cluster", con=con, if_exists="replace", index=False)
