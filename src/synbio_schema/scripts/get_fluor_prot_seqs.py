import json

import requests
import yaml
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re

fp_url = "https://www.fpbase.org/api/proteins/?format=json"
fasta_out = "data/fpbase.fasta"
json_out = "resources/fpbase.json"
yaml_out = "resources/fpbase.yaml"

fp = requests.get(fp_url)

fp_lod = fp.json()

fp_dict = {}
for fp in fp_lod:
    fp_dict[fp['slug']] = fp

with open(json_out, 'w') as outfile:
    json.dump(fp_dict, outfile)

with open(yaml_out, 'w') as outfile:
    yaml.dump(fp_dict, outfile)

with open(fasta_out, "w") as f_out:
    for seqs in fp_lod:
        if seqs['slug'] and seqs['seq']:
            tidy = re.sub(r'\s+', '', seqs["seq"])
            sr = SeqRecord(Seq(tidy), str(seqs["slug"]), "", "")
            r = SeqIO.write(sr, f_out, "fasta")
            if r != 1:
                print("Error while writing sequence:  " + sr.id)
