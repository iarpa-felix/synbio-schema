import logging

import click
import click_log

import json

import requests
import yaml
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re

logger = logging.getLogger(__name__)
click_log.basic_config(logger)


@click.command()
@click_log.simple_verbosity_option(logger)
@click.option("--fp_url", default="https://www.fpbase.org/api/proteins/?format=json")
@click.option("--fasta_out", type=click.Path(), default="data/fpbase.fasta")
@click.option("--json_out", type=click.Path(), default="resources/fpbase.json")
@click.option("--yaml_out", type=click.Path(), default="resources/fpbase.yaml")
def cli(fp_url: str, fasta_out: str, json_out: str, yaml_out: str):

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
                    logger.warning("Error while writing sequence:  " + sr.id)


if __name__ == "__main__":
    cli()
