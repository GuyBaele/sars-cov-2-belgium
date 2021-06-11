import argparse
import random
from typing import List

import pandas as pd


def subsample_lineages(lineages: List[str], number: int):

    metadataFilename = "data/metadata.tsv"

    df = pd.read_csv(metadataFilename, sep="\t")
    df = df[df["pango_lineage"].isin(lineages)]

    strains = []
    for index, row in df.iterrows():
        if row["country"] != "Belgium":
            strains.append(row["strain"])

    newStrains = random.sample(strains, number)

    return newStrains


def print_strains(strains):

    for strain in strains:
        print(strain)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--lineages", nargs="+", required=True)
    parser.add_argument("--number", type=int, required=True)
    args = parser.parse_args()

    subsampledStrains = subsample_lineages(args.lineages, args.number)

    print_strains(subsampledStrains)
