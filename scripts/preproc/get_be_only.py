import argparse
from Bio import SeqIO
import pandas as pd

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sequences", required=True, help="FASTA file of sequences")
    parser.add_argument("--metadata", required=True, type=str, help="nextmeta")
    parser.add_argument("--output", required=True, help="outfile")

    args = parser.parse_args()

    bel = pd.read_csv(args.metadata,sep='\t').query("country == 'Belgium'")['strain'].to_list()
    bel_seq=[]
    for r in SeqIO.parse(args.sequences,'fasta'):
    	if r.id in bel:
    		bel_seq.append(r)
    SeqIO.write(bel_seq,args.output,'fasta')

