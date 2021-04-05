from Bio import SeqIO
import random
import subprocess
from tqdm import tqdm

def main():
    size = 100
    full_aln = "results/aligned-filtered.fasta"
    smaller_aln = "results/belgium/TEST_FASTA.fasta"
    records = []

    with open(full_aln, "r") as handle:
        print(f"Processing {full_aln}")
        call = ["grep", "-c", "\">\"", full_aln]
        lines = subprocess.Popen(" ".join(call), shell=True, stdout=subprocess.PIPE)
        nlines = int(lines.stdout.read().strip())
        for record in tqdm(SeqIO.parse(handle, "fasta"), desc=f"Nextstrain fasta import", total=nlines):
            n = random.random()
            if (record.id in ['Wuhan/Hu-1/2019', 'Wuhan/WH01/2019'] or n < (size/nlines)):
                records.append(record)

    print(f"Writing {smaller_aln}")
    with open(smaller_aln, "w") as output_handle:
        SeqIO.write(records, output_handle, "fasta")

def custom_subsample():
    master_fasta = "data/ALL_SEQUENCES.fasta"
    seq_list_file = "data/sequence_lists/p1.txt"
    aln_out_file = "results/belgium/TEST_FASTA_unaligned.fasta"

    with open(seq_list_file, "r") as f:
        seq_list = set([line.strip() for line in f.readlines()])

    print(f"Read {len(seq_list)} sequences to be added to fasta")

    print(f"Processing {master_fasta}")
    call = ["grep", "-c", "\">\"", master_fasta]
    lines = subprocess.Popen(" ".join(call), shell=True, stdout=subprocess.PIPE)
    nlines = int(lines.stdout.read().strip())
    with open(master_fasta, "r") as handle:
        records = []
        for record in tqdm(SeqIO.parse(handle, "fasta"), desc=f"Nextstrain fasta import", total=nlines):
            if record.id in seq_list:
                records.append(record)

    print(f"Writing {aln_out_file}")
    with open(aln_out_file, "w") as output_handle:
        SeqIO.write(records, output_handle, "fasta")


    call = ["grep", "-c", "\">\"", aln_out_file]
    lines = subprocess.Popen(" ".join(call), shell=True, stdout=subprocess.PIPE)
    nlines = int(lines.stdout.read().strip())
    print(f"Final fasta size: {nlines} sequences")



if __name__ == '__main__':
    # main()
    custom_subsample()
