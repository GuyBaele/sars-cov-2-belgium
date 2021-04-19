import argh
from Bio import SeqIO
import random
import subprocess
from tqdm import tqdm

# def main():
#     size = 100
#     full_aln = "results/aligned-filtered.fasta"
#     smaller_aln = "results/belgium/TEST_FASTA.fasta"
#     records = []
#
#     with open(full_aln, "r") as handle:
#         print(f"Processing {full_aln}")
#         call = ["grep", "-c", "\">\"", full_aln]
#         lines = subprocess.Popen(" ".join(call), shell=True, stdout=subprocess.PIPE)
#         nlines = int(lines.stdout.read().strip())
#         for record in tqdm(SeqIO.parse(handle, "fasta"), desc=f"Nextstrain fasta import", total=nlines):
#             n = random.random()
#             if (record.id in ['Wuhan/Hu-1/2019', 'Wuhan/WH01/2019'] or n < (size/nlines)):
#                 records.append(record)
    #
    # print(f"Writing {smaller_aln}")
    # with open(smaller_aln, "w") as output_handle:
    #     SeqIO.write(records, output_handle, "fasta")

def create_unaligned_fasta(build, excludes=None):
    master_fasta = "data/ALL_SEQUENCES.fasta"
    seq_list_file = f"data/sequence_lists/{build}.txt"
    aln_out_file = f"results/belgium/{build}_unaligned.fasta"
    metadata_file = "data/ALL_METADATA.tsv"

    if not excludes:
        excludes = set([])
    meta_sequences = set([])
    with open(metadata_file, "r") as f:
        for line in f.readlines():
            line = line.split("\t")
            meta_sequences.add(line[1])
            if len(line[5]) < 10: # Exclude sequences with bad dates
                excludes.add(line[1])
                print(f"{line[1]}: {line[5]}")

    with open(seq_list_file, "r") as f:
        seq_list = set([])
        for line in f.readlines():
            if line:
                if line.startswith("#"):
                    pass
                else:
                    seq_list.add(line.strip())

    print(f"Read {len(seq_list)} sequences to be added to fasta")

    print(f"Processing {master_fasta}")
    call = ["grep", "-c", "\">\"", master_fasta]
    lines = subprocess.Popen(" ".join(call), shell=True, stdout=subprocess.PIPE)
    nlines = int(lines.stdout.read().strip())
    with open(master_fasta, "r") as handle:
        records = []
        for record in tqdm(SeqIO.parse(handle, "fasta"), desc=f"Nextstrain fasta import", total=nlines):
            if record.id in seq_list and record.id not in excludes and record.id in meta_sequences:
                records.append(record)

    print(f"Writing {aln_out_file}")
    with open(aln_out_file, "w") as output_handle:
        SeqIO.write(records, output_handle, "fasta")


    call = ["grep", "-c", "\">\"", aln_out_file]
    lines = subprocess.Popen(" ".join(call), shell=True, stdout=subprocess.PIPE)
    nlines = int(lines.stdout.read().strip())
    print(f"Final fasta size: {nlines} sequences")

    return aln_out_file

def create_aligned_fasta(build, unaligned):
    '''Use nextalign to create an alignment
    '''

    call = [ "nextalign",
             "--verbose",
             f"--sequences={unaligned}",
             "--reference=data/source_files/reference.fasta",
             "--output-dir=results/belgium",
             f"--output-basename={build}"
           ]

    print("Buildng alignment with:")
    print(" ".join(call))
    subprocess.call(call)

def concat_two_fastas(base, additional, excludes=None, meta_sequences=None):
    """
    Merge one fasta into another
    """
    if not excludes:
        excludes = set([])
    if not meta_sequences:
        meta_sequences = set([])
    with open(base, "r") as b:
        records = [ record for record in SeqIO.parse(b, "fasta") ]
    with open(additional, "r") as a:
        for record in SeqIO.parse(a, "fasta"):
            if record.id not in excludes and record.id in meta_sequences:
                records.append(record)
            else:
                print(f"Excluded: {record.id}")
    with open(base, "w") as o:
        SeqIO.write(records, o, "fasta")
    return base


def create_alignment(build, extra_sequences=None):
    """
    Create an alignment based on a list of sequence names and (optionally) a secondary fasta.
    """

    with open("defaults/exclude.txt", "r") as f:
        excludes = set([line.strip() for line in f.readlines()])

    unaligned = create_unaligned_fasta(build, excludes)
    if extra_sequences:
        unaligned = concat_two_fastas(unaligned, extra_sequences, excludes)
    aligned = create_aligned_fasta(build, unaligned)


if __name__ == '__main__':
    argh.dispatch_command(create_alignment)
    # creaste_alignment(build, extra_sequences)
