import argparse
import subprocess
import pandas as pd

def main(vars, id, debug):
    all_fastas = concat_metadata(vars,id,debug)
    concat_fastas(all_fastas, debug)

def concat_metadata(vars,id,debug):
    all_fasta = []
    b = "nextmeta_base.tsv"
    l = "latest_metadata.tsv"
    for i in range(len(vars)):
        # the first time through build off the base nextmeta file
        prefix = f"BEL_{vars[i]}_gisaid_hcov-19_{id}"
        f = f"{prefix}.fasta"
        all_fasta.append(f"updated_{f}")
        p = f"{prefix}_patient_meta.tsv"
        s = f"{prefix}_sequence_meta.tsv"
        if i == 0:
            inmeta = "nextmeta_base.tsv"
            outmeta = "latest_metadata.tsv"
            call = create_cm_call(f,s,p,b,l)
            execute_call(call, debug)
        else:
            call = create_cm_call(f,s,p,l,l)
            execute_call(call,debug)

    pd.read_csv(l,sep='\t').to_excel("latest_metadata.xlsx", index=None, header=True)

    return all_fasta

def concat_fastas(fastas, debug):
    print("Concatenating fastas\n")
    if not debug:
        with open("latest_sequences.fasta", 'w') as outfile:
            for fname in fastas:
                with open(fname) as infile:
                    for line in infile:
                        if line.startswith('>'):
                            print(line[1:].strip('\n'))
                        outfile.write(line)

def create_cm_call(fasta,seq_meta,patient_meta,
                nextmeta,outfile):
    call = [
        "python",
        "convert_metadata.py",
        "--fasta", fasta,
        "--seq_meta", seq_meta,
        "--patient_meta", patient_meta,
        "--nextmeta", nextmeta,
        "--outfile", outfile
        ]
    return call

def execute_call(call, debug=False):
    print(" ".join(call))
    if not debug:
        subprocess.call(call)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--id", required=True, help="Date ID listed in gisaid downloads. E.g. 2020_01_01_01")
    parser.add_argument("--variants", nargs='*', default=["v1", "v2", "69del"], help="variants included")
    parser.add_argument("--debug", action="store_true", default=False, help="Only print calls, don't run")
    args = parser.parse_args()
    main(args.variants, args.id, args.debug)
