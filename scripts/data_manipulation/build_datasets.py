import sys,os
import datetime as dt
import numpy as np
import pandas as pd
from Bio import SeqIO

def concat_fasta(base_fname, fasta_dir, o_fname):
    records = []
    record_ids = set()
    with open(base_fname, "r") as handle:
        print(f"Processing {base_fname}")
        for record in SeqIO.parse(handle, "fasta"):
            if record.id not in record_ids:
                records.append(record)
                record_ids.add(record.id)
            else:
                print(f"Duplicate record ID: {record.id}, skipping.")
    print(f"Added {len(record_ids)} records")
    for fname in os.listdir(fasta_dir):
        if fname.endswith(".fasta"):
            print(f"Processing {fname}")
            added = 0
            with open(f"{fasta_dir}/{fname}", "r") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    if record.id not in record_ids:
                        records.append(record)
                        record_ids.add(record.id)
                        added += 1
                    else:
                        print(f"Duplicate record ID: {record.id}, skipping.")
            print(f"Added {added} records")
    print(f"Writing {o_fname}")
    with open(o_fname, "w") as output_handle:
        SeqIO.write(records, output_handle, "fasta")

    # Transform records into a dictionary keyed on id
    def transform_records(records):
        return { record.id: record for record in records }

    records = transform_records(records)

    return record_ids, records

def concat_metadata(base_fname, meta_dir, o_fname, record_ids, records):
    metadata = pd.read_csv(base_fname, sep='\t', header=0)
    print(f"Metadata original rows: {len(metadata)}")
    renames = {
                "sequence name" : "strain",
                "Town" : "location",
                "Acc.Number" : "gisaid_epi_isl",
                "Sex" : "sex",
                "Age" : "age",
                "sample date": "date"
               }
    drop_cols = ["#", "ZIP"]

    for file in os.listdir(meta_dir):
        if file.endswith(".xlsx"):
            new_meta = pd.read_excel(f"{meta_dir}/{file}")
            new_meta = new_meta.rename(columns=renames)
            new_meta = new_meta.drop(columns=drop_cols)
            new_meta["region"] = "Europe"
            new_meta["country"] = "Belgium"
            new_meta["virus"] = "ncov"
            new_meta["segment"] = "genome"
            new_meta["host"] = "Human"
            new_meta["length"] = np.nan
            new_meta["date_submitted"] = "2020-12-31"
            # Some things need to happen to every row individually
            # 1) remove year from sequence name (to match fasta) TODO: make this smarter
            # 2) set country (if it isn't Belgium)
            # 3) set sequence length
            drop_rows = []
            for (index, row) in new_meta.iterrows():
                # strain name fix. I know this sucks
                new_meta.at[index,"date"] = pd.to_datetime(row["date"]).strftime("%Y-%m-%d")
                s = row["strain"]
                if s.endswith("/2020") or s.endswith("/2019"):
                    s = "/".join(s.split("/")[:-1])
                row["strain"] = s
                # set country
                c = s.split("/")[1]
                if c != "Belgium":
                    row["Country"] = c
                    # print(f"Updated {s}'s country to {c}")
                # set length
                if row["strain"] in records.keys():
                    new_meta.at[index,"length"] = int(len(records[row["strain"]].seq))
                else:
                    drop_rows.append(index)
                new_meta.at[index,"date_submitted"] = new_meta.at[index,"date"]
            new_meta = new_meta.drop(index=drop_rows)
            for item in set(metadata.columns).difference(set(new_meta.columns)):
                new_meta[item] = "?"
            print(new_meta)
            metadata = pd.concat([metadata, new_meta])
            print(f"New metadata length: {len(metadata)}")
    print(metadata)
    print(f"Writing {o_fname}")
    metadata.to_csv(o_fname, sep='\t', index=False)

def main():
    FASTA_BASE = "data/sequences_2020-12-04_07-42.fasta"
    FASTA_DIR = "data/raw"
    OUTPUT_FASTA = "data/ALL_SEQUENCES.fasta"

    METADATA_BASE = "data/metadata_2020-12-05_11-42.tsv"
    METADATA_DIR = "data/raw"
    OUTPUT_META_FNAME = "data/ALL_METADATA.tsv"

    def make_dummy():
        ids = set()
        records = {}
        with open(OUTPUT_FASTA, "r") as i:
            for record in SeqIO.parse(i, "fasta"):
                ids.add(record.id)
                records[record.id] = record
        return ids, records

    (record_ids, records) = concat_fasta(FASTA_BASE, FASTA_DIR, OUTPUT_FASTA)
    # (record_ids, records) = make_dummy() # for fast debugging
    concat_metadata(METADATA_BASE, METADATA_DIR, OUTPUT_META_FNAME, record_ids, records)

if __name__ == '__main__':
    main()
