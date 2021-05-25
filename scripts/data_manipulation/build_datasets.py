import os
import subprocess
import datetime as dt
import numpy as np
import json
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm
import random


def main():
    """The main process to follow for incorporating metadata files
    """
    # There should be two files that come from GISAID:
    #   1) A dated metadata tsv
    #   2) A dated sequences fasta
    # These files can be found throught:
    # GISAID
    #   -> EpiCoV
    #     -> Downloads
    #        -> Genomic epidemiology
    #          -> "FASTA" and "metadata" links
    # After being downloaded and extracted with `gunzip`
    # they can be renamed/relocated to the paths shown below
    gisaid_metadata = "data/metadata.tsv"
    gisaid_fasta = "data/sequences.fasta"

    # We expect to have a directory full of data (both sequence and metadata)
    # which is not on GISAID
    non_gisaid_dir = "data/non_gisaid"

    # Define names of the updated sequence and metadata files
    # that we want at the end of the pipeline
    OUTPUT_FASTA = "data/ALL_SEQUENCES.fasta"
    OUTPUT_META_FNAME = "data/ALL_METADATA.tsv"

    # We will take a set of other metadata files from other sources as well
    # they
    # should
    # go
    # here

    ##################
    #  Main process  #
    ##################
    # First, concatenate all the fasta files into one master fasta
    # This gives us two outputs:
    #   record_ids: set of all the record names (i.e. fasta headers)
    #     that are included in the dataset
    #   records: dictionary mapping the record_ids to their associated sequence
    #     TODO: Change this to be just sequence length, since that is all we need
    (record_ids, records) = concat_and_write_fasta(gisaid_fasta,
                                                   non_gisaid_dir,
                                                   OUTPUT_FASTA)

    # Second, concatenate all the associated metadata
    # This is a bit of a mess
    concat_and_write_metadata(gisaid_metadata,
                              non_gisaid_dir,
                              OUTPUT_META_FNAME,
                              record_ids, records)


def concat_and_write_fasta(base_fname, fasta_dir, o_fname):
    """
    Take a single fasta (containing multiple GISAID records) and add a set of other fasta records
        stored in a given directory to that fasta. Write the output to a new file.

    Return both a set of unique record IDs and a dictionary of the records
    """
    # Initialize empty lists
    records = []
    record_ids = set()
    duplicates = []
    # Read the gisaid fasta
    with open(base_fname, "r") as handle:
        print(f"Processing {base_fname}")
        call = ["grep", "-c", "\">\"", base_fname]
        lines = subprocess.Popen(" ".join(call), shell=True, stdout=subprocess.PIPE)
        nlines = int(lines.stdout.read().strip())
        for record in tqdm(SeqIO.parse(handle, "fasta"), desc="Nextstrain fasta import", total=nlines):
            # Check if a sequence with the same name already exists in the dataset
            if record.id not in record_ids:
                # If not add the record to the master group of records
                records.append(record)
                # Also keep track of the sequence names that have been processed
                record_ids.add(record.id)
            else:
                # If it already exists, warn the user
                duplicates.append(f"WARNING: Duplicate record ID {record.id}, skipping.")
    print(f"Added {len(record_ids)} records")

    # Now, process each of the files in the directory that contains non-gisaid fasta files
    for fname in os.listdir(fasta_dir):
        # Note: some of the files may be metadata files, we only care about fastas now
        if fname.endswith(".fasta"):
            print(f"Processing {fname}")
            # Keep track of how many new sequences were added from the file for debugging
            added = 0
            with open(f"{fasta_dir}/{fname}", "r") as handle:
                for record in tqdm(SeqIO.parse(handle, "fasta"), desc=f"Importing {fname}"):
                    # Use the same logic as we did handling the gisaid fasta above
                    if record.id not in record_ids:
                        records.append(record)
                        record_ids.add(record.id)
                        added += 1
                    else:
                        duplicates.append(f"Duplicate record ID: {record.id}, skipping.")
            print(f"Added {added} records")

    # Write the output fasta
    print(f"Writing {o_fname}")
    with open(o_fname, "w") as output_handle:
        SeqIO.write(records, output_handle, "fasta")
    print(f"Writing duplicates to results/duplicate_sequence.txt")
    with open("results/duplicate_sequences.txt", "w") as output_handle:
        for line in duplicates:
            output_handle.write(f"{line}\n")

    # Transform records into a dictionary keyed on id, as that will be easier to handle later
    transform_records = lambda records: { record.id: record for record in records }
    records = transform_records(records) # probably an abuse of names to rename records to a new datatype, but whatever

    return record_ids, records

def concat_and_write_metadata(base_fname, meta_dir, o_fname, record_ids, records):
    """
    IMPORTANT: This function absolutely sucks. I'll try to
    break it apart into more sub-functions with time

    This function takes multiple metadata spreadsheets and
    sticks them together.
    Then it appropriately fixes Belgian samples so they
    behave the way that we want them to.
    It then writes to a new file.

    Also it prints out a bunch of stuff that needs to be fixed later on
    """

    # Define some things that will be used later.
    # There is some inconsistency in how headers are labeled,
    # `renames` maps between those
    renames = {
                "sequence name": "strain",
                "Town": "location",
                "Acc.Number": "gisaid_epi_isl",
                "Sex": "sex",
                "Age": "age",
                "sample date": "date"
               }
    # `drop_cols` is the column names that we will take out of the final merge
    # Note: Maybe include "ZIP"?
    drop_cols = ["#"]

    # First, we read in the GISAID metadata file
    metadata = pd.read_csv(base_fname, sep='\t', header=0)
    print(f"Metadata original rows: {len(metadata)}")

    # Second, look at every file in the
    for file in tqdm(os.listdir(meta_dir), desc="Reading metadata files"):
        # Only deal with excel spreadsheets for now
        if file.endswith(".xlsx"):
            # Make a new dataframe
            new_meta = pd.read_excel(f"{meta_dir}/{file}", engine='openpyxl')
            # Rename the columns appropriately so that the stuff we want matches
            new_meta = new_meta.rename(columns=renames)
            # Slam in some "reasonable" assumptions:
            new_meta["region"] = "Europe" # our belgian sequences are probably european
            new_meta["country"] = "Belgium" # they are also probably Belgian (some are French; we deal with that later)
            new_meta["virus"] = "ncov" # our ncov sequences are _hopefully_ ncov
            new_meta["segment"] = "genome" # full genome
            new_meta["host"] = "Human" # they all come from human hosts
            new_meta["length"] = np.nan # We are just filling in an empty column for sequence lenght, dealt with later
            # These aren't from GISAID, but they need a date to avoid gumming up the works. We use today's date
            new_meta["date_submitted"] = dt.date.today().strftime("%Y-%m-%d")

            # Some things need to happen to every row (i.e. new sequence) individually
            # 1) remove year from sequence name (to match fasta) TODO: make this smarter
            # 2) set country (if it isn't Belgium)
            # 3) set sequence length
            drop_rows = []
            for (index, row) in tqdm(new_meta.iterrows(), total=len(new_meta)):
                # strain name fix. I know this sucks
                try:
                    new_meta.at[index,"date"] = pd.to_datetime(row["date"]).strftime("%Y-%m-%d")
                    row["strain"] = fix_strain_name(row["strain"])
                    # fix country
                    row["country"] = fix_country_from_strain_name(row["strain"])
                except:
                    drop_rows.append(index)
                # set length for each sequence, if it doesn't have a length for some reason indicate it should be dropped
                if row["strain"] in records.keys():
                    new_meta.at[index,"length"] = int(len(records[row["strain"]].seq))
                else:
                    drop_rows.append(index)
                # I don't know why this next line exists but it seems to be necessary for things to work
                # I'll try to figure out why later if it becomes an issue
                new_meta.at[index,"date_submitted"] = new_meta.at[index,"date"]
                # determine division
            # Indicate missing data for columns for which we don't have data
            for item in set(metadata.columns).difference(set(new_meta.columns)):
                new_meta[item] = "?"

            metadata = pd.concat([metadata, new_meta])
            metadata = metadata.reset_index(drop=True)

            print(f"New metadata length: {len(metadata)}")


    # Build the big strain name-zip dictionary
    strain_zip = build_strain_to_zip()
    epi_isl_zip = build_isl_to_zip()

    # Read the mapping fils we need
    mmap = read_muni_map()
    zpro,zmun = get_zip_location_map()
    my_manual_fixes = read_manual_fix_map()
    loc_fixes = fix_location_map()
    lonelyboys = set()

    c = 0
    for (index, row) in metadata.iterrows():
        try:
            if len(row["date"]) <= 9:
                drop_rows.append(index)
        except:
            drop_rows.append(index)
        if metadata.at[index,"country"] == "Belgium":
            if metadata.at[index,"country_exposure"] == "?":
                metadata.at[index,"country_exposure"] = "Belgium"
            # This identifies any sequences without a "location"
            if isinstance(row["location"],float):
                if row["division"] != "Belgium":
                    # Trickle down
                    metadata.at[index,"location"] = metadata.at[index,"division"]
                    metadata.at[index,"division"] = "?"

            # Set ZIP:
            if metadata.at[index,"gisaid_epi_isl"] in epi_isl_zip.keys():
                metadata.at[index,"ZIP"] = epi_isl_zip[metadata.at[index,"gisaid_epi_isl"]]
            elif metadata.at[index,"strain"] in strain_zip.keys():
                metadata.at[index,"ZIP"] = strain_zip[metadata.at[index,"strain"]]
            zip = str(metadata.at[index,"ZIP"])
            loc = metadata.at[index,"location"]
            # Fix location names
            if loc in loc_fixes.keys():
                loc = loc_fixes[loc]
                metadata.at[index,"location"] = loc

            metadata.at[index,"location_exposure"] = loc
            metadata.at[index,"region_exposure"] = "Europe"


            if zip in zpro.keys():
                metadata.at[index,"division"] = zpro[zip]
                metadata.at[index,"division_exposure"] = zpro[zip]
                metadata.at[index,"location"] = zmun[zip]
                metadata.at[index,"location_exposure"] = zmun[zip]
                c += 1
            elif loc in mmap.keys():
                metadata.at[index,"division"] = mmap[loc]
                metadata.at[index,"division_exposure"] = mmap[loc]
            elif loc in my_manual_fixes.keys():
                metadata.at[index,"division"] = my_manual_fixes[loc]
                metadata.at[index,"division_exposure"] = my_manual_fixes[loc]
            else:
                lonelyboys.add(loc)

            fix_liege(metadata,index)

    print(f"Set {c} locations based on ZIP")

    print("Los Lonely Boys:")
    for thing in lonelyboys:
        print(thing)

    # Before we write, drop all the filenames
    metadata = metadata.drop(index=drop_rows)
    metadata = metadata.drop(columns=drop_cols)
    # Drop duplicates
    metadata = metadata.drop_duplicates(subset="strain", ignore_index=True).reset_index()

    # metadata = coarse_downsample(metadata)
    # print(metadata)


    print(f"Writing {o_fname}")
    metadata.to_csv(o_fname, sep='\t', index=False)

def coarse_downsample(df):
    p=0.0 # drop European, non-belgian sequences
    p1=0.0 # drop DK and UK sequences
    p2=0.0 # drop non-european sequences
    force_includes = read_includes()
    print(f"Started downsampling with {len(df.index)} rows.")
    drops = []
    for index,row in df.iterrows():
        if df.at[index,"country"] != "Belgium":
            n = random.random()
            if df.at[index,"strain"] not in force_includes:
                if df.at[index,"country"] in ["Denmark", "United Kingdom"]:
                    if (n<p1):
                        drops.append(index)
                elif df.at[index,"region"] != "Europe":
                    if n<p2:
                        drops.append(index)
                elif (n < p):
                    drops.append(index)
            if not df.at[index,"date"]:
                drops.append(index)
            elif not df.at[index,"strain"]:
                drops.append(index)
            elif not df.at[index,"date_submitted"]:
                drops.append(index)

    print(f"Attempting to remove {len(drops)} rows.")
    df = df.drop(index=drops).reset_index() # drop the noted sequences
    print(f"Final dataset of {len(df.index)} rows.")
    return df

def read_includes():
    inclf = "defaults/include.txt"
    incl = set([])
    with open(inclf,'r') as f:
        for line in f.readlines():
            line=line.strip('\n')
            incl.add(line)
    return incl

def fix_strain_name(s):
    """
    This can be expanded later if we need it
    """
    # Cast to str
    s = str(s)
    # Remove trailing dates from strain names
    if s.endswith("/2020") or s.endswith("/2019"):
        s = "/".join(s.split("/")[:-1])
    # Remove leading ""
    return s

def fix_location_map():
    m = {}
    fixfname = "data/source_files/municipalities_name_fixes.csv"
    with open(fixfname,'r') as f:
        for line in f.readlines():
            l = line.strip('\n').split(',')
            k = l[0]
            v = l[1]
            m[k] = v
    return m

def fix_country_from_strain_name(s):
    """
    Pull a country from a strain name
    """
    c = s.split("/")[1]

    return c

def build_strain_to_zip():
    m = {}
    liege_file = "data/zip_codes/SARS-CoV-2_ULiegeSeq_211220.xlsx"
    liege_file2 = "data/zip_codes/SARS-CoV-2_ULiegeSeq_011220.csv"
    df = pd.read_excel(liege_file, engine='openpyxl').rename(columns={"virus name": "strain","Postal code": "ZIP"})
    df2 = pd.read_csv(liege_file2).rename(columns={"sequence_ID": "strain"})

    # df = pd.concat([df,df2])
    df = pd.concat([df,df2],ignore_index=True,verify_integrity=True)

    def sf(s):
        if s.startswith("hCoV-19"):
            s = s[8:]
        return s

    for i,r in df.iterrows():
        k = df.at[i,"strain"]
        v = df.at[i,"ZIP"]
        k = sf(str(k)).strip()
        try:
            int(str(v.strip()))
            m[k] = str(v)
        except:
            pass
    return m

def build_isl_to_zip():
    r = {}
    # Add other files here
    gfile = "data/zip_codes/PostCodes_2020-12-29.xlsx"
    df = pd.concat([pd.read_excel(gfile,sheet_name=0, engine='openpyxl'),pd.read_excel(gfile,sheet_name=1, engine='openpyxl')])
    for i,row in df.iterrows():
        s = str(df.at[i,"GISAID_ID"])
        if s.startswith("EPI"):
            # print(s)
            # print(df.at[i])
            r[s] = str(df.at[i,"Postcode"]).strip()
    return r


def read_muni_map(case_json="data/epi/COVID19BE_CASES_MUNI_CUM.json"):
    """Parse a set of files mapping municipality names to their province.

    Keyword arguments:
    case_json -- a string that gives the path (relative to project root)
                     of the case json file that is being read

    Output:
    map -- a dictionary keying all named municipalities in Belgium to their
               province. Each municipality will be the key two times:
               once in Dutch and once in French.

    {"Leuven" : "VlaamsBrabant",
     "Louvain" : "VlaamsBrabant",
      ...
    }
    """

    print("Creating a map of municipalities to their province.")
    print("This may take several minutes.")

    map = {}  # Initialize the final dictionary to be returned

    # Use the json module to load the full file into memory :(
    with open(case_json, "r") as f:
        data = json.load(f)

    # Add a small function that will clean up municipalities with parentheses
    fixit = lambda x: x.split("(")[0][:-1] if '(' in x else x

    # TODO: Handle all these poorly caught exceptions properly
    # Set both dutch and french names
    for item in tqdm(data, desc="Reading municipalities"):
        # Add the Dutch municipality name
        try:
            map[fixit(item["TX_DESCR_NL"])] = item["PROVINCE"]
        except Exception as e:
            print(f"WARNING: {e}")
            pass
        # Add the French municipality name
        try:
            map[fixit(item["TX_DESCR_FR"])] = item["PROVINCE"]
        except Exception as e:
            print(f"WARNING: {e}")
            pass

    with open("data/source_files/municipalities_to_provinces.csv", 'r') as f:
        for line in f.readlines():
            try:
                line = line.strip('\n').split(',')
                map[line[0]] = line[1]
            except Exception as e:
                print(f"WARNING: {e}")
                pass
    return map


def read_manual_fix_map():
    fname = "data/source_files/municipalities_to_provinces.csv"
    m = {}
    with open(fname, "r") as f:
        for line in f.readlines():
            try:
                line = line.strip('\n').split(',')
                k = line[0]
                v = line[1]
                m[k] = v
            except Exception as e:
                print(f"WARNING: {e}")
                pass
    return m


def fix_liege(df,i):
    """
    Add diacritic marks to Liège
    """
    geo_fixes = ["location", "location_exposure", "division", "division_exposure"]
    for gf in geo_fixes:
        if df.at[i, gf] == "Liege":
            df.at[i, gf] = "Liège"


def get_zip_location_map():
    """make dictionaries taking zip code to province and municipality
    """
    bmap = pd.read_csv("../Belgium-Geographic-Data/dist/metadata/be-dictionary.csv",
                       error_bad_lines=False, encoding="ISO-8859-1")
    bmap["PostCode"] = bmap["PostCode"].astype(int, errors='ignore')
    pro = {}
    mun = {}

    fn = {"Vlaams-Brabant": "VlaamsBrabant",
          "Brabant Wallon": "BrabantWallon",
          "West-Vlaanderen": "WestVlaanderen",
          "Oost-Vlaanderen": "OostVlaanderen",
          "Liège": "Liège"}

    myfix = lambda n: fn[n] if n in fn.keys() else n

    for index, row in bmap.iterrows():
        zip = str(bmap.at[index, "PostCode"])
        if zip not in pro.keys():
            pro[zip] = myfix(bmap.at[index, "Province"])
            mun[zip] = bmap.at[index, "Municipality"]
    return pro, mun


if __name__ == "__main__":
    # print(build_strain_to_zip())
    # sys.exit()
    main()
