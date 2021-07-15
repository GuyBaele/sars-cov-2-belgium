"""build_datasets.py

"""
from __future__ import print_function

import datetime as dt
import json
# import multiprocessing as mp
import os
import random
import subprocess
import sys

import numpy as np
import pandas as pd
from argh import dispatch_command  # type: ignore
from Bio import SeqIO  # type: ignore
from redis_cache import cache_it  # type: ignore
from tqdm import tqdm  # type: ignore

CACHE_HOURS = 3


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

    sequence_names = read_all_sequence_lists()
    exclude_names = read_excludes()

    ##################
    #  Main process  #
    ##################
    # First, concatenate all the fasta files into one master fasta
    # This gives us two outputs:
    #   recordIDs: set of all the record names (i.e. fasta headers)
    #     that are included in the dataset
    #   records: dictionary mapping the recordIDs to their associated sequence
    #     TODO: Change this to be just sequence length
    #           since that is all we really need
    (recordIDs, records) = concat_and_write_fasta(
        gisaid_fasta, non_gisaid_dir, OUTPUT_FASTA, sequence_names, exclude_names
    )

    print(f"recordIDs: {len(recordIDs)}")
    print(f"records: {len(records)}")

    # (recordIDs, records) = bypass_fasta_prep(OUTPUT_FASTA)

    # Second, concatenate all the associated metadata
    # This is a bit of a mess
    concat_and_write_metadata(
        gisaid_metadata, non_gisaid_dir, OUTPUT_META_FNAME, recordIDs, records
    )


def bypass_fasta_prep(fastaFile):
    """Save a bunch of time during teting."""
    recIDs = set()
    recDict = {}
    with open(fastaFile, "r") as f:
        for record in tqdm(
            SeqIO.parse(f, "fasta"),
            desc="Reading fasta",
            total=count_lines_in_fasta(fastaFile),
        ):
            recIDs.add(record.id)
            recDict[record.id] = record

    return (recIDs, recDict)


# @cache_it(limit=100000, expire=60*60*CACHE_HOURS)
def concat_and_write_fasta(baseFname, fastaDir, oFname, sequence_set, exclude_set):
    """
    Take a single fasta (containing multiple GISAID records) and add a
        set of other fasta records stored in a given directory to that
        fasta. Write the output to a new file.

    Return both a set of unique record IDs and a dictionary of the records
    """

    # Initialize empty lists for outputs
    records = []
    recordIDs = set()
    duplicates = []

    @cache_it(limit=1_000_000, expire=60 * 60 * CACHE_HOURS)
    def check_record_validity(id):
        """A little helper to make check the following:

        1. The record is in our master sequence set
        2. The record is not flagged to be excluded
        """
        if id in sequence_set:
            if id in exclude_set:
                eprint(f"Excluding {id}")
                return False
            else:
                return True
        return False

    # Read the gisaid fasta
    nLines = count_lines_in_fasta(baseFname)
    print(f"Reading the base GISAID fasta: {baseFname}")
    with open(baseFname, "r") as handle:
        for record in tqdm(
            SeqIO.parse(handle, "fasta"), desc="Nextstrain fasta import", total=nLines
        ):
            # Check if a sequence with the same name already exists
            if record.id not in recordIDs:
                if check_record_validity(record.id):
                    records.append(record)
                    # Keep track of the sequence names that have been processed
                    recordIDs.add(record.id)
            else:
                # If it already exists, warn the user
                duplicates.append(f"WARNING: Duplicate record ID {record.id}.")
    print(f"Added {len(recordIDs)} records")

    process_non_gisaid_fastas(fastaDir, records, recordIDs, duplicates)

    print(f"Final dataset size (in sequences): {len(records)}")

    # Write the output fasta
    print(f"Writing {oFname}")
    with open(oFname, "w") as output_handle:
        SeqIO.write(records, output_handle, "fasta")

    # Write the list of duplicates, we care for debugging issues
    print("Writing duplicates to results/duplicate_sequence.txt")
    with open("results/duplicate_sequences.txt", "w") as output_handle:
        for line in duplicates:
            output_handle.write(f"{line}\n")

    # Transform records into a dictionary keyed on id,
    # as that will be easier to handle later
    # NOTE: This fucking sucks.
    new_records = {record.id: record for record in records}

    return recordIDs, new_records


def process_non_gisaid_fastas(fastaDir, records, recordIDs, duplicates):
    # NOTE: The following logic is more or less deprecated, as we don't really
    #         use additional fastas at this point and just pull things from
    #         GISAID instead. That said, I'm keeping it in for now.
    # TODO: Check if everything works correctly without doing this, as it will
    #         clean up the whole process quite a bit
    # Now, process each of the files in the directory that
    # contains non-gisaid fastas
    for fname in os.listdir(fastaDir):
        # Note: some of the files may be metadata files,
        # we only care about fastas for now
        if fname.endswith(".fasta"):
            print(f"Processing {fname}")
            # Keep track of how many sequences we add from additional files
            added = 0
            with open(f"{fastaDir}/{fname}", "r") as handle:
                for record in tqdm(
                    SeqIO.parse(handle, "fasta"), desc=f"Importing {fname}"
                ):
                    # Use the same logic as we did handling the gisaid fasta
                    if record.id not in recordIDs:
                        records.append(record)
                        recordIDs.add(record.id)
                        added += 1
                    else:
                        duplicates.append(f"Duplicate record ID: {record.id}.")
            print(f"Added {added} records")


def concat_and_write_metadata(baseFname, metaDir, oFname, recordIDs, records):
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
        "sample date": "date",
    }
    # `drop_cols` is the column names that we will take out of the final merge
    # Note: Maybe include "ZIP"?
    drop_cols = ["#"]

    # First, we read in the GISAID metadata file
    metadata = pd.read_csv(baseFname, sep="\t", header=0)
    print(f"Metadata original rows: {len(metadata)}")

    def reduce_metadata(df, ids):
        """Reduce a potentially massive metadata file to only include listed entries."""
        print(f"Length of metada before redution: {len(df)}.")
        newMeta = df[df["strain"].isin(list(ids))]
        newMeta = df[df["date"].apply(lambda x: len(str(x)) == 10)]
        print(f"Length of metadata after reduction: {len(newMeta)}.")

        return newMeta

    dropRows = []
    # Reduce the metadata dataFrame so that it is more reasonable to work with
    metadata = reduce_metadata(metadata, recordIDs)
    # Second, look at every file in the
    print(
        "Completed base metadata file import, now processing all excel spreadsheets in {metaDir}"
    )
    for file in tqdm(os.listdir(metaDir), desc="Reading metadata files"):
        # Only deal with excel spreadsheets for now
        if file.endswith(".xlsx"):
            # Make a new dataframe
            newMeta = pd.read_excel(f"{metaDir}/{file}", engine="openpyxl")
            # Rename the columns appropriately
            newMeta = newMeta.rename(columns=renames)
            # Slam in some "reasonable" assumptions:
            # our belgian sequences are probably european
            newMeta["region"] = "Europe"
            # they are also probably Belgian (some are French)
            newMeta["country"] = "Belgium"
            # our ncov sequences are _hopefully_ ncov
            newMeta["virus"] = "ncov"
            # full genome
            newMeta["segment"] = "genome"
            # they all come from human hosts
            newMeta["host"] = "Human"
            # We are just filling in an empty column for sequence lenght
            newMeta["length"] = np.nan
            # These aren't from GISAID, but they need a date to avoid
            # gumming up the works. We use today's date
            newMeta["date_submitted"] = dt.date.today().strftime("%Y-%m-%d")

            newMeta = reduce_metadata(newMeta, recordIDs)
            # Some things need to happen to every sequence individually
            # 1) remove year from sequence name (to match fasta)
            # 2) set country (if it isn't Belgium)
            # 3) set sequence length
            for (index, row) in tqdm(
                newMeta.iterrows(), total=len(newMeta), desc=f"Processing {file}"
            ):
                # strain name fix. I know this sucks
                try:
                    newDate = pd.to_datetime(row["date"]).strftime("%Y-%m-%d")
                    newMeta.at[index, "date"] = newDate
                    row["strain"] = fix_strain_name(row["strain"])
                    # fix country
                    row["country"] = fix_country_from_strain_name(row["strain"])
                except Exception as e:
                    with open("logs/build_datasets_warnings.log", "a") as logFile:
                        logFile.write(f"WARNING: {e}.\n")
                    dropRows.append(index)
                # set length for each sequence, if it doesn't have a length
                # for some reason indicate it should be dropped
                if row["strain"] in records.keys():
                    newMeta.at[index, "length"] = int(len(records[row["strain"]].seq))
                else:
                    dropRows.append(index)
                # I don't know why this next line exists but it seems to be necessary for things to work
                # I'll try to figure out why later if it becomes an issue
                newMeta.at[index, "date_submitted"] = newMeta.at[index, "date"]
                # determine division
            # Indicate missing data for columns for which we don't have data
            for item in set(metadata.columns).difference(set(newMeta.columns)):
                newMeta[item] = "?"

            metadata = pd.concat([metadata, newMeta])
            metadata = metadata.reset_index(drop=True)

            print(f"New metadata length: {len(metadata)}")

    # Build the big strain name-zip dictionary
    strainNameToZip = build_strain_to_zip()
    epiIslToZip = build_isl_to_zip()

    # Read the mapping fils we need
    munMap = read_muni_map()
    zipCodesToProvinces, zipCodesToMunicipalities = get_zip_location_map()
    myManualFixes = read_manual_fix_map()
    locFixes = fix_location_map()
    lonelyBoys = set()

    for (index, row) in tqdm(
        metadata.iterrows(), desc="Applying location fixes", total=len(metadata)
    ):
        # Not a location fix, but while we are looking at the individual rows
        # we also should drop anything that we don't want.
        try:
            if len(row["date"]) <= 9:
                dropRows.append(index)
                continue
        except Exception as e:
            with open("logs/build_datasets_warnings.log", "a") as logFile:
                logFile.write(f"WARNING: {e}.\n")
            dropRows.append(index)
            continue
        metadata = apply_location_corrections(
            metadata,
            index,
            row,
            epiIslToZip,
            strainNameToZip,
            locFixes,
            zipCodesToProvinces,
            zipCodesToMunicipalities,
            munMap,
            myManualFixes,
            lonelyBoys,
        )

    print("Los Lonely Boys:")
    for thing in lonelyBoys:
        print(thing)

    # Before we write, drop all the filenames
    metadata = metadata.drop(index=dropRows)
    metadata = metadata.drop(columns=drop_cols)
    # Drop duplicates
    metadata = metadata.drop_duplicates(
        subset="strain", ignore_index=True
    ).reset_index()

    # metadata = coarse_downsample(metadata)
    # print(metadata)
    print(f"Writing {oFname}")
    metadata.to_csv(oFname, sep="\t", index=False)


def spotcheck(df: pd.DataFrame, r: pd.Series, note: str) -> None:
    try:
        s = "Belgium/rega-4590/2021"
        # Europe / Belgium / Vilvoorde
        if r["strain"] == s:
            print(f"{note}: {r['location']}")
            print(df[df["strain"] == s]["location"])
    except:
        pass


def apply_location_corrections(
    metadata: pd.DataFrame,
    index: int,
    row: pd.Series,
    epiIslToZip: dict,
    strainNameToZip: dict,
    locFixes: dict,
    zipCodesToProvinces: dict,
    zipCodesToMunicipalities: dict,
    munMap: dict,
    myManualFixes: dict,
    lonelyBoys: set,
) -> pd.DataFrame:
    """NOTE: this function fucking sucks
    """
    if metadata.at[index, "country"] == "Belgium":
        if metadata.at[index, "country_exposure"] == "?":
            metadata.at[index, "country_exposure"] = "Belgium"
        # This identifies any sequences without a "location"
        if isinstance(row["location"], float):
            if row["division"] != "Belgium":
                # Trickle down
                metadata.at[index, "location"] = metadata.at[index, "division"]
                metadata.at[index, "division"] = "?"
        spotcheck(metadata, row, "1")

        # Set ZIP:
        if metadata.at[index, "gisaid_epi_isl"] in epiIslToZip.keys():
            metadata.at[index, "ZIP"] = epiIslToZip[
                metadata.at[index, "gisaid_epi_isl"]
            ]
        elif metadata.at[index, "strain"] in strainNameToZip.keys():
            metadata.at[index, "ZIP"] = strainNameToZip[metadata.at[index, "strain"]]
            spotcheck(metadata, row, "2")
        zip = str(metadata.at[index, "ZIP"])
        loc = metadata.at[index, "location"]
        # Fix location names
        if loc in locFixes.keys():
            loc = locFixes[loc]
            metadata.at[index, "location"] = loc
            spotcheck(metadata, row, "3")
        metadata.at[index, "location_exposure"] = loc
        metadata.at[index, "region_exposure"] = "Europe"

        if zip in zipCodesToProvinces.keys():
            metadata.at[index, "division"] = zipCodesToProvinces[zip]
            metadata.at[index, "division_exposure"] = zipCodesToProvinces[zip]
            metadata.at[index, "location"] = zipCodesToMunicipalities[zip]
            metadata.at[index, "location_exposure"] = zipCodesToMunicipalities[zip]
            spotcheck(metadata, row, "4")
        elif loc in munMap.keys():
            metadata.at[index, "division"] = munMap[loc]
            metadata.at[index, "division_exposure"] = munMap[loc]
            spotcheck(metadata, row, "5")
        elif loc in myManualFixes.keys():
            metadata.at[index, "division"] = myManualFixes[loc]
            metadata.at[index, "division_exposure"] = myManualFixes[loc]
            spotcheck(metadata, row, "6")
        else:
            lonelyBoys.add(loc)
            spotcheck(metadata, row, "7")

        fix_liege(metadata, index)
        spotcheck(metadata, row, "8")

    return metadata


@cache_it(limit=1000, expire=60 * 60 * CACHE_HOURS)
def read_all_sequence_lists():
    """Read all the .txt files in the sequence list directory

    This creates a sort of "master" list of all sequences that we ever might use
    """
    # Name of the directory we care about
    seqListDir = "data/sequence_lists/"
    print(f"Creating a master sequence list from {seqListDir}.")

    # Empty set to store our output
    allSeqs = set([])

    for fname in os.listdir(seqListDir):
        if fname.endswith(".txt"):
            with open(f"{seqListDir}{fname}", "r") as f:
                for line in f.readlines():
                    # Remove \n characters
                    line = line.strip()
                    allSeqs.add(line)
    print(f"Sequence list initialized with {len(allSeqs)} sequences.")
    return allSeqs


@cache_it(limit=1000, expire=60 * 60 * CACHE_HOURS)
def read_excludes():
    """Read the exclude list to give us the set of what we should ignore."""
    # Name of the file we are reading
    excludeFile = "defaults/exclude.txt"
    print(f"Creating a master exclude list from {excludeFile}")

    # Empty set to store our outuput
    exclude = set([])

    with open(excludeFile, "r") as f:
        for line in f.readlines():
            line = line.strip()
            exclude.add(line)
    print(f"Initialized exclude list with {len(exclude)} sequences.")
    return exclude


@cache_it(limit=1000, expire=60 * 60 * CACHE_HOURS)
def count_lines_in_fasta(fname):
    print(f"Processing {fname} for total fasta entries.")
    call = ["grep", "-c", '">"', fname]
    lines = subprocess.Popen(" ".join(call), shell=True, stdout=subprocess.PIPE)
    nLines = int(lines.stdout.read().strip())
    print(f"Found {nLines} fasta entries.")
    return nLines


def coarse_downsample(df):
    p = 0.0  # drop European, non-belgian sequences
    p1 = 0.0  # drop DK and UK sequences
    p2 = 0.0  # drop non-european sequences
    force_includes = read_includes()
    print(f"Started downsampling with {len(df.index)} rows.")
    drops = []
    for index, row in df.iterrows():
        if df.at[index, "country"] != "Belgium":
            n = random.random()
            if df.at[index, "strain"] not in force_includes:
                if df.at[index, "country"] in ["Denmark", "United Kingdom"]:
                    if n < p1:
                        drops.append(index)
                elif df.at[index, "region"] != "Europe":
                    if n < p2:
                        drops.append(index)
                elif n < p:
                    drops.append(index)
            if not df.at[index, "date"]:
                drops.append(index)
            elif not df.at[index, "strain"]:
                drops.append(index)
            elif not df.at[index, "date_submitted"]:
                drops.append(index)

    print(f"Attempting to remove {len(drops)} rows.")
    df = df.drop(index=drops).reset_index()  # drop the noted sequences
    print(f"Final dataset of {len(df.index)} rows.")
    return df


def read_includes():
    inclf = "defaults/include.txt"
    incl = set([])
    with open(inclf, "r") as f:
        for line in f.readlines():
            line = line.strip("\n")
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
    with open(fixfname, "r") as f:
        for line in f.readlines():
            l = line.strip("\n").split(",")
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
    df = pd.read_excel(liege_file, engine="openpyxl").rename(
        columns={"virus name": "strain", "Postal code": "ZIP"}
    )
    df2 = pd.read_csv(liege_file2).rename(columns={"sequence_ID": "strain"})

    # df = pd.concat([df,df2])
    df = pd.concat([df, df2], ignore_index=True, verify_integrity=True)

    def sf(s):
        if s.startswith("hCoV-19"):
            s = s[8:]
        return s

    for index, row in df.iterrows():
        strainName = row["strain"].strip()
        zipCode = str(row["ZIP"])
        try:
            int(zipCode.strip())
            m[strainName] = zipCode
        except Exception as e:
            # print(f"Wah: {e}")
            # print(strainName)
            pass
    return m


def build_isl_to_zip():
    r = {}
    # Add other files here
    gfile = "data/zip_codes/PostCodes_2020-12-29.xlsx"
    df = pd.concat(
        [
            pd.read_excel(gfile, sheet_name=0, engine="openpyxl"),
            pd.read_excel(gfile, sheet_name=1, engine="openpyxl"),
        ]
    )
    for i, row in df.iterrows():
        s = str(df.at[i, "GISAID_ID"])
        if s.startswith("EPI"):
            # print(s)
            # print(df.at[i])
            r[s] = str(df.at[i, "Postcode"]).strip()
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
    fixit = lambda x: x.split("(")[0][:-1] if "(" in x else x

    # TODO: Handle all these poorly caught exceptions properly
    # Set both dutch and french names
    for item in tqdm(data, desc="Reading municipalities"):
        # Add the Dutch municipality name
        try:
            map[fixit(item["TX_DESCR_NL"])] = item["PROVINCE"]
        except Exception as e:
            with open("logs/build_datasets_warnings.log", "a") as logFile:
                logFile.write(f"WARNING: {e}.\n")
        # Add the French municipality name
        try:
            map[fixit(item["TX_DESCR_FR"])] = item["PROVINCE"]
        except Exception as e:
            with open("logs/build_datasets_warnings.log", "a") as logFile:
                logFile.write(f"WARNING: {e}.\n")

    with open("data/source_files/municipalities_to_provinces.csv", "r") as f:
        for line in f.readlines():
            try:
                line = line.strip("\n").split(",")
                map[line[0]] = line[1]
            except Exception as e:
                with open("logs/build_datasets_warnings.log", "a") as logFile:
                    logFile.write("WARNING: {e}.\n")
    return map


def read_manual_fix_map():
    fname = "data/source_files/municipalities_to_provinces.csv"
    m = {}
    with open(fname, "r") as f:
        for line in f.readlines():
            try:
                line = line.strip("\n").split(",")
                k = line[0]
                v = line[1]
                m[k] = v
            except Exception as e:
                with open("logs/build_datasets_warnings.log", "a") as logFile:
                    logFile.write(f"WARNING: {e}.\n")

    return m


def fix_liege(df, i):
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
    bmap = pd.read_csv(
        "../Belgium-Geographic-Data/dist/metadata/be-dictionary.csv",
        error_bad_lines=False,
        encoding="ISO-8859-1",
    )
    bmap["PostCode"] = bmap["PostCode"].astype(int, errors="ignore")
    pro = {}
    mun = {}

    fn = {
        "Vlaams-Brabant": "VlaamsBrabant",
        "Brabant Wallon": "BrabantWallon",
        "West-Vlaanderen": "WestVlaanderen",
        "Oost-Vlaanderen": "OostVlaanderen",
        "Liège": "Liège",
    }

    myfix = lambda n: fn[n] if n in fn.keys() else n

    for index, row in bmap.iterrows():
        try:
            zip = str(int(bmap.at[index, "PostCode"]))
        except:
            continue
        if zip not in pro.keys():
            pro[zip] = myfix(bmap.at[index, "Province"])
            mun[zip] = bmap.at[index, "Municipality"]
    return pro, mun


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


if __name__ == "__main__":
    # print(build_strain_to_zip())
    # sys.exit()
    dispatch_command(main)
