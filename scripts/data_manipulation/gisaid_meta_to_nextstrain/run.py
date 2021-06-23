"""run.py

This module serves as a wrapper for convert_metadata.py,
calling it for each variant in a user-defined list.

Example:
    If we want to process only B.1.617 files, we simply run:

    $ python run.py --variantList B.1.617
"""
import argparse  # TODO: convert from argparse to argh
import subprocess
from typing import List

import pandas as pd
from tqdm import tqdm  # type: ignore


def main(variantList: List[str], debug: bool) -> None:
    """Call concatenate_metadata and concatenate_fastas.
    """
    all_fastas = concatenate_metadata(variantList, debug)
    concatenate_fastas(all_fastas, debug)


def concatenate_metadata(variantList: List[str], debug: bool) -> List[str]:
    """For each desired variant, call `convert_metadata.py`

    Args:
        variantList: all the desired varants for which there are downloaded
                        GISAID metadata files
        debug: flag for whether commands should be printed or executed

    Returns:
        A list of all the modified fasta filenames.
    """
    allFastas = []  # initialize empty list to store return values
    baseFname = "nextmeta_base.tsv"  # tsv with the expected headers in order
    outputFname = "latest_metadata.tsv"  # output file
    for i in tqdm(range(len(variantList))):
        prefix = variantList[i]
        fastaFname = f"{prefix}.fasta"
        patientMetadataFname = f"{prefix}_patient_meta.tsv"
        sequenceMetadataFname = f"{prefix}_sequence_meta.tsv"
        if i == 0:  # the first time through build off the base nextMeta file
            call = create_cm_call(
                fastaFname,
                sequenceMetadataFname,
                patientMetadataFname,
                baseFname,
                outputFname,
            )
            execute_call(call, debug)
        else:
            call = create_cm_call(
                fastaFname,
                sequenceMetadataFname,
                patientMetadataFname,
                outputFname,
                outputFname,
            )
            execute_call(call, debug)
        # After the main functionality has finished, track which
        # files have been processed so far.
        # NOTE: the "updated_" prefix is added as part of the execution
        #       of `convert_metadata.py`
        # TODO: make sure this isn't a lie
        allFastas.append(f"updated_{fastaFname}")

    # The final step here is to convert the tsv we just created to an excel
    # formatted spreadsheet, so that it plays nicely with downstream steps
    pd.read_csv(outputFname, sep="\t").to_excel(
        "latest_metadata.xlsx", index=None, header=True
    )

    return allFastas


def concatenate_fastas(fastaList: List[str], debug: bool):
    """Iterate combine all fastas in fastaList into one."""
    print(f"Concatenating {fastaList}\n")
    outFileName = "latest_sequences.fasta"
    if not debug:
        with open(outFileName, "w") as outFileHandle:
            for fname in fastaList:
                with open(fname) as infile:
                    for line in infile:
                        if line.startswith(">"):
                            print(line[1:].strip("\n"))
                        outFileHandle.write(line)


def create_cm_call(
    fasta: str, seqMeta: str, patientMeta: str, nextMeta: str, outFile: str
) -> List[str]:
    """Create bash-like call to execute convert_metadata.py

    Args:
        fasta: TODO: finish these annotations after
                        rereading `convert_metadata.py`
        seqMeta:
        patientMeta:
        nextMeta:
        outFile:

    Returns:
        A a list that encodes a bash-like call that can be fed into
        the subprocess.call() command to be executed.

        E.g:
            $ ls -lah
        would be encoded as:
            ["ls", "-lah"]

        Here, the encoded call is:
            $ python convert_metadata.py ARGS
    """
    call = [
        "python",
        "convert_metadata.py",
        "--fasta",
        fasta,
        "--seqMeta",
        seqMeta,
        "--patientMeta",
        patientMeta,
        "--nextMeta",
        nextMeta,
        "--outFile",
        outFile,
    ]

    return call


def execute_call(call: List[str], debug: bool = False) -> None:
    """Print a call to stdout and execute that call if not in debug mode.

    Args:
        call: A bash-like command represented as a list baseFname
        debug: Flag for debug mode; True => the call will
            only be printed, not executed
    """
    print(" ".join(call))
    if not debug:
        subprocess.call(call)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--variants",
        nargs="*",
        default=["P.1", "B.1.617"],
        help="variants to be added to the analysis",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        default=False,
        help="Only print calls, don't run",
    )
    args = parser.parse_args()
    main(args.variants, args.debug)

# Local Variables:
# conda-env-name-for-buffer: nextstrain
# End:
