"""extract_nextstrain_seq_list.py
"""
from argh import dispatch_command  # type: ignore


def print_tsv_strains(nextstrainFilename: str) -> None:
    """Read a nextstrain metadata tsv and print
          all the strains it contains to stdout.
    """

    with open(nextstrainFilename, "r") as f:
        for line in f.readlines():
            print(line.split("\t")[0])


if __name__ == "__main__":
    dispatch_command(print_tsv_strains)
