import argh
import sys, os

def rename_files(filestem):
    suffixes = ["", "_tip-frequencies", "_root-sequence"]

    for suffix in suffixes:
        old_fname = f"auspice/ncov_belgium{suffix}.json"
        new_fname = f"auspice/sars-cov-2-belgium_{filestem}{suffix}.json"
        print(f"mv {old_fname} {new_fname}")
        os.rename(old_fname, new_fname)

if __name__ == '__main__':
    argh.dispatch_command(rename_files)
