import sys, os
import subprocess
from argh import dispatch_command


def grep_file(filename):
    greps = [
        "rega-7871",
        "rega-7877",
        "rega-7872",
        "rega-7884",
        "rega-7862",
        "rega-7868",
        "rega-7224",
        "rega-7216",
        "rega-7222",
        "rega-7214",
        "rega-7209",
        "rega-7212",
        "rega-7210",
        "rega-7213",
        "rega-7211",
        "rega-7208",
        "rega-7207",
        "rega-7221",
        "rega-7219",
        "rega-7225",
        "rega-7220",
        "rega-7223",
        "Aalst-OLVZ-8078888",
        "Aalst-OLVZ-8078886"
    ]

    success = 0
    fails = []
    for thing in greps:
        call = ["grep", "-c", f"\"{thing}\"", filename]
        print(" ".join(call))
        try:
            out = subprocess.check_output(" ".join(call), shell=True)
            print(f"{thing} succeeded.")
            success += 1
        except:
            print(f"{thing} not found.")
            fails.append(thing)
    if success == 24:
        print("All checks passed!")
    else:
        print("The following items were not in the dataset:")
        print(fails)

if __name__ == '__main__':
    dispatch_command(grep_file)
