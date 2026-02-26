from dlhandler import unzip_fastas

import argparse
import os
from pathlib import Path

parser = argparse.ArgumentParser()

parser.add_argument("bucketsPath")

args = parser.parse_args()

bucketsPath = Path(args.bucketsPath)

items = os.listdir(bucketsPath)

for item in items:
    curPath = bucketsPath / item
    unzip_fastas(curPath)