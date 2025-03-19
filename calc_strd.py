#!/usr/bin/env python3

import argparse
import sys
import pandas as pd
import itertools

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input')
    parser.add_argument('-o', '--output', default=None)
    parser.add_argument('--skip_loci', nargs='+', default=None)
    return parser.parse_args()

args = parse_args()
mlva_df = pd.read_csv(args.input, index_col=0)
skip_loci = args.skip_loci if args.skip_loci else []
loci = list(mlva_df.columns)
keep_loci = [l for l in loci if l not in skip_loci]
strds = []
for i1, i2 in itertools.combinations(mlva_df.index, 2):
    r = [i1, i2]
    r1 = mlva_df.loc[i1]
    r2 = mlva_df.loc[i2]
    strd = 0
    lv = 0
    for l in keep_loci:
        l1 = r1[l]
        l2 = r2[l]
        if l1 == 'DC' or l2 == 'DC':
            continue
        try:
            l1 = int(l1)
            l2 = int(l2)
        except:
            continue
        if l1 != l2:
            lv += 1
            strd += abs(l1 - l2)
    r = [i1, i2, strd, lv]
    strds.append(r)
strd_df = pd.DataFrame(strds, columns=['Isolate1', 'Isolate2', 'STRD', 'LV'])
output = args.output if args.output else sys.stdout
strd_df.to_csv(output, index=False)
