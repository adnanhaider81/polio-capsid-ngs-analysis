#!/usr/bin/env python3
# Compute McNemar test and paired t test from provided CSV tables
import argparse, pandas as pd
from statsmodels.stats.contingency_tables import mcnemar
from scipy import stats

ap = argparse.ArgumentParser()
ap.add_argument('--contingency', help='CSV with columns method1, method2 where values are 0 or 1')
ap.add_argument('--paired_id', help='CSV with columns id, identity_ngs, identity_sanger')
args = ap.parse_args()

if args.contingency:
    df = pd.read_csv(args.contingency)
    b = ((df['method1']==1) & (df['method2']==0)).sum()
    c = ((df['method1']==0) & (df['method2']==1)).sum()
    res = mcnemar([[0, b],[c, 0]], exact=True)
    print('McNemar exact p-value:', res.pvalue)

if args.paired_id:
    df = pd.read_csv(args.paired_id)
    t, p = stats.ttest_rel(df['identity_ngs'], df['identity_sanger'])
    print('Paired t-test t-statistic:', t, 'p-value:', p)
