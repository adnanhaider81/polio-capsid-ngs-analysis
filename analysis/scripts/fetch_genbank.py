#!/usr/bin/env python3
import argparse
import os
from pathlib import Path
import time

from Bio import Entrez


def main():
    ap = argparse.ArgumentParser(description='Fetch GenBank FASTA by accession list')
    ap.add_argument('--email', required=False, help='Entrez email, or set env NCBI_EMAIL')
    ap.add_argument('--api_key', required=False, help='NCBI API key, or set env NCBI_API_KEY')
    ap.add_argument('--acc', required=True, help='Text file with one accession per line')
    ap.add_argument('--out_fasta', required=True)
    args = ap.parse_args()

    email = args.email or os.getenv('NCBI_EMAIL')
    if not email:
        raise SystemExit('Set --email or env NCBI_EMAIL to comply with NCBI usage policy')
    Entrez.email = email
    api_key = args.api_key or os.getenv('NCBI_API_KEY')
    if api_key:
        Entrez.api_key = api_key

    accs = []
    with open(args.acc) as handle:
        for line in handle:
            line = line.strip()
            if line and not line.startswith('#'):
                accs.append(line)

    out_path = Path(args.out_fasta)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open('w') as out:
        for acc in accs:
            h = Entrez.efetch(db='nucleotide', id=acc, rettype='fasta', retmode='text')
            out.write(h.read())
            h.close()
            time.sleep(0.34)
    print(f'Wrote {out_path}')

if __name__ == '__main__':
    main()
