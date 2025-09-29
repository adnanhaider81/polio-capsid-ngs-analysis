#!/usr/bin/env python3
# Summarize AA changes in VP1, VP2, VP3 antigenic sites relative to a reference
import argparse, os
from Bio import Entrez, SeqIO
from Bio.Seq import Seq

AG1_VP1 = (88, 106)
AG2_VP1 = (221, 226)
AG2_VP2 = (165, 174)
AG3_VP1 = (287, 292)
AG3_VP3B = (58, 60)
AG3_VP3C = (70, 76)

def fetch_rec(acc, email, api_key=None):
    Entrez.email = email
    if api_key: Entrez.api_key = api_key
    h = Entrez.efetch(db='nucleotide', id=acc, rettype='gb', retmode='text')
    rec = SeqIO.read(h, 'genbank'); h.close()
    return rec

def cds_coords(rec, gene_name):
    for f in rec.features:
        if f.type == 'CDS':
            gene = ' '.join(f.qualifiers.get('gene', [])).upper()
            product = ' '.join(f.qualifiers.get('product', [])).lower()
            if gene_name in gene or gene_name in product:
                return int(f.location.start), int(f.location.end), int(f.location.strand or 1)
    raise SystemExit(f'CDS not found for {gene_name}')

ap = argparse.ArgumentParser(description='Antigenic site AA comparison in VP1, VP2, VP3')
ap.add_argument('--email', required=False)
ap.add_argument('--api_key', required=False)
ap.add_argument('--ref_acc', required=True, help='Reference accession with annotated capsid proteins, for example V01150.1')
ap.add_argument('--genomes', required=True, help='FASTA with one or more sample genomes')
ap.add_argument('--out_tsv', required=True)
a = ap.parse_args()

email = a.email or os.getenv('NCBI_EMAIL')
if not email: raise SystemExit('Set --email or env NCBI_EMAIL')
api_key = a.api_key or os.getenv('NCBI_API_KEY')

ref = fetch_rec(a.ref_acc, email, api_key)

coords = {
    'VP1': cds_coords(ref, 'VP1'),
    'VP2': cds_coords(ref, 'VP2'),
    'VP3': cds_coords(ref, 'VP3'),
}

def nt_to_aa(seq_nt, start, end, strand):
    nt = seq_nt[start:end]
    if strand != 1:
        nt = nt.reverse_complement()
    return nt.translate()

ref_aa = {}
for gene, (s,e,strand) in coords.items():
    ref_aa[gene] = nt_to_aa(ref.seq, s, e, strand)

def slice_region(seq, start, end):
    # 1-based positions in paper, convert to 0-based indices
    return seq[start-1:end]

regions = {
    'Ag1_VP1': ('VP1', AG1_VP1),
    'Ag2_VP1': ('VP1', AG2_VP1),
    'Ag2_VP2': ('VP2', AG2_VP2),
    'Ag3_VP1': ('VP1', AG3_VP1),
    'Ag3_VP3b': ('VP3', AG3_VP3B),
    'Ag3_VP3c': ('VP3', AG3_VP3C),
}

with open(a.out_tsv, 'w') as out:
    out.write('sample	region	pos	ref	alt	gene
')
    for rec in SeqIO.parse(a.genomes, 'fasta'):
        # naive split of polyprotein into capsid ORFs using reference coords
        aa = {}
        for gene, (s,e,strand) in coords.items():
            aa[gene] = nt_to_aa(rec.seq, s, e, strand)
        for rname, (gene, (spos, epos)) in regions.items():
            ref_slice = slice_region(ref_aa[gene], spos, epos)
            samp_slice = slice_region(aa[gene], spos, epos)
            L = min(len(ref_slice), len(samp_slice))
            for i in range(L):
                if ref_slice[i] != samp_slice[i]:
                    out.write(f'{rec.id}	{rname}	{spos+i}	{ref_slice[i]}	{samp_slice[i]}	{gene}
')
