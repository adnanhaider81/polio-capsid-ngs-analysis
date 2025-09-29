# Comparative Evaluation of Whole-Capsid NGS and Traditional GPLN Methods in Post-Culture Polio Surveillance - Pakistan and Afghanistan Study

Reproducible analysis and workflows aligned to the manuscript "Comparative Evaluation of Whole-Capsid Next-Generation Sequencing and Traditional Methods in Post-Culture Polio Surveillance: Pakistan and Afghanistan Study (2021-2024)". The repository implements the Illumina whole-capsid approach after culture and provides comparable outputs to conventional ITD plus Sanger VP1.

## Program summary
This repository contains one end to end pipeline that can reproduce the analyses reported in the manuscript. It supports the following components under one roof:
- Inputs
  - Paired end FASTQ from cultured isolates prepared for whole capsid sequencing on Illumina MiSeq 2x150.
- Quality control and trimming
  - fastp for adapter and quality trimming with default Q30 sliding window and minimum length 50.
- Mapping and primer trimming
  - BWA MEM mapping to references for WPV1 and Sabin strains.
  - BAMClipper to remove primer sequences when amplicon primers are used.
  - Picard MarkDuplicates to mark PCR duplicates.
- Coverage and consensus
  - BEDTools genomecov to compute depth profiles.
  - FreeBayes variant calling followed by BCFtools consensus with a depth mask. Default minimum depth 10.
- Antigenic site analysis
  - Translate VP1, VP2, VP3. Summarize amino acid changes in canonical antigenic sites Ag1, Ag2, Ag3. Output tidy TSV reports per sample.
- Phylogeny and context
  - Fetch selected context references from GenBank.
  - MAFFT alignment of masked consensuses plus context.
  - IQ-TREE ML tree with 1000 ultrafast bootstraps. ModelFinder can be toggled on.
- Statistics
  - Optional scripts to compute McNemar test and paired t test when you provide detection tables and per site identities from Sanger vs NGS.
- Optional de novo assembly
  - SPAdes assembly and contig QC, BLAST for sanity checking and for identifying best reference when needed.

The workflow is implemented in Snakemake and can be run on a laptop or a small server. No clinical data are included.

## Requirements
- Python 3.11 or newer
- Option A: pip and virtualenv
- Option B: conda or mamba
- Snakemake 7 or newer for the full pipeline
- System tools installed through conda: fastp, bwa, samtools, bcftools, freebayes, bedtools, picard, bamclipper, mafft, iqtree, spades, seqkit, entrez-direct

### NCBI usage note
Set a contact email once per shell for E-utilities. Optional API key improves rate limits.
```bash
export NCBI_EMAIL="you@example.com"
export NCBI_API_KEY="xxxxxxxxxxxxxxxxxxxxxxxxxxxx"   # optional
```

## Quick verification
```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r env/requirements.txt
python analysis/scripts/example_qc_plot.py --in data-example/example_counts.tsv --out results-example/example_plot.png
```

## One command run
```bash
export NCBI_EMAIL="you@example.com"
conda env create -f env/environment.yml
conda activate polio-capsid-env
snakemake -s workflow/Snakefile -c 4 --printshellcmds
```

## Configuration
Edit `config/config.yaml`. Minimal example:
```yaml
pairs:
  - sample: PV1_001
    r1: data-private/PV1_001_R1.fastq.gz
    r2: data-private/PV1_001_R2.fastq.gz

references:
  wpv1: KY941935.1
  sabin1: V01150.1
  sabin2: AY184220.1
  sabin3: AY184221.1

context_acc:
  - KY941935.1
  - V01150.1
  - AY184220.1
  - AY184221.1

params:
  threads: 4
  min_len: 50
  min_qual: 20
  min_depth_consensus: 10
  use_model_finder: false
  iqtree_model: GTR+G+I
  bootstrap: 1000
```

## Outputs
- `results/consensus/<sample>.fa` - masked consensus per sample
- `results/consensus/all_consensus.fasta` - combined consensuses
- `results/aln/wg_alignment.fasta` - alignment used for phylogeny
- `results/iqtree/wg.treefile` - ML tree
- `results/coverage/<sample>.depth.txt` - per base depth
- `results/mutations/<sample>_ag_sites.tsv` - amino acid changes in Ag1, Ag2, Ag3
- Optional de novo outputs under `results/spades/<sample>/` with QC tables and BLAST hits

## How to cite
- Software: Haider SA. Polio whole-capsid post-culture NGS analysis pipeline. Version 1.0.0. 2025. GitHub repository.

## References
- Chen S, Zhou Y, Chen Y, et al. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics. 2018.
- Li H. Aligning sequence reads with BWA-MEM. 2013. arXiv:1303.3997.
- Danecek P, et al. Twelve years of SAMtools and BCFtools. GigaScience. 2021.
- Garrison E, Marth G. Haplotype-based variant detection from short-read sequencing. arXiv:1207.3907. FreeBayes.
- Quinlan AR, Hall IM. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics. 2010.
- Picard Toolkit. Broad Institute. 2019.
- Au CH, Ho DN, Kwong A, et al. BAMClipper. Sci Rep. 2017.
- Katoh K, Standley DM. MAFFT multiple sequence alignment software version 7. Mol Biol Evol. 2013.
- Minh BQ, et al. IQ-TREE 2. Mol Biol Evol. 2020.
- Bankevich A, et al. SPAdes. J Comput Biol. 2012.
- Cock PJ, et al. Biopython. Bioinformatics. 2009.
- Kans J. Entrez Programming Utilities Help. NCBI.
- Köster J, Rahmann S. Snakemake. Bioinformatics. 2012.
- Minor PD, Ferguson M, Evans DM, et al. Antigenic structure of polioviruses. J Gen Virol. 1986.
- Hogle J, Chow M, Filman D. Three-dimensional structure of poliovirus at 2.9 Å resolution. Science. 1985.

## Contributing
See `CONTRIBUTING.md`. Open an issue for questions. Do not commit restricted data.

## License
MIT. See `LICENSE` for details.
