# SeqMetrics

[![CI](https://github.com/zavolanlab/SeqMetrics/actions/workflows/ci.yml/badge.svg)](https://github.com/zavolanlab/SeqMetrics/actions/workflows/ci.yml)

SeqMetrics computes per-sequence metrics for coding sequences, combining:

- nucleotide composition
- tissue-specific codon usage metrics such as CAI and fraction of optimal codons
- species-specific tRNA adaptation index from tRNA gene copy tables
- protein-level properties after translation


The package exposes a Python API and a `seqmetrics` command-line interface.

## Installation

Install from a local checkout:

```bash
pip install .
```

For development:

```bash
pip install -e .[test]
pytest
python -m build
```

## CLI usage

Basic analysis from a CDS FASTA file:

```bash
seqmetrics data/coding_regions.fasta > results.tsv
```

With tissue-specific codon usage and species-specific tRNA counts:

```bash
seqmetrics data/coding_regions.fasta \
  --tissue brain \
  --species H_sapiens \
  --tissue-usage-file helpers/tissue_codon_usage.tsv \
  --trna-counts helpers/hg19-tRNAs-confidence-set.out \
  --output results.tsv
```

## Python usage

```python
from seqmetrics import CodonUsageTable, SequenceAnalyzer, TRNAWeightTable
from seqmetrics.codon_definitions import SYN_CODONS_BY_AA

tissue_usage = CodonUsageTable.load_from_table(
    "helpers/tissue_codon_usage.tsv",
    syn_codons_by_aa=SYN_CODONS_BY_AA,
)

trna_table = TRNAWeightTable.from_trna_gene_file(
    species="H_sapiens",
    path="helpers/hg19-tRNAs-confidence-set.out",
)

analyzer = SequenceAnalyzer(
    tissue_usage=tissue_usage,
    trna_table=trna_table,
    default_tissue="brain",
)

rows = analyzer.analyze_fasta("data/coding_regions.fasta")
```

The resulting rows are plain dictionaries and can be passed directly to
`pandas.DataFrame` if needed by downstream analysis code.

## Output columns

Each analyzed sequence returns:

- identifiers and description
- nucleotide composition fractions (`A`, `C`, `G`, `T`)
- tissue-specific metrics (`tissue`, `tissue_cai`, `frac_opt_codons`)
- species-specific tAI (`species`, `species_tai`)
- protein properties (`length_aa`, `mw`, `pI`, `gravy`, `aromaticity`)
- secondary-structure propensity fractions (`helix_frac`, `sheet_frac`, `coil_frac`)
- amino-acid composition columns (`aa_A` through `aa_Y`) as percentages

## Reference data

The bundled helper files are based on these sources:

- Tissue-dependent codon usage table: https://dnahive.fda.gov/dna.cgi?cmd=codon_usage&id=537&mode=tisspec. PMID: 31982380. DOI: 10.1016/j.jmb.2020.01.011
- tRNA wobble-pairing efficiency parameters: https://pmc.ncbi.nlm.nih.gov/articles/PMC521650/
