# Drakenstein Nested Case-Control analysis

This repository contains the raw tables and Python entry points for the
MaAsLin2 analyses. The analysis scripts read the CSV files in `raw_data/`
directly; the loader and its dependencies are now included in the repository.

## Requirements

- Python 3.9 or newer
- R and MaAsLin2 1.26.0 (or a compatible MaAsLin2 release)
- Python packages listed in `requirements.txt`

Install the R package from Bioconductor:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install("Maaslin2")
```

Install the Python dependencies in the environment that will run the scripts:

```sh
python -m pip install -r requirements.txt
```

Make sure `Rscript` is on `PATH` and points to the R installation where
MaAsLin2 is installed. Set the `RSCRIPT` environment variable to the full
executable path if needed.

## Main case-control analysis

Run from the repository root:

```sh
python stats_ana/masslin2.py --level asv
```

Choose `asv`, `genus`, or `phylum` with `--level`. The model uses case status,
specimen collection age, and antibiotic use as fixed effects, with TSS
normalization and log transformation. MaAsLin2's full output and a
case-status result table are written under `results/Section_A/`; plots are
written alongside them. The model is called through `Rscript`, so `rpy2` is
not required. Use `--no-plots` to produce model results only. Set
`--output-dir` to choose another output location. Each run needs a fresh output
directory so existing results are preserved.

## Pathogen-stratified analyses

Run the deposited subset analyses for all listed pathogen targets and both
positive and negative strata:

```sh
python stats_ana/masslin2_subsets.py
```

Select one or more targets by repeating `--target`, for example:

```sh
python stats_ana/masslin2_subsets.py --target RSV_A_B --target-value Positive
```

Use `--target-value Negative` or `--target-value both` to choose strata, and
`--level` to choose ASV, genus, or phylum data. Outputs are written under
`results/Section_B/` by default.

## Repository files

- `raw_data/ASVtable_323sets_final.csv`: ASV count table
- `raw_data/metadata_323sets_final.csv`: specimen and clinical metadata
- `raw_data/labdata_323sets_final.csv`: laboratory measurements, including
  16S qPCR load
- `raw_data/taxdat_323sets_final.csv`: taxonomy assigned to each ASV
- `utils/load_data.py`: validates specimen and ASV identifiers, aligns tables,
  and aggregates counts to genus and phylum levels
