# Files Description

## Raw Data

1. `raw_data/ASVtable_323sets_final.csv` - Feature table
2. `raw_data/metadata_323sets_final.csv` - Metadata table
3. `raw_data/labdata_323sets_final.csv` - Labdata table
4. `raw_data/taxdat_323sets_final.csv` - Taxonomic annotations for features

## Codes

- `stats_ana/masslin2.py` - Script to run Microbiome Multivariable Associations with Linear Models (MaAsLin2) to identify ASVs which are differentially abundant between cases and controls.
- `stats_ana/masslin2_subsets.py` - Script to run MaAsLin2 on subsets comprising individuals infected with viruses or bacteria.
