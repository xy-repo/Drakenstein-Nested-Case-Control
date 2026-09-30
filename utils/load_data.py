"""Load the raw tables deposited with this repository.

The original analysis scripts imported this module, but it was not included in
the repository.  Tables are kept in the orientation expected by those scripts:
features by rows and specimens by columns.
"""

from pathlib import Path

import pandas as pd


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATA_DIR = REPOSITORY_ROOT / "raw_data"


class dataset:
    """Provide aligned ASV, genus, phylum, metadata, and laboratory tables."""

    def __init__(self, data_dir=None):
        data_dir = Path(data_dir) if data_dir else DEFAULT_DATA_DIR

        self.asv = pd.read_csv(
            data_dir / "ASVtable_323sets_final.csv", index_col=0
        )
        self.meta = pd.read_csv(
            data_dir / "metadata_323sets_final.csv", dtype={"16S_ID": str}
        ).set_index("16S_ID")
        self.lab = pd.read_csv(
            data_dir / "labdata_323sets_final.csv", index_col=0
        )
        tax = pd.read_csv(
            data_dir / "taxdat_323sets_final.csv", index_col=0
        )

        if not self.asv.index.is_unique or not self.asv.columns.is_unique:
            raise ValueError("ASV table contains duplicate feature or specimen IDs")
        if not self.meta.index.is_unique or not self.lab.index.is_unique:
            raise ValueError("Metadata or laboratory table contains duplicate specimen IDs")

        specimen_ids = self.asv.columns
        for name, table in (("metadata", self.meta), ("laboratory", self.lab)):
            missing = specimen_ids.difference(table.index)
            if len(missing):
                raise ValueError(
                    f"{name} table is missing {len(missing)} ASV-table specimens"
                )

        missing_taxonomy = self.asv.index.difference(tax.index)
        if len(missing_taxonomy):
            raise ValueError(
                f"Taxonomy table is missing {len(missing_taxonomy)} ASVs"
            )

        self.meta = self.meta.loc[specimen_ids].copy()
        self.lab = self.lab.loc[specimen_ids].copy()
        tax = tax.loc[self.asv.index]

        self.genus = self._aggregate(self.asv, tax["Genus"])
        self.phylum = self._aggregate(self.asv, tax["Phylum"])

    @staticmethod
    def _aggregate(asv_table, labels):
        """Sum ASV counts within taxonomy labels, preserving sample order."""
        grouped = asv_table.groupby(labels.rename("taxon"), sort=False).sum()
        grouped.index.name = "taxon"
        return grouped
