"""Run MaAsLin2 case-control analyses within pathogen-positive/negative subsets."""

import argparse
from pathlib import Path
import sys

import pandas as pd

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from utils.load_data import dataset
from utils.run_maaslin import run_maaslin


CASE_TARGETS = [
    "Adenovirus",
    "Influenza_A_B_C",
    "Parainfluenza_1_2_3_4",
    "RSV_A_B",
    "Rhinovirus",
    "Cytomegalovirus",
    "S_pneumoniae",
    "H_influenzae_spp",
    "S_aureus",
    "M_catarrhalis",
    "H_influenzae_B",
    "k_pneumoniae",
]
MIN_PREVALENCE = 0.1


def run_subset(level, target, target_value, output_dir):
    """Fit one pathogen-stratified model and save its MaAsLin2 results."""
    ds = dataset()
    if level not in {"asv", "genus", "phylum"}:
        raise ValueError("level must be one of: asv, genus, phylum")
    if target not in ds.meta.columns:
        raise ValueError(f"Target {target!r} is not a metadata column")
    if target_value not in {"Positive", "Negative"}:
        raise ValueError("target_value must be Positive or Negative")

    subset_meta = ds.meta.loc[ds.meta[target].eq(target_value)].copy()
    if subset_meta.empty:
        print(f"Skipping {target}={target_value}: no matching specimens")
        return None
    features = getattr(ds, level).loc[:, subset_meta.index]
    metadata = subset_meta[
        ["age_days_npp_collection", target, "antibiotics_use_at_collection_analysis", "case"]
    ].copy()
    metadata["age_days_npp_collection"] = pd.to_numeric(
        metadata["age_days_npp_collection"], errors="raise"
    )
    metadata["case"] = pd.Categorical(metadata["case"], categories=["Case", "Control"])
    metadata["antibiotics_use_at_collection_analysis"] = pd.Categorical(
        metadata["antibiotics_use_at_collection_analysis"]
    )

    out = Path(output_dir).expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)
    raw_out = out / "maaslin2"
    if raw_out.exists():
        raise FileExistsError(
            f"MaAsLin2 output already exists: {raw_out}. Choose a new --output-dir "
            "to preserve the existing results."
        )

    fixed_effects = ["case", "age_days_npp_collection"]
    antibiotic_values = metadata["antibiotics_use_at_collection_analysis"].dropna().unique()
    reference = None
    if len(antibiotic_values) > 1:
        fixed_effects.append("antibiotics_use_at_collection_analysis")
        reference = ("antibiotics_use_at_collection_analysis", "no")

    run_maaslin(
        features=features,
        metadata=metadata,
        output=raw_out,
        fixed_effects=fixed_effects,
        min_prevalence=MIN_PREVALENCE,
        reference=reference,
    )

    results = pd.read_csv(raw_out / "all_results.tsv", sep="\t")
    case_results = results.loc[results["metadata"].eq("case")].copy()
    case_results.to_csv(out / f"{level}_case_status_results.csv", index=False)
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--level", choices=["asv", "genus", "phylum"], default="asv")
    parser.add_argument(
        "--target",
        action="append",
        choices=CASE_TARGETS,
        help="Metadata pathogen target to analyze; repeat to select several (default: all)",
    )
    parser.add_argument(
        "--target-value",
        choices=["Positive", "Negative", "both"],
        default="both",
        help="Analyze positive, negative, or both target strata (default: %(default)s)",
    )
    parser.add_argument(
        "--output-dir",
        default="results/Section_B",
        help="Directory for this run's output (default: %(default)s)",
    )
    args = parser.parse_args()

    ds = dataset()
    targets = args.target or CASE_TARGETS
    unknown = [target for target in targets if target not in ds.meta.columns]
    if unknown:
        parser.error("Metadata is missing target columns: " + ", ".join(unknown))
    values = ["Positive", "Negative"] if args.target_value == "both" else [args.target_value]

    for target in targets:
        for value in values:
            target_dir = Path(args.output_dir) / f"{target}_{value}_{args.level}"
            run_subset(args.level, target, value, target_dir)


if __name__ == "__main__":
    main()
