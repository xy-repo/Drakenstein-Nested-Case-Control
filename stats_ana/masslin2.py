"""Run the deposited case-control MaAsLin2 analysis."""

import argparse
from pathlib import Path
import sys

import pandas as pd

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from utils.load_data import dataset
from utils.run_maaslin import run_maaslin


def run(level="asv", output_dir="results/Section_A/Q6_Masslin2_case_control", plots=True):
    """Fit the case-control model and save MaAsLin2 output and plots."""
    ds = dataset()
    if level not in {"asv", "genus", "phylum"}:
        raise ValueError("level must be one of: asv, genus, phylum")

    features = getattr(ds, level).copy()
    metadata = ds.meta[
        ["case", "age_days_npp_collection", "antibiotics_use_at_collection_analysis"]
    ].copy()
    metadata["case"] = pd.Categorical(metadata["case"], categories=["Case", "Control"])
    metadata["age_days_npp_collection"] = pd.to_numeric(
        metadata["age_days_npp_collection"], errors="raise"
    )
    metadata["antibiotics_use_at_collection_analysis"] = pd.Categorical(
        metadata["antibiotics_use_at_collection_analysis"]
    )

    out = Path(output_dir).expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)
    raw_out = out / f"{level}_maaslin2"
    if raw_out.exists():
        raise FileExistsError(
            f"MaAsLin2 output already exists: {raw_out}. Choose a new --output-dir "
            "to preserve the existing results."
        )

    run_maaslin(
        features=features,
        metadata=metadata,
        output=raw_out,
        fixed_effects=[
            "case",
            "age_days_npp_collection",
            "antibiotics_use_at_collection_analysis",
        ],
        min_prevalence=0.1,
        reference=("antibiotics_use_at_collection_analysis", "no"),
    )

    results = pd.read_csv(raw_out / "all_results.tsv", sep="\t")
    case_results = results.loc[results["metadata"].eq("case")].copy()
    case_results.to_csv(out / f"{level}_case_status_results.csv", index=False)
    if not plots:
        return results

    import matplotlib.pyplot as plt
    import numpy as np
    import seaborn as sns
    from scipy.stats import spearmanr

    significant = case_results.loc[case_results["qval"].le(0.25)]
    relative = features.div(features.sum(axis=0).replace(0, np.nan), axis=1)
    palette = {"Case": "#E97B69", "Control": "#8FA3D1"}
    plot_dir = out / f"{level}_plots"
    scatter_dir = out / f"{level}_scatter_plots"
    plot_dir.mkdir(exist_ok=True)
    scatter_dir.mkdir(exist_ok=True)

    for feature in significant["feature"].drop_duplicates():
        plot_data = metadata[["case", "age_days_npp_collection"]].copy()
        plot_data[feature] = relative.loc[feature, plot_data.index]

        ax = sns.violinplot(data=plot_data, x="case", y=feature, hue="case", palette=palette, legend=False)
        sns.stripplot(data=plot_data, x="case", y=feature, color="black", alpha=0.25, ax=ax)
        ax.figure.tight_layout()
        ax.figure.savefig(plot_dir / f"{feature}.svg")
        plt.close(ax.figure)

        lab = ds.lab["templateDNA_16SqPCR_copies"].reindex(plot_data.index)
        plot_data["16S_copies_log10"] = np.log10(lab.where(lab > 0))
        plot_data = plot_data.dropna(subset=[feature, "16S_copies_log10"])
        if plot_data.empty:
            continue
        for x_col, suffix in (("16S_copies_log10", "16S_copy"), ("age_days_npp_collection", "age")):
            if plot_data[x_col].nunique() < 2:
                continue
            grid = sns.lmplot(
                data=plot_data,
                x=x_col,
                y=feature,
                hue="case",
                order=2,
                scatter_kws={"s": 10},
                palette=palette,
            )
            corr = spearmanr(plot_data[x_col], plot_data[feature]).statistic
            grid.fig.suptitle(f"Spearman = {corr:.3f}")
            grid.fig.tight_layout()
            grid.savefig(scatter_dir / f"{feature}_{suffix}.svg")
            plt.close(grid.fig)
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--level", choices=["asv", "genus", "phylum"], default="asv")
    parser.add_argument(
        "--output-dir",
        default="results/Section_A/Q6_Masslin2_case_control",
        help="Directory for this run's output (default: %(default)s)",
    )
    parser.add_argument("--no-plots", action="store_true", help="Save model results without plots")
    args = parser.parse_args()
    run(level=args.level, output_dir=args.output_dir, plots=not args.no_plots)


if __name__ == "__main__":
    main()
