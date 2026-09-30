"""Invoke the R MaAsLin2 package without requiring the rpy2 bridge."""

import os
from pathlib import Path
import shutil
import subprocess
import tempfile


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
R_RUNNER = REPOSITORY_ROOT / "stats_ana" / "run_maaslin.R"


def run_maaslin(features, metadata, output, fixed_effects, min_prevalence=0.1,
                reference=None):
    """Write aligned input tables and call MaAsLin2 through Rscript."""
    rscript = os.environ.get("RSCRIPT") or shutil.which("Rscript")
    if not rscript:
        raise RuntimeError("Rscript is not on PATH; install R and make Rscript available.")

    output = Path(output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        raise FileExistsError(f"MaAsLin2 output already exists: {output}")
    with tempfile.TemporaryDirectory(prefix="maaslin2-input-") as temp_dir:
        temp_dir = Path(temp_dir)
        feature_file = temp_dir / "features.tsv"
        metadata_file = temp_dir / "metadata.tsv"
        features.to_csv(feature_file, sep="\t", index=True, index_label="feature")
        metadata.to_csv(metadata_file, sep="\t", index=True, index_label="sample")
        command = [
            rscript,
            str(R_RUNNER),
            str(feature_file),
            str(metadata_file),
            str(output),
            ",".join(fixed_effects),
            str(min_prevalence),
            reference[0] if reference else "",
            reference[1] if reference else "",
        ]
        try:
            subprocess.run(command, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as exc:
            details = (exc.stderr or exc.stdout or "").strip()
            raise RuntimeError(
                "MaAsLin2 failed. Check the R error below and confirm that the "
                f"Maaslin2 package is installed.\n{details[-6000:]}"
            ) from exc
    print(f"MaAsLin2 completed: {output}")
