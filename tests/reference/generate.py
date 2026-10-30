"""Regenerate test fixtures by executing the published EMPOWER-1.0 R code."""

import argparse
import hashlib
from pathlib import Path
import shutil
import subprocess
import tempfile
import zipfile


ROOT = Path(__file__).resolve().parents[2]
ARCHIVE_SHA256 = "733c54529e357cb5d5561b9dc95d5880a00b39d4bcc7e2224f3662a8a2eb990f"


def main() -> None:
    """Verify the official archive and generate fixtures with base R."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--archive",
        type=Path,
        default=ROOT / "data/reference/empower-1.0/gmd-8-2231-2015-supplement.zip",
    )
    parser.add_argument("--output", type=Path, default=ROOT / "tests/fixtures/empower")
    args = parser.parse_args()
    if hashlib.sha256(args.archive.read_bytes()).hexdigest() != ARCHIVE_SHA256:
        raise ValueError("Supplement SHA-256 does not match the published reference")

    with tempfile.TemporaryDirectory(prefix="empower-reference-") as temporary:
        work = Path(temporary)
        with zipfile.ZipFile(args.archive) as archive:
            for name in (
                "EMPOWER1.r",
                "NPZD_parms.txt",
                "NPZD_extra.txt",
                "stations_forcing.txt",
            ):
                (work / name).write_bytes(archive.read(f"EMPOWER/{name}"))
        generated = work / "generated"
        generated.mkdir()
        subprocess.run(
            ["Rscript", str(Path(__file__).with_suffix(".R")), str(generated)],
            cwd=work,
            check=True,
        )
        args.output.mkdir(parents=True, exist_ok=True)
        for path in generated.iterdir():
            shutil.copyfile(path, args.output / path.name)
    print(f"Reference fixtures saved to {args.output}")


if __name__ == "__main__":
    main()
