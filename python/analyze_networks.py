"""Analyze every exported coarse/fine network, preserving its directory hierarchy."""
import argparse
from pathlib import Path
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("network_root", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--kind", choices=["curvature", "homology", "all"], default="all")
    parser.add_argument("--exe", type=Path)
    parser.add_argument("--timeout", type=float, default=300)
    parser.add_argument("--shift", type=float)
    parser.add_argument("--filtration", choices=["signed", "magnitude"], default="signed")
    args = parser.parse_args()
    if args.kind != "curvature" and args.exe is None:
        parser.error("--exe is required for homology")
    files = sorted(args.network_root.rglob("edges.csv"))
    if not files:
        parser.error("No edges.csv files found under network_root")
    scripts = Path(__file__).resolve().parent
    for edges in files:
        relative = edges.parent.relative_to(args.network_root)
        print(f"Analyzing {relative}", flush=True)
        if args.kind in ("curvature", "all"):
            subprocess.run([sys.executable, str(scripts/"curvature.py"), str(edges),
                            "--output", str(args.output/"curvature"/relative)], check=True)
        if args.kind in ("homology", "all"):
            for sign in ("all", "positive", "negative"):
                output = args.output/"homology"/relative/sign
                command = [sys.executable, str(scripts/"homology.py"), str(edges),
                           "--exe", str(args.exe), "--output", str(output),
                           "--edge-filter", sign, "--filtration", args.filtration,
                           "--timeout", str(args.timeout)]
                if args.shift is not None:
                    command += ["--shift", str(args.shift)]
                subprocess.run(command, check=True)
                subprocess.run([sys.executable, str(scripts/"plot_barcodes.py"),
                                str(output/"barcodes.csv"), "--output", str(output/"barcodes.pdf")], check=True)


if __name__ == "__main__":
    main()
