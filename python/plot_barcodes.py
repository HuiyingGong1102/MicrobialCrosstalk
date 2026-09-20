"""Plot GLMY barcodes on the recorded filtration scale; -1 means infinite death."""
import argparse
import csv
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("barcodes", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    with args.barcodes.open(encoding="utf-8", newline="") as stream:
        rows = [{"dimension": int(r["dimension"]), "birth": float(r["birth"]),
                 "death": float(r["death"])} for r in csv.DictReader(stream)]
    finite = [r[key] for r in rows for key in ("birth", "death") if r[key] != -1]
    left = min([0, *finite]); right = max([1, *finite]); right += (right-left)*0.1
    fig, axes = plt.subplots(3, 1, figsize=(8, 7), constrained_layout=True)
    for dim, ax in enumerate(axes):
        selected = sorted([r for r in rows if r["dimension"] == dim], key=lambda r: r["birth"])
        for i, r in enumerate(selected):
            end = right if r["death"] == -1 else r["death"]
            ax.hlines(i, r["birth"], end)
            if r["death"] == -1:
                ax.plot(end, i, marker=">", color="black")
        ax.set(title=f"H{dim}", xlim=(left, right), ylabel="Interval")
    axes[-1].set_xlabel("Filtration threshold (including recorded shift)")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output)
    plt.close(fig)


if __name__ == "__main__":
    main()
