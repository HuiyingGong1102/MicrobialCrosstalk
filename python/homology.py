"""Run the supplied external GLMY executable in an isolated working directory."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import tempfile
from network_io import read_network, write_csv


def build_input(nodes, edges, edge_filter="all", filtration="signed", shift=None):
    if edge_filter not in ("all", "positive", "negative") or filtration not in ("signed", "magnitude"):
        raise ValueError("Unknown edge filter or filtration")
    selected = [e for e in edges if edge_filter == "all" or
                e["sign"] == (1 if edge_filter == "positive" else -1)]
    values = [e["signed_weight"] if filtration == "signed" else e["weight"] for e in selected]
    shift = max(0.0, -min(values, default=0.0)) if shift is None else shift
    if not math.isfinite(shift) or shift < 0 or any(w + shift < 0 for w in values):
        raise ValueError("Shift produces negative filtration values")
    mapping = {node: i + 1 for i, node in enumerate(nodes)}
    lines = [",".join(map(str, mapping.values()))]
    lines += [f'({mapping[e["source"]]},{mapping[e["target"]]},{w+shift:.17g})'
              for e, w in zip(selected, values)]
    # Preserve the original executable's dialogue: maximum dimension 4, JSON yes.
    return "\n".join(lines) + "\n#\n4\ny\n\n\n", mapping, shift


def run_glmy(nodes, edges, executable, output, edge_filter, filtration, shift, timeout):
    text, mapping, shift = build_input(nodes, edges, edge_filter, filtration, shift)
    output = Path(output); output.mkdir(parents=True, exist_ok=True)
    executable = Path(executable).resolve(strict=True)
    (output / "input.txt").write_text(text, encoding="ascii")
    metadata = dict(node_map=mapping, edge_filter=edge_filter, filtration=filtration,
                    shift=shift, vertex_birth=0,
                    executable_sha256=hashlib.sha256(executable.read_bytes()).hexdigest())
    (output / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    with tempfile.TemporaryDirectory(prefix="glmy_") as work:
        try:
            process = subprocess.run([str(executable)], input=text.encode("ascii"),
                                     capture_output=True, cwd=work, timeout=timeout, check=False)
        except subprocess.TimeoutExpired as error:
            (output / "stdout.log").write_bytes(error.stdout or b"")
            (output / "stderr.log").write_bytes(error.stderr or b"")
            raise RuntimeError("GLMY timed out; increase --timeout or reduce graph size") from error
        (output / "stdout.log").write_bytes(process.stdout)
        (output / "stderr.log").write_bytes(process.stderr)
        if process.returncode != 0:
            raise RuntimeError(f"GLMY failed with exit code {process.returncode}")
        result = Path(work) / "homology.json"
        if not result.exists():
            raise RuntimeError("GLMY did not write homology.json; inspect stdout.log")
        bars = json.loads(result.read_text(encoding="utf-8-sig"))
        (output / "homology.json").write_text(json.dumps(bars, indent=2), encoding="utf-8")
    rows = []
    for dimension, intervals in bars.items():
        for interval in intervals:
            if len(interval) != 2:
                raise ValueError("Expected [birth, death] intervals in GLMY JSON")
            birth, death = interval
            rows.append(dict(dimension=int(dimension), birth=birth, death=death,
                             infinite=(death == -1)))
    write_csv(output / "barcodes.csv", rows, ["dimension", "birth", "death", "infinite"])
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("edges", type=Path)
    parser.add_argument("--nodes", type=Path)
    parser.add_argument("--exe", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--edge-filter", choices=["all", "positive", "negative"], default="all")
    parser.add_argument("--filtration", choices=["signed", "magnitude"], default="signed")
    parser.add_argument("--shift", type=float)
    parser.add_argument("--timeout", type=float, default=300)
    args = parser.parse_args()
    nodes, edges = read_network(args.edges, args.nodes)
    run_glmy(nodes, edges, args.exe, args.output, args.edge_filter,
             args.filtration, args.shift, args.timeout)


if __name__ == "__main__":
    main()
