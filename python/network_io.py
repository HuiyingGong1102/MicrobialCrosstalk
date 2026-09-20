"""Shared, strict edge/node CSV contract for both network scales."""
import csv
import math
from pathlib import Path


def read_network(edge_path, node_path=None):
    edge_path = Path(edge_path)
    edges, seen = [], set()
    with edge_path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        if not {"source", "target", "weight", "edge_type"} <= set(reader.fieldnames or []):
            raise ValueError("Expected source,target,weight,edge_type columns")
        for row in reader:
            u, v = row["source"], row["target"]
            w, kind = float(row["weight"]), int(row["edge_type"])
            if not u or not v or u == v or (u, v) in seen:
                raise ValueError("Missing node ID, self-loop, or duplicate directed edge")
            if not math.isfinite(w) or w <= 0 or kind not in (1, 2):
                raise ValueError("Weights must be positive and edge_type must be 1 or 2")
            signed = w if kind == 1 else -w
            if "signed_weight" in row and not math.isclose(float(row["signed_weight"]), signed):
                raise ValueError("signed_weight disagrees with edge_type")
            edges.append(dict(source=u, target=v, weight=w, sign=1 if kind == 1 else -1,
                              edge_type=kind, signed_weight=signed))
            seen.add((u, v))
    node_path = Path(node_path) if node_path else edge_path.with_name("nodes.csv")
    if not node_path.exists():
        raise FileNotFoundError("nodes.csv is required to preserve isolated nodes")
    with node_path.open(encoding="utf-8-sig", newline="") as stream:
        nodes = [row["node"] for row in csv.DictReader(stream)]
    if not nodes or len(nodes) != len(set(nodes)) or any(not n for n in nodes):
        raise ValueError("Node IDs must be nonempty and unique")
    if not {n for pair in seen for n in pair} <= set(nodes):
        raise ValueError("Edge endpoint missing from nodes.csv")
    return nodes, edges


def write_csv(path, rows, fields):
    with Path(path).open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
