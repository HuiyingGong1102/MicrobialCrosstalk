"""Directed signed-control Forman and LLY on the undirected unweighted projection."""
import argparse
import json
from pathlib import Path
import numpy as np
import networkx as nx
from scipy.optimize import linprog
from network_io import read_network, write_csv


def curvature(nodes, edges):
    graph = nx.Graph()
    graph.add_nodes_from(nodes)
    graph.add_edges_from((e["source"], e["target"]) for e in edges)
    distances = dict(nx.all_pairs_shortest_path_length(graph))
    pair_rows, lly_by_pair = [], {}
    for u, v in graph.edges():
        support = list(dict.fromkeys([u, v, *graph[u], *graph[v]]))
        n = len(support)
        index = {node: i for i, node in enumerate(support)}
        a, b = np.zeros(n), np.zeros(n)
        a[index[u]], b[index[v]] = 0.5, 0.5
        for z in graph[u]:
            a[index[z]] += 0.5 / graph.degree(u)
        for z in graph[v]:
            b[index[z]] += 0.5 / graph.degree(v)
        cost = np.array([[distances[x][y] for y in support] for x in support])
        constraints = np.vstack([np.kron(np.eye(n), np.ones((1, n))),
                                 np.kron(np.ones((1, n)), np.eye(n))])
        solution = linprog(cost.ravel(), A_eq=constraints, b_eq=np.r_[a, b],
                           bounds=(0, None), method="highs")
        if not solution.success:
            raise RuntimeError(solution.message)
        lly = 2 * (1 - solution.fun)
        pair = tuple(sorted((u, v)))
        lly_by_pair[pair] = lly
        pair_rows.append(dict(pair_u=pair[0], pair_v=pair[1], LLY=lly,
                              OR_alpha_0_5=1-solution.fun, W1_alpha_0_5=solution.fun))
    edge_rows = []
    for e in edges:
        u, v = e["source"], e["target"]
        adjacent = [a for a in edges if (a["target"] == u and a["source"] != v)
                    or (a["source"] == v and a["target"] != u)]
        unsigned = sum(np.sqrt(e["weight"] / a["weight"]) for a in adjacent)
        signed = sum(e["sign"] * a["sign"] * np.sqrt(e["weight"] / a["weight"]) for a in adjacent)
        edge_rows.append(dict(e, Directed_Forman=2-unsigned, SC_Forman=2-signed,
                              LLY=lly_by_pair[tuple(sorted((u, v)))]))
    node_rows = []
    for node in nodes:
        incident = [e["SC_Forman"] for e in edge_rows if node in (e["source"], e["target"])]
        pair_values = [v for pair, v in lly_by_pair.items() if node in pair]
        node_rows.append(dict(node=node, in_degree=sum(e["target"] == node for e in edges),
                              out_degree=sum(e["source"] == node for e in edges),
                              LLY_degree=graph.degree(node),
                              SC_Forman_mean=float(np.mean(incident)) if incident else None,
                              LLY_mean=float(np.mean(pair_values)) if pair_values else None))
    return edge_rows, node_rows, pair_rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("edges", type=Path)
    parser.add_argument("--nodes", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    nodes, edges = read_network(args.edges, args.nodes)
    edge_rows, node_rows, pairs = curvature(nodes, edges)
    args.output.mkdir(parents=True, exist_ok=True)
    write_csv(args.output / "edge_curvature.csv", edge_rows,
              ["source", "target", "weight", "edge_type", "sign", "Directed_Forman", "SC_Forman", "LLY"])
    write_csv(args.output / "node_curvature.csv", node_rows,
              ["node", "in_degree", "out_degree", "LLY_degree", "SC_Forman_mean", "LLY_mean"])
    write_csv(args.output / "pair_curvature.csv", pairs,
              ["pair_u", "pair_v", "LLY", "OR_alpha_0_5", "W1_alpha_0_5"])
    summary = dict(nodes=len(nodes), directed_edges=len(edges), undirected_pairs=len(pairs),
                   mean_SC_Forman=float(np.mean([e["SC_Forman"] for e in edge_rows])) if edges else None,
                   mean_LLY=float(np.mean([p["LLY"] for p in pairs])) if pairs else None,
                   LLY_projection="undirected, unweighted, unsigned", alpha=0.5)
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
