#!/usr/bin/env python3
"""Spot label permutation null for stdlib_marker_lr_pilot candidates."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
import random
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Set, Tuple

PROJECT = Path(os.environ.get('PROJECT', str(Path(__file__).resolve().parents[2])))
RUN_ID = os.environ.get('RUN_ID', 'manual')
CANDIDATES = PROJECT / 'results/task_c_pilot_baseline/stdlib_marker_lr_pilot_candidates.tsv'
MANIFEST = PROJECT / 'docs/minimal_input_manifest.tsv'
ANNOT_MANIFEST = PROJECT / 'docs/spot_annotation_manifest.tsv'
LR_PANEL = PROJECT / 'docs/pilot_lr_panel.tsv'
OUT_DIR = PROJECT / 'results/task_d_null_models'
RUN_DIR = PROJECT / 'runs' / RUN_ID
COMPARTMENTS = ['tumor_like', 'immune', 'stroma_fibroblast', 'endothelial']


def rel(path: Path) -> str:
    return str(path.relative_to(PROJECT))


def open_text(path: Path):
    if path.suffix == '.gz':
        return gzip.open(path, 'rt', errors='replace')
    return path.open('rt', errors='replace')


def read_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write_tsv(path: Path, rows: List[Dict[str, str]], fields: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def read_barcodes(path: Path) -> List[str]:
    with open_text(path) as handle:
        return [line.rstrip('\n').split('\t', 1)[0] for line in handle if line.strip()]


def read_annotations(path: Path) -> Dict[str, str]:
    out = {}
    with gzip.open(path, 'rt', errors='replace', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        for row in reader:
            out[row['barcode']] = row['spot_compartment']
    return out


def read_spatial(path: Path) -> Dict[str, Tuple[int, int, int]]:
    coords = {}
    with open_text(path) as handle:
        for idx, line in enumerate(handle):
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split(',')
            if idx == 0 and any(token in parts for token in ['barcode','in_tissue','array_row','array_col']):
                continue
            if len(parts) < 4:
                continue
            try:
                coords[parts[0]] = (int(float(parts[1])), int(float(parts[2])), int(float(parts[3])))
            except ValueError:
                continue
    return coords


def read_lr_genes() -> Set[str]:
    genes = set()
    for row in read_tsv(LR_PANEL):
        genes.add(row['ligand'].upper())
        genes.add(row['receptor'].upper())
    return genes


def read_features(path: Path, wanted: Set[str]) -> Dict[int, str]:
    index_to_gene = {}
    with open_text(path) as handle:
        for idx, line in enumerate(handle, start=1):
            fields = line.rstrip('\n').split('\t')
            names = []
            if fields:
                names.append(fields[0].upper())
            if len(fields) > 1:
                names.append(fields[1].upper())
            for name in names:
                if name in wanted:
                    index_to_gene[idx] = name
                    break
    return index_to_gene


def read_selected_entries(path: Path, index_to_gene: Dict[int, str]):
    entries = []
    shape = None
    with open_text(path) as handle:
        for line in handle:
            if line.startswith('%'):
                continue
            fields = line.strip().split()
            if len(fields) < 3:
                continue
            if shape is None:
                shape = (int(fields[0]), int(fields[1]), int(fields[2]))
                continue
            gene_idx = int(fields[0])
            gene = index_to_gene.get(gene_idx)
            if gene is None:
                continue
            entries.append((gene, int(fields[1]) - 1, float(fields[2])))
    return shape, entries


def adjacency_counts(barcodes: List[str], labels: List[str], coords: List[Tuple[int, int, int] | None]) -> Counter:
    pos_to_indices = defaultdict(list)
    valid = []
    for i, coord in enumerate(coords):
        if coord is None or labels[i] not in COMPARTMENTS:
            continue
        in_tissue, row, col = coord
        if in_tissue != 1:
            continue
        pos_to_indices[(row, col)].append(i)
        valid.append(i)
    counts = Counter()
    seen = set()
    for i in valid:
        c1 = labels[i]
        _, row, col = coords[i]
        for dr in [-1, 0, 1]:
            for dc in [-2, -1, 0, 1, 2]:
                if dr == 0 and dc == 0:
                    continue
                for j in pos_to_indices.get((row + dr, col + dc), []):
                    if i == j:
                        continue
                    pair = (i, j)
                    if pair in seen:
                        continue
                    seen.add(pair)
                    c2 = labels[j]
                    counts[(c1, c2)] += 1
    return counts


def expression_by_comp(entries, labels: List[str], comp_counts: Counter):
    gene_comp_sum = defaultdict(lambda: Counter())
    gene_comp_spots = defaultdict(lambda: defaultdict(set))
    for gene, spot_idx, value in entries:
        comp = labels[spot_idx]
        if comp not in COMPARTMENTS:
            continue
        gene_comp_sum[gene][comp] += value
        gene_comp_spots[gene][comp].add(spot_idx)
    stats = {}
    genes = set(gene_comp_sum) | set(gene_comp_spots)
    for gene in genes:
        stats[gene] = {}
        for comp in COMPARTMENTS:
            n = comp_counts.get(comp, 0)
            mean = gene_comp_sum[gene][comp] / n if n else 0.0
            pct = len(gene_comp_spots[gene][comp]) / n if n else 0.0
            stats[gene][comp] = (mean, pct)
    return stats


def prepare_sample(row, annot_row, wanted_genes):
    barcodes = read_barcodes(PROJECT / row['barcodes_path'])
    annotations = read_annotations(PROJECT / annot_row['annotation_path'])
    labels = [annotations.get(b, 'missing') for b in barcodes]
    coords_by_barcode = read_spatial(PROJECT / row['spatial_positions_path'])
    coords = [coords_by_barcode.get(b) for b in barcodes]
    index_to_gene = read_features(PROJECT / row['features_path'], wanted_genes)
    shape, entries = read_selected_entries(PROJECT / row['expression_path'], index_to_gene)
    return barcodes, labels, coords, entries


def process_sample(sample_id, sample_candidates, manifest_by_sample, annot_by_sample, wanted_genes, n_perm, seed):
    row = manifest_by_sample[sample_id]
    annot = annot_by_sample[sample_id]
    barcodes, labels, coords, entries = prepare_sample(row, annot, wanted_genes)
    rng = random.Random(seed + sum(ord(c) for c in sample_id))
    shuffled_positions = list(range(len(labels)))
    null_scores = defaultdict(list)
    null_adj = defaultdict(list)
    null_ligand_mean = defaultdict(list)
    null_receptor_mean = defaultdict(list)
    for _ in range(n_perm):
        rng.shuffle(shuffled_positions)
        perm_labels = [labels[i] for i in shuffled_positions]
        comp_counts = Counter(x for x in perm_labels if x in COMPARTMENTS)
        expr = expression_by_comp(entries, perm_labels, comp_counts)
        adj = adjacency_counts(barcodes, perm_labels, coords)
        for cand in sample_candidates:
            sender = cand['sender_compartment']
            receiver = cand['receiver_compartment']
            ligand = cand['ligand']
            receptor = cand['receptor']
            ligand_mean = expr.get(ligand, {}).get(sender, (0.0, 0.0))[0]
            receptor_mean = expr.get(receptor, {}).get(receiver, (0.0, 0.0))[0]
            adjacency = adj.get((sender, receiver), 0)
            score = ligand_mean * receptor_mean * math.log1p(adjacency)
            key = id(cand)
            null_scores[key].append(score)
            null_adj[key].append(adjacency)
            null_ligand_mean[key].append(ligand_mean)
            null_receptor_mean[key].append(receptor_mean)
    out_rows = []
    for cand in sample_candidates:
        key = id(cand)
        observed = float(cand['pilot_score'])
        scores = null_scores[key]
        ge = sum(1 for score in scores if score >= observed)
        empirical_p = (ge + 1) / (len(scores) + 1)
        mean_score = sum(scores) / len(scores) if scores else 0.0
        ratio = observed / mean_score if mean_score > 0 else ''
        out = dict(cand)
        out.update({
            'null_model': 'spot_compartment_label_permutation',
            'null_permutations': str(n_perm),
            'null_mean_score': f'{mean_score:.8g}',
            'null_ge_observed_n': str(ge),
            'empirical_p_ge': f'{empirical_p:.8g}',
            'observed_to_null_mean_ratio': f'{ratio:.8g}' if ratio != '' else '',
            'null_mean_adjacency_pairs': f'{(sum(null_adj[key]) / len(null_adj[key])) if null_adj[key] else 0.0:.8g}',
            'null_mean_ligand_mean_sender': f'{(sum(null_ligand_mean[key]) / len(null_ligand_mean[key])) if null_ligand_mean[key] else 0.0:.8g}',
            'null_mean_receptor_mean_receiver': f'{(sum(null_receptor_mean[key]) / len(null_receptor_mean[key])) if null_receptor_mean[key] else 0.0:.8g}',
        })
        out_rows.append(out)
    return out_rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--permutations', type=int, default=30)
    parser.add_argument('--seed', type=int, default=20260709)
    args = parser.parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    candidates = read_tsv(CANDIDATES)
    manifest = {r['sample_id']: r for r in read_tsv(MANIFEST) if r.get('status') == 'ready'}
    annot = {r['sample_id']: r for r in read_tsv(ANNOT_MANIFEST) if r.get('status') == 'done'}
    wanted_genes = read_lr_genes()
    by_sample = defaultdict(list)
    for cand in candidates:
        by_sample[cand['sample_id']].append(cand)
    all_rows = []
    for sample_id in sorted(by_sample):
        all_rows.extend(process_sample(sample_id, by_sample[sample_id], manifest, annot, wanted_genes, args.permutations, args.seed))
    fields = list(all_rows[0].keys())
    out_path = OUT_DIR / 'stdlib_marker_lr_label_null.tsv'
    write_tsv(out_path, all_rows, fields)
    pvals = [float(r['empirical_p_ge']) for r in all_rows]
    summary = {
        'run_id': RUN_ID,
        'null_model': 'spot_compartment_label_permutation',
        'source_method': 'stdlib_marker_lr_pilot',
        'candidate_rows': len(all_rows),
        'samples': len(by_sample),
        'permutations': args.permutations,
        'seed': args.seed,
        'p_le_0_05': sum(1 for p in pvals if p <= 0.05),
        'p_le_0_10': sum(1 for p in pvals if p <= 0.10),
        'output': rel(out_path),
    }
    (OUT_DIR / 'stdlib_marker_lr_label_null_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    (RUN_DIR / 'label_null_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
