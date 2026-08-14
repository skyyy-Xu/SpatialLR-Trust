#!/usr/bin/env python3
"""Spatial coordinate permutation null for stdlib_marker_lr_pilot candidates."""

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
from typing import Dict, List, Tuple

PROJECT = Path(os.environ.get('PROJECT', str(Path(__file__).resolve().parents[2])))
RUN_ID = os.environ.get('RUN_ID', 'manual')
CANDIDATES = PROJECT / 'results/task_c_pilot_baseline/stdlib_marker_lr_pilot_candidates.tsv'
MANIFEST = PROJECT / 'docs/minimal_input_manifest.tsv'
ANNOT_MANIFEST = PROJECT / 'docs/spot_annotation_manifest.tsv'
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


def adjacency_from_assigned_coords(barcodes: List[str], annotations: Dict[str, str], assigned_coords: List[Tuple[int, int, int] | None]) -> Counter:
    pos_to_indices = defaultdict(list)
    valid = []
    for i, coord in enumerate(assigned_coords):
        if coord is None:
            continue
        comp = annotations.get(barcodes[i])
        if comp not in COMPARTMENTS:
            continue
        in_tissue, row, col = coord
        if in_tissue != 1:
            continue
        pos_to_indices[(row, col)].append(i)
        valid.append(i)
    counts = Counter()
    seen = set()
    for i in valid:
        c1 = annotations[barcodes[i]]
        _, r, c = assigned_coords[i]
        for dr in [-1, 0, 1]:
            for dc in [-2, -1, 0, 1, 2]:
                if dr == 0 and dc == 0:
                    continue
                for j in pos_to_indices.get((r + dr, c + dc), []):
                    if i == j:
                        continue
                    pair = (i, j)
                    if pair in seen:
                        continue
                    seen.add(pair)
                    c2 = annotations[barcodes[j]]
                    counts[(c1, c2)] += 1
    return counts


def prepare_sample(row: Dict[str, str], annot_row: Dict[str, str]):
    barcodes = read_barcodes(PROJECT / row['barcodes_path'])
    annotations = read_annotations(PROJECT / annot_row['annotation_path'])
    coords_by_barcode = read_spatial(PROJECT / row['spatial_positions_path'])
    coords = [coords_by_barcode.get(b) for b in barcodes]
    coord_pool = [coord for coord in coords if coord is not None]
    return barcodes, annotations, coords, coord_pool


def process_sample(sample_id: str, sample_candidates: List[Dict[str, str]], manifest_by_sample, annot_by_sample, n_perm: int, seed: int):
    row = manifest_by_sample[sample_id]
    annot = annot_by_sample[sample_id]
    barcodes, annotations, coords, coord_pool = prepare_sample(row, annot)
    rng = random.Random(seed + sum(ord(c) for c in sample_id))
    null_scores = defaultdict(list)
    null_adj = defaultdict(list)
    valid_indices = [i for i, coord in enumerate(coords) if coord is not None]
    for perm in range(n_perm):
        shuffled_pool = coord_pool[:]
        rng.shuffle(shuffled_pool)
        assigned = [None] * len(coords)
        for idx, coord in zip(valid_indices, shuffled_pool):
            assigned[idx] = coord
        counts = adjacency_from_assigned_coords(barcodes, annotations, assigned)
        for cand in sample_candidates:
            key = (cand['sender_compartment'], cand['receiver_compartment'])
            adj = counts.get(key, 0)
            multiplier = float(cand['ligand_mean_sender']) * float(cand['receptor_mean_receiver'])
            null_adj[id(cand)].append(adj)
            null_scores[id(cand)].append(multiplier * math.log1p(adj))
    out_rows = []
    for cand in sample_candidates:
        observed = float(cand['pilot_score'])
        scores = null_scores[id(cand)]
        adjs = null_adj[id(cand)]
        ge = sum(1 for score in scores if score >= observed)
        empirical_p = (ge + 1) / (len(scores) + 1)
        mean_score = sum(scores) / len(scores) if scores else 0.0
        mean_adj = sum(adjs) / len(adjs) if adjs else 0.0
        ratio = observed / mean_score if mean_score > 0 else ''
        out = dict(cand)
        out.update({
            'null_model': 'spatial_coordinate_permutation',
            'null_permutations': str(n_perm),
            'null_mean_score': f'{mean_score:.8g}',
            'null_ge_observed_n': str(ge),
            'empirical_p_ge': f'{empirical_p:.8g}',
            'observed_to_null_mean_ratio': f'{ratio:.8g}' if ratio != '' else '',
            'null_mean_adjacency_pairs': f'{mean_adj:.8g}',
        })
        out_rows.append(out)
    return out_rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--permutations', type=int, default=50)
    parser.add_argument('--seed', type=int, default=20260709)
    args = parser.parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    candidates = read_tsv(CANDIDATES)
    manifest = {r['sample_id']: r for r in read_tsv(MANIFEST) if r.get('status') == 'ready'}
    annot = {r['sample_id']: r for r in read_tsv(ANNOT_MANIFEST) if r.get('status') == 'done'}
    by_sample = defaultdict(list)
    for cand in candidates:
        by_sample[cand['sample_id']].append(cand)
    all_rows = []
    for sample_id in sorted(by_sample):
        all_rows.extend(process_sample(sample_id, by_sample[sample_id], manifest, annot, args.permutations, args.seed))
    fields = list(all_rows[0].keys())
    out_path = OUT_DIR / 'stdlib_marker_lr_spatial_null.tsv'
    write_tsv(out_path, all_rows, fields)
    pvals = [float(r['empirical_p_ge']) for r in all_rows]
    summary = {
        'run_id': RUN_ID,
        'null_model': 'spatial_coordinate_permutation',
        'source_method': 'stdlib_marker_lr_pilot',
        'candidate_rows': len(all_rows),
        'samples': len(by_sample),
        'permutations': args.permutations,
        'seed': args.seed,
        'p_le_0_05': sum(1 for p in pvals if p <= 0.05),
        'p_le_0_10': sum(1 for p in pvals if p <= 0.10),
        'output': rel(out_path),
    }
    (OUT_DIR / 'stdlib_marker_lr_spatial_null_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    (RUN_DIR / 'spatial_null_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
