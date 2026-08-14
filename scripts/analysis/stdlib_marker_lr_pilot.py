#!/usr/bin/env python3
"""Stdlib-only pilot ligand-receptor baseline for SpatialLR-Trust.

This creates standardized LR candidate tables from coarse spot compartments,
selected LR pairs, raw mtx counts, and a simple array-coordinate adjacency proxy.
It is explicitly not CellChat/COMMOT/CytoSignal.
"""

from __future__ import annotations

import csv
import gzip
import json
import math
import os
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple

PROJECT = Path(os.environ.get('PROJECT', str(Path(__file__).resolve().parents[2])))
RUN_ID = os.environ.get('RUN_ID', 'manual')
MANIFEST = PROJECT / 'docs/minimal_input_manifest.tsv'
ANNOT_MANIFEST = PROJECT / 'docs/spot_annotation_manifest.tsv'
LR_PANEL = PROJECT / 'docs/pilot_lr_panel.tsv'
OUT_DIR = PROJECT / 'results/task_c_pilot_baseline'
RUN_DIR = PROJECT / 'runs' / RUN_ID
COMPARTMENTS = ['tumor_like', 'immune', 'stroma_fibroblast', 'endothelial']
METHOD = 'stdlib_marker_lr_pilot'


def rel(path: Path) -> str:
    return str(path.relative_to(PROJECT))


def open_text(path: Path):
    if path.suffix == '.gz':
        return gzip.open(path, 'rt', errors='replace')
    return path.open('rt', errors='replace')


def read_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def read_lr_panel() -> List[Dict[str, str]]:
    rows = read_tsv(LR_PANEL)
    return [{**row, 'ligand': row['ligand'].upper(), 'receptor': row['receptor'].upper()} for row in rows]


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
            hit = next((name for name in names if name in wanted), None)
            if hit:
                index_to_gene[idx] = hit
    return index_to_gene


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


def read_mtx_selected(path: Path, index_to_gene: Dict[int, str], barcode_to_comp: Dict[int, str]):
    gene_comp_sum = defaultdict(lambda: Counter())
    gene_comp_spots = defaultdict(lambda: defaultdict(set))
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
            if gene_idx not in index_to_gene:
                continue
            spot_idx = int(fields[1]) - 1
            comp = barcode_to_comp.get(spot_idx)
            if comp not in COMPARTMENTS:
                continue
            value = float(fields[2])
            gene = index_to_gene[gene_idx]
            gene_comp_sum[gene][comp] += value
            gene_comp_spots[gene][comp].add(spot_idx)
    return shape, gene_comp_sum, gene_comp_spots


def adjacency_counts(barcodes: List[str], annotations: Dict[str, str], coords: Dict[str, Tuple[int, int, int]]):
    # Uses array row/column local neighborhoods. Counts ordered sender->receiver compartment pairs.
    pos_to_indices = defaultdict(list)
    valid_indices = []
    for i, barcode in enumerate(barcodes):
        comp = annotations.get(barcode)
        coord = coords.get(barcode)
        if comp not in COMPARTMENTS or coord is None:
            continue
        in_tissue, row, col = coord
        if in_tissue != 1:
            continue
        pos_to_indices[(row, col)].append(i)
        valid_indices.append(i)
    counts = Counter()
    seen_pairs = set()
    for i in valid_indices:
        b1 = barcodes[i]
        c1 = annotations[b1]
        _, r, c = coords[b1]
        for dr in [-1, 0, 1]:
            for dc in [-2, -1, 0, 1, 2]:
                if dr == 0 and dc == 0:
                    continue
                for j in pos_to_indices.get((r + dr, c + dc), []):
                    if i == j:
                        continue
                    pair = (i, j)
                    if pair in seen_pairs:
                        continue
                    seen_pairs.add(pair)
                    c2 = annotations[barcodes[j]]
                    counts[(c1, c2)] += 1
    return counts


def process_sample(row: Dict[str, str], annot_row: Dict[str, str], lr_rows: List[Dict[str, str]]):
    wanted = {lr['ligand'] for lr in lr_rows} | {lr['receptor'] for lr in lr_rows}
    matrix = PROJECT / row['expression_path']
    features = PROJECT / row['features_path']
    barcodes_path = PROJECT / row['barcodes_path']
    spatial_path = PROJECT / row['spatial_positions_path']
    annot_path = PROJECT / annot_row['annotation_path']
    barcodes = read_barcodes(barcodes_path)
    annotations = read_annotations(annot_path)
    coords = read_spatial(spatial_path)
    comp_counts = Counter(annotations.get(b, 'missing') for b in barcodes)
    barcode_to_comp = {i: annotations.get(b) for i, b in enumerate(barcodes)}
    index_to_gene = read_features(features, wanted)
    shape, gene_comp_sum, gene_comp_spots = read_mtx_selected(matrix, index_to_gene, barcode_to_comp)
    adj = adjacency_counts(barcodes, annotations, coords)
    candidates = []
    for lr in lr_rows:
        ligand = lr['ligand']
        receptor = lr['receptor']
        for sender in COMPARTMENTS:
            sender_n = comp_counts.get(sender, 0)
            if sender_n <= 0:
                continue
            ligand_sum = gene_comp_sum[ligand][sender]
            ligand_mean = ligand_sum / sender_n if sender_n else 0.0
            ligand_pct = len(gene_comp_spots[ligand][sender]) / sender_n if sender_n else 0.0
            for receiver in COMPARTMENTS:
                receiver_n = comp_counts.get(receiver, 0)
                if receiver_n <= 0:
                    continue
                receptor_sum = gene_comp_sum[receptor][receiver]
                receptor_mean = receptor_sum / receiver_n if receiver_n else 0.0
                receptor_pct = len(gene_comp_spots[receptor][receiver]) / receiver_n if receiver_n else 0.0
                adjacency = adj.get((sender, receiver), 0)
                score = ligand_mean * receptor_mean * math.log1p(adjacency)
                if ligand_mean <= 0 or receptor_mean <= 0 or adjacency <= 0:
                    continue
                candidates.append({
                    'dataset': row['dataset'],
                    'sample_id': row['sample_id'],
                    'cancer': row['cancer'],
                    'method': METHOD,
                    'sender_compartment': sender,
                    'receiver_compartment': receiver,
                    'ligand': ligand,
                    'receptor': receptor,
                    'pathway': lr['pathway'],
                    'ligand_mean_sender': f'{ligand_mean:.8g}',
                    'receptor_mean_receiver': f'{receptor_mean:.8g}',
                    'ligand_pct_sender': f'{ligand_pct:.8g}',
                    'receptor_pct_receiver': f'{receptor_pct:.8g}',
                    'adjacency_pairs': str(adjacency),
                    'pilot_score': f'{score:.8g}',
                    'annotation_method': annot_row['annotation_method'],
                    'notes': 'Pilot stdlib candidate; not CellChat/COMMOT/CytoSignal.',
                })
    sample_summary = {
        'dataset': row['dataset'],
        'sample_id': row['sample_id'],
        'cancer': row['cancer'],
        'matrix_rows': str(shape[0] if shape else ''),
        'matrix_cols': str(shape[1] if shape else ''),
        'matrix_nnz': str(shape[2] if shape else ''),
        'candidate_n': str(len(candidates)),
        'genes_detected_in_panel': str(len(set(index_to_gene.values()))),
        'tumor_like_spots': str(comp_counts.get('tumor_like', 0)),
        'immune_spots': str(comp_counts.get('immune', 0)),
        'stroma_fibroblast_spots': str(comp_counts.get('stroma_fibroblast', 0)),
        'endothelial_spots': str(comp_counts.get('endothelial', 0)),
        'ambiguous_or_low_signal_spots': str(comp_counts.get('ambiguous_or_low_signal', 0)),
    }
    return candidates, sample_summary


def write_tsv(path: Path, rows: List[Dict[str, str]], fields: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    manifest = [row for row in read_tsv(MANIFEST) if row.get('status') == 'ready']
    annot = {row['sample_id']: row for row in read_tsv(ANNOT_MANIFEST) if row.get('status') == 'done'}
    lr_rows = read_lr_panel()
    all_candidates = []
    summaries = []
    for row in manifest:
        if row['sample_id'] not in annot:
            raise RuntimeError(f"No annotation row for {row['sample_id']}")
        candidates, summary = process_sample(row, annot[row['sample_id']], lr_rows)
        all_candidates.extend(candidates)
        summaries.append(summary)
    cand_fields = ['dataset','sample_id','cancer','method','sender_compartment','receiver_compartment','ligand','receptor','pathway','ligand_mean_sender','receptor_mean_receiver','ligand_pct_sender','receptor_pct_receiver','adjacency_pairs','pilot_score','annotation_method','notes']
    summary_fields = ['dataset','sample_id','cancer','matrix_rows','matrix_cols','matrix_nnz','candidate_n','genes_detected_in_panel','tumor_like_spots','immune_spots','stroma_fibroblast_spots','endothelial_spots','ambiguous_or_low_signal_spots']
    cand_path = OUT_DIR / 'stdlib_marker_lr_pilot_candidates.tsv'
    sample_path = OUT_DIR / 'stdlib_marker_lr_pilot_sample_summary.tsv'
    write_tsv(cand_path, all_candidates, cand_fields)
    write_tsv(sample_path, summaries, summary_fields)
    summary = {
        'run_id': RUN_ID,
        'method': METHOD,
        'samples': len(summaries),
        'lr_pairs': len(lr_rows),
        'candidates': len(all_candidates),
        'candidate_table': rel(cand_path),
        'sample_summary': rel(sample_path),
    }
    (OUT_DIR / 'stdlib_marker_lr_pilot_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    (RUN_DIR / 'pilot_lr_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
