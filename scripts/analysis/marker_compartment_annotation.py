#!/usr/bin/env python3
"""Marker-based spot compartment annotation for mtx Visium inputs.

This is a conservative stdlib-only bridge until a Scanpy/Seurat/deconvolution
environment is available. It emits coarse spot compartments, not cell types.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Set

PROJECT = Path(__file__).resolve().parents[2]
MARKERS = {
    'tumor_like': {
        'EPCAM','KRT8','KRT18','KRT19','KRT7','KRT5','KRT14','KRT20','CDH1','MUC1','CLDN4','CEACAM5','CEACAM6','VIL1','CDX2',
        'MLANA','PMEL','TYR','MITF','SOX10','S100B','MART1'
    },
    'immune': {'PTPRC','CD3D','CD3E','CD3G','CD2','TRAC','MS4A1','CD79A','CD79B','NKG7','GNLY','LST1','LYZ','C1QA','C1QB','FCGR3A'},
    'stroma_fibroblast': {'COL1A1','COL1A2','COL3A1','COL6A1','COL6A2','DCN','LUM','FAP','ACTA2','TAGLN','PDGFRA','PDGFRB','THY1'},
    'endothelial': {'PECAM1','VWF','KDR','FLT1','ENG','PLVAP','CLDN5','RAMP2','ESAM'},
}
LABELS = list(MARKERS)
AMBIG = 'ambiguous_or_low_signal'


def open_text(path: Path):
    if path.suffix == '.gz':
        return gzip.open(path, 'rt', errors='replace')
    return path.open('rt', errors='replace')


def rel(path: Path) -> str:
    return str(path.relative_to(PROJECT))


def read_manifest(path: Path) -> List[Dict[str, str]]:
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def select_rows(rows: List[Dict[str, str]], mode: str) -> List[Dict[str, str]]:
    ready = [r for r in rows if r.get('status') == 'ready']
    if mode == 'all':
        return ready
    wanted = {
        'GSE300445': 'GSM9060732_AdjIII-0019',
        'GSE283052': 'GSM8655157_P02',
        'GSE292299': 'GSM8855716_NSCLC_P14',
    }
    return [r for r in ready if wanted.get(r['dataset']) == r['sample_id']]


def read_features(path: Path):
    gene_to_indices: Dict[str, Set[int]] = defaultdict(set)
    index_to_labels: Dict[int, Set[str]] = defaultdict(set)
    marker_hits = defaultdict(set)
    marker_to_labels = defaultdict(set)
    for label, genes in MARKERS.items():
        for gene in genes:
            marker_to_labels[gene.upper()].add(label)
    with open_text(path) as handle:
        for idx, line in enumerate(handle, start=1):
            fields = line.rstrip('\n').split('\t')
            candidates = []
            if fields:
                candidates.append(fields[0])
            if len(fields) > 1:
                candidates.append(fields[1])
            for candidate in candidates:
                gene = candidate.upper()
                if gene in marker_to_labels:
                    gene_to_indices[gene].add(idx)
                    for label in marker_to_labels[gene]:
                        index_to_labels[idx].add(label)
                        marker_hits[label].add(gene)
    return index_to_labels, marker_hits


def read_barcodes(path: Path) -> List[str]:
    with open_text(path) as handle:
        return [line.rstrip('\n').split('\t', 1)[0] for line in handle if line.strip()]


def init_scores(n: int):
    return {label: [0.0] * n for label in LABELS}


def parse_matrix_scores(path: Path, index_to_labels: Dict[int, Set[str]], spot_n: int):
    scores = init_scores(spot_n)
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
            if gene_idx not in index_to_labels:
                continue
            spot_idx = int(fields[1]) - 1
            value = float(fields[2])
            for label in index_to_labels[gene_idx]:
                scores[label][spot_idx] += value
    return shape, scores


def label_spots(scores: Dict[str, List[float]], barcodes: List[str], ratio: float):
    rows = []
    counts = Counter()
    for i, barcode in enumerate(barcodes):
        vals = [(label, scores[label][i]) for label in LABELS]
        vals.sort(key=lambda x: x[1], reverse=True)
        top_label, top_score = vals[0]
        second_score = vals[1][1] if len(vals) > 1 else 0.0
        if top_score <= 0:
            assigned = AMBIG
            confidence = 'low'
        elif second_score <= 0 or top_score >= ratio * second_score:
            assigned = top_label
            confidence = 'medium' if second_score == 0 else 'medium'
        else:
            assigned = AMBIG
            confidence = 'low'
        counts[assigned] += 1
        out = {
            'barcode': barcode,
            'spot_compartment': assigned,
            'annotation_confidence': confidence,
            'top_score': f'{top_score:.6g}',
            'second_score': f'{second_score:.6g}',
        }
        for label in LABELS:
            out[f'{label}_marker_score'] = f'{scores[label][i]:.6g}'
        rows.append(out)
    return rows, counts


def write_annotation(path: Path, rows: List[Dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = ['barcode','spot_compartment','annotation_confidence','top_score','second_score'] + [f'{label}_marker_score' for label in LABELS]
    with gzip.open(path, 'wt', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def annotate_one(row: Dict[str, str], out_root: Path, ratio: float) -> Dict[str, str]:
    matrix = PROJECT / row['expression_path']
    features = PROJECT / row['features_path']
    barcodes_path = PROJECT / row['barcodes_path']
    index_to_labels, marker_hits = read_features(features)
    barcodes = read_barcodes(barcodes_path)
    shape, scores = parse_matrix_scores(matrix, index_to_labels, len(barcodes))
    annotations, counts = label_spots(scores, barcodes, ratio)
    out_path = out_root / row['dataset'] / f"{row['sample_id']}_spot_compartments.tsv.gz"
    write_annotation(out_path, annotations)
    total_hits = sum(len(v) for v in marker_hits.values())
    return {
        'dataset': row['dataset'],
        'sample_id': row['sample_id'],
        'cancer': row['cancer'],
        'spots': str(len(barcodes)),
        'matrix_rows': str(shape[0] if shape else ''),
        'matrix_cols': str(shape[1] if shape else ''),
        'matrix_nnz': str(shape[2] if shape else ''),
        'marker_genes_detected_total': str(total_hits),
        'tumor_like_markers_detected': str(len(marker_hits['tumor_like'])),
        'immune_markers_detected': str(len(marker_hits['immune'])),
        'stroma_fibroblast_markers_detected': str(len(marker_hits['stroma_fibroblast'])),
        'endothelial_markers_detected': str(len(marker_hits['endothelial'])),
        'tumor_like_spots': str(counts['tumor_like']),
        'immune_spots': str(counts['immune']),
        'stroma_fibroblast_spots': str(counts['stroma_fibroblast']),
        'endothelial_spots': str(counts['endothelial']),
        'ambiguous_or_low_signal_spots': str(counts[AMBIG]),
        'annotation_path': rel(out_path),
        'status': 'done',
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', choices=['pilot', 'all'], default='pilot')
    parser.add_argument('--run-id', required=True)
    parser.add_argument('--ratio', type=float, default=1.5)
    args = parser.parse_args()
    manifest = read_manifest(PROJECT / 'docs/minimal_input_manifest.tsv')
    rows = select_rows(manifest, args.mode)
    out_root = PROJECT / 'data/processed/spot_annotations' / args.run_id
    result_root = PROJECT / 'results/task_b_annotation'
    result_root.mkdir(parents=True, exist_ok=True)
    summaries = [annotate_one(row, out_root, args.ratio) for row in rows]
    summary_path = result_root / f'marker_compartment_{args.mode}_summary.tsv'
    fields = list(summaries[0].keys())
    with summary_path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(summaries)
    json_path = result_root / f'marker_compartment_{args.mode}_summary.json'
    json_path.write_text(json.dumps({'run_id': args.run_id, 'mode': args.mode, 'samples': len(summaries), 'summary_tsv': rel(summary_path)}, indent=2) + '\n')
    print(json.dumps({'mode': args.mode, 'samples': len(summaries), 'summary_tsv': rel(summary_path), 'annotation_root': rel(out_root)}, indent=2))


if __name__ == '__main__':
    main()
