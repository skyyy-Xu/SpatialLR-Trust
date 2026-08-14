#!/usr/bin/env python3
"""Build a pilot receiver-compartment pathway target activation proxy."""
from __future__ import annotations
import argparse
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
TARGET_PANEL = PROJECT / 'docs/pilot_pathway_target_panel.tsv'
SCORE_V3 = PROJECT / 'scores/spatiallr_trust_score_pilot_v3.tsv'
OUT_DIR = PROJECT / 'results/task_e_target_activation'
RUN_DIR = PROJECT / 'runs' / RUN_ID
COMPARTMENTS = ['tumor_like', 'immune', 'stroma_fibroblast', 'endothelial', 'ambiguous_or_low_signal']
KEY_FIELDS = ['dataset', 'sample_id', 'cancer', 'sender_compartment', 'receiver_compartment', 'ligand', 'receptor', 'pathway']


def rel(path: Path) -> str:
    return str(path.relative_to(PROJECT))


def open_text(path: Path):
    if path.suffix == '.gz':
        return gzip.open(path, 'rt', errors='replace')
    return path.open('rt', errors='replace')


def read_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write_tsv(path: Path, rows: Iterable[Dict[str, object]], fields: List[str]) -> None:
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', extrasaction='ignore')
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


def read_target_panel() -> Tuple[Dict[str, List[str]], Set[str]]:
    by_pathway: Dict[str, List[str]] = defaultdict(list)
    wanted: Set[str] = set()
    for row in read_tsv(TARGET_PANEL):
        gene = row['gene'].upper()
        pathway = row['pathway']
        if gene not in by_pathway[pathway]:
            by_pathway[pathway].append(gene)
        wanted.add(gene)
    return dict(by_pathway), wanted


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


def read_mtx_targets(path: Path, index_to_gene: Dict[int, str], barcode_to_comp: Dict[int, str]):
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
            gene = index_to_gene.get(gene_idx)
            if not gene:
                continue
            spot_idx = int(fields[1]) - 1
            comp = barcode_to_comp.get(spot_idx)
            if comp not in COMPARTMENTS:
                continue
            value = float(fields[2])
            gene_comp_sum[gene][comp] += value
            gene_comp_spots[gene][comp].add(spot_idx)
    return shape, gene_comp_sum, gene_comp_spots


def choose_samples(rows: List[Dict[str, str]], sample_ids: List[str], max_samples: int) -> List[Dict[str, str]]:
    ready = [row for row in rows if row.get('status', 'ready') == 'ready']
    if sample_ids:
        wanted = set(sample_ids)
        selected = [row for row in ready if row['sample_id'] in wanted]
        missing = sorted(wanted - {row['sample_id'] for row in selected})
        if missing:
            raise SystemExit('Missing sample_id(s): ' + ', '.join(missing))
        return selected
    return ready[:max_samples]


def process_sample(row: Dict[str, str], annot_row: Dict[str, str], pathway_targets: Dict[str, List[str]], wanted: Set[str]):
    barcodes = read_barcodes(PROJECT / row['barcodes_path'])
    annotations = read_annotations(PROJECT / annot_row['annotation_path'])
    comp_counts = Counter(annotations.get(barcode, 'missing') for barcode in barcodes)
    barcode_to_comp = {i: annotations.get(barcode) for i, barcode in enumerate(barcodes)}
    index_to_gene = read_features(PROJECT / row['features_path'], wanted)
    shape, gene_comp_sum, gene_comp_spots = read_mtx_targets(PROJECT / row['expression_path'], index_to_gene, barcode_to_comp)
    out = []
    for pathway, genes in sorted(pathway_targets.items()):
        detected = [gene for gene in genes if gene in {g for g in gene_comp_sum.keys()}]
        detected_n = len(detected)
        comp_means = {}
        comp_pct = {}
        for comp in COMPARTMENTS:
            spots = comp_counts.get(comp, 0)
            denom = max(spots * max(detected_n, 1), 1)
            total = sum(gene_comp_sum[gene][comp] for gene in detected)
            positive = set()
            for gene in detected:
                positive.update(gene_comp_spots[gene][comp])
            comp_means[comp] = total / denom
            comp_pct[comp] = len(positive) / spots if spots else 0.0
        sorted_vals = sorted(comp_means.values())
        for comp in COMPARTMENTS:
            other_spots = sum(comp_counts.get(c, 0) for c in COMPARTMENTS if c != comp)
            other_total = sum(gene_comp_sum[gene][c] for gene in detected for c in COMPARTMENTS if c != comp)
            other_denom = max(other_spots * max(detected_n, 1), 1)
            other_mean = other_total / other_denom
            rank = sum(1 for val in sorted_vals if val <= comp_means[comp])
            percentile = rank / len(sorted_vals) if sorted_vals else 0.0
            enrichment = (comp_means[comp] + 1e-9) / (other_mean + 1e-9)
            support = 1 if detected_n >= 2 and percentile >= 0.75 and enrichment >= 1.25 else 0
            out.append({
                'dataset': row['dataset'],
                'sample_id': row['sample_id'],
                'cancer': row['cancer'],
                'pathway': pathway,
                'receiver_compartment': comp,
                'target_genes_in_panel': len(genes),
                'target_genes_detected': detected_n,
                'target_genes_detected_list': ','.join(detected),
                'receiver_spots': comp_counts.get(comp, 0),
                'receiver_target_mean': f'{comp_means[comp]:.8g}',
                'other_compartment_target_mean': f'{other_mean:.8g}',
                'receiver_target_enrichment': f'{enrichment:.8g}',
                'receiver_target_positive_spot_fraction': f'{comp_pct[comp]:.8g}',
                'receiver_target_percentile': f'{percentile:.8g}',
                'target_activation_support': support,
                'matrix_rows': shape[0] if shape else '',
                'matrix_cols': shape[1] if shape else '',
                'matrix_nnz': shape[2] if shape else '',
            })
    return out


def join_score(score_rows: List[Dict[str, str]], activity_rows: List[Dict[str, object]]):
    activity = {(r['dataset'], r['sample_id'], r['pathway'], r['receiver_compartment']): r for r in activity_rows}
    out = []
    for row in score_rows:
        key = (row['dataset'], row['sample_id'], row['pathway'], row['receiver_compartment'])
        act = activity.get(key, {})
        support = int(act.get('target_activation_support', 0) or 0)
        new = dict(row)
        for field in ['target_genes_in_panel', 'target_genes_detected', 'target_genes_detected_list', 'receiver_target_mean', 'other_compartment_target_mean', 'receiver_target_enrichment', 'receiver_target_positive_spot_fraction', 'receiver_target_percentile']:
            new[field] = act.get(field, '')
        new['target_activation_support'] = support
        new['target_activation_note'] = 'Pilot receiver-compartment target expression proxy; not validated pathway activation.'
        out.append(new)
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--max-samples', type=int, default=999)
    parser.add_argument('--sample-id', action='append', default=[])
    parser.add_argument('--output-prefix', default='target_activation_all')
    args = parser.parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    manifest = choose_samples(read_tsv(MANIFEST), args.sample_id, args.max_samples)
    annot_rows = {(r['dataset'], r['sample_id']): r for r in read_tsv(ANNOT_MANIFEST)}
    pathway_targets, wanted = read_target_panel()
    activity_rows = []
    failures = []
    for row in manifest:
        try:
            annot = annot_rows[(row['dataset'], row['sample_id'])]
            activity_rows.extend(process_sample(row, annot, pathway_targets, wanted))
        except Exception as exc:
            failures.append({'dataset': row.get('dataset', ''), 'sample_id': row.get('sample_id', ''), 'error': repr(exc)})
    score_rows = read_tsv(SCORE_V3)
    selected_samples = {r['sample_id'] for r in manifest}
    if args.sample_id or args.max_samples < 999:
        score_rows = [r for r in score_rows if r['sample_id'] in selected_samples]
    joined_rows = join_score(score_rows, activity_rows)
    activity_path = OUT_DIR / f'{args.output_prefix}_activity_by_receiver.tsv'
    joined_path = OUT_DIR / f'{args.output_prefix}_score_v3_with_target.tsv'
    failure_path = OUT_DIR / f'{args.output_prefix}_failures.tsv'
    summary_path = OUT_DIR / f'{args.output_prefix}_summary.json'
    activity_fields = ['dataset', 'sample_id', 'cancer', 'pathway', 'receiver_compartment', 'target_genes_in_panel', 'target_genes_detected', 'target_genes_detected_list', 'receiver_spots', 'receiver_target_mean', 'other_compartment_target_mean', 'receiver_target_enrichment', 'receiver_target_positive_spot_fraction', 'receiver_target_percentile', 'target_activation_support', 'matrix_rows', 'matrix_cols', 'matrix_nnz']
    write_tsv(activity_path, activity_rows, activity_fields)
    joined_fields = list(score_rows[0].keys()) + ['target_genes_in_panel', 'target_genes_detected', 'target_genes_detected_list', 'receiver_target_mean', 'other_compartment_target_mean', 'receiver_target_enrichment', 'receiver_target_positive_spot_fraction', 'receiver_target_percentile', 'target_activation_support', 'target_activation_note'] if score_rows else []
    write_tsv(joined_path, joined_rows, joined_fields)
    write_tsv(failure_path, failures, ['dataset', 'sample_id', 'error'])
    high_rows = [r for r in joined_rows if r.get('confidence_tier_pilot_v3') == 'high_pilot_v3']
    summary = {
        'run_id': RUN_ID,
        'samples_requested': len(manifest),
        'samples_failed': len(failures),
        'activity_rows': len(activity_rows),
        'joined_score_rows': len(joined_rows),
        'target_supported_rows': sum(int(r.get('target_activation_support', 0) or 0) for r in joined_rows),
        'high_rows': len(high_rows),
        'high_target_supported_rows': sum(int(r.get('target_activation_support', 0) or 0) for r in high_rows),
        'outputs': {
            'activity_by_receiver': rel(activity_path),
            'score_v3_with_target': rel(joined_path),
            'failures': rel(failure_path),
            'summary': rel(summary_path),
        },
        'note': 'Pilot receiver-compartment target expression proxy. This is not validated pathway activation or experimental evidence.',
    }
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + '\n')
    manifest_path = RUN_DIR / 'outputs_manifest.tsv'
    with manifest_path.open('w') as handle:
        handle.write('path\ttype\tdescription\n')
        for path, typ, desc in [
            (activity_path, 'tsv', 'Pathway target activity proxy by sample, pathway, and receiver compartment'),
            (joined_path, 'tsv', 'Pilot score v3 rows with receiver target activation proxy columns'),
            (failure_path, 'tsv', 'Target activation sample failures'),
            (summary_path, 'json', 'Target activation summary'),
            (manifest_path, 'tsv', 'Run output manifest'),
        ]:
            handle.write(f'{rel(path)}\t{typ}\t{desc}\n')
    print(json.dumps(summary, indent=2, sort_keys=True))

if __name__ == '__main__':
    main()
