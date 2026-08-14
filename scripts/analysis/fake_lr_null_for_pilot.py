#!/usr/bin/env python3
"""Expression-matched fake ligand-receptor null for stdlib_marker_lr_pilot."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import random
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

PROJECT = Path(os.environ.get('PROJECT', str(Path(__file__).resolve().parents[2])))
RUN_ID = os.environ.get('RUN_ID', 'manual')
CANDIDATES = PROJECT / 'results/task_c_pilot_baseline/stdlib_marker_lr_pilot_candidates.tsv'
LR_PANEL = PROJECT / 'docs/pilot_lr_panel.tsv'
OUT_DIR = PROJECT / 'results/task_d_null_models'
RUN_DIR = PROJECT / 'runs' / RUN_ID


def rel(path: Path) -> str:
    return str(path.relative_to(PROJECT))


def read_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write_tsv(path: Path, rows: List[Dict[str, str]], fields: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def build_real_pairs() -> set[Tuple[str, str]]:
    return {(r['ligand'].upper(), r['receptor'].upper()) for r in read_tsv(LR_PANEL)}


def nearest_genes(pool: Dict[str, float], target_gene: str, target_mean: float) -> List[str]:
    genes = [g for g in pool if g != target_gene]
    genes.sort(key=lambda g: (abs(pool[g] - target_mean), g))
    return genes


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--draws', type=int, default=100)
    parser.add_argument('--seed', type=int, default=20260709)
    args = parser.parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    rows = read_tsv(CANDIDATES)
    real_pairs = build_real_pairs()
    ligand_pool = defaultdict(dict)   # (sample,sender)->ligand mean
    receptor_pool = defaultdict(dict) # (sample,receiver)->receptor mean
    for r in rows:
        ligand_pool[(r['sample_id'], r['sender_compartment'])][r['ligand']] = float(r['ligand_mean_sender'])
        receptor_pool[(r['sample_id'], r['receiver_compartment'])][r['receptor']] = float(r['receptor_mean_receiver'])
    rng = random.Random(args.seed)
    out_rows = []
    skipped = 0
    for idx, r in enumerate(rows):
        sample = r['sample_id']
        sender = r['sender_compartment']
        receiver = r['receiver_compartment']
        ligand = r['ligand']
        receptor = r['receptor']
        ligand_mean = float(r['ligand_mean_sender'])
        receptor_mean = float(r['receptor_mean_receiver'])
        adjacency = float(r['adjacency_pairs'])
        observed = float(r['pilot_score'])
        ligands = nearest_genes(ligand_pool[(sample, sender)], ligand, ligand_mean)
        receptors = nearest_genes(receptor_pool[(sample, receiver)], receptor, receptor_mean)
        fake_scores = []
        fake_pairs = []
        attempts = 0
        while len(fake_scores) < args.draws and attempts < args.draws * 50:
            attempts += 1
            if not ligands or not receptors:
                break
            # Bias toward expression-matched alternatives by sampling from the closest 5 when available.
            l_choices = ligands[: min(5, len(ligands))]
            r_choices = receptors[: min(5, len(receptors))]
            fake_l = rng.choice(l_choices)
            fake_r = rng.choice(r_choices)
            if (fake_l, fake_r) in real_pairs:
                continue
            fake_l_mean = ligand_pool[(sample, sender)][fake_l]
            fake_r_mean = receptor_pool[(sample, receiver)][fake_r]
            fake_scores.append(fake_l_mean * fake_r_mean * math.log1p(adjacency))
            fake_pairs.append((fake_l, fake_r))
        if not fake_scores:
            skipped += 1
            mean_score = 0.0
            ge = 0
            empirical_p = ''
            ratio = ''
            fake_n = 0
        else:
            fake_n = len(fake_scores)
            mean_score = sum(fake_scores) / fake_n
            ge = sum(1 for s in fake_scores if s >= observed)
            empirical_p = (ge + 1) / (fake_n + 1)
            ratio = observed / mean_score if mean_score > 0 else ''
        out = dict(r)
        out.update({
            'null_model': 'expression_matched_fake_lr_pair',
            'fake_draws_requested': str(args.draws),
            'fake_draws_used': str(fake_n),
            'null_mean_score': f'{mean_score:.8g}',
            'null_ge_observed_n': str(ge),
            'empirical_p_ge': f'{empirical_p:.8g}' if empirical_p != '' else '',
            'observed_to_null_mean_ratio': f'{ratio:.8g}' if ratio != '' else '',
            'fake_pair_examples': ';'.join(f'{a}-{b}' for a,b in fake_pairs[:5]),
        })
        out_rows.append(out)
    fields = list(out_rows[0].keys())
    out_path = OUT_DIR / 'stdlib_marker_lr_fake_pair_null.tsv'
    write_tsv(out_path, out_rows, fields)
    pvals = [float(r['empirical_p_ge']) for r in out_rows if r['empirical_p_ge']]
    summary = {
        'run_id': RUN_ID,
        'null_model': 'expression_matched_fake_lr_pair',
        'source_method': 'stdlib_marker_lr_pilot',
        'candidate_rows': len(out_rows),
        'draws_requested': args.draws,
        'skipped_no_fake_pairs': skipped,
        'p_le_0_05': sum(1 for p in pvals if p <= 0.05),
        'p_le_0_10': sum(1 for p in pvals if p <= 0.10),
        'output': rel(out_path),
    }
    (OUT_DIR / 'stdlib_marker_lr_fake_pair_null_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    (RUN_DIR / 'fake_lr_null_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
