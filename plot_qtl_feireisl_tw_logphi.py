#!/usr/bin/env python3
"""Plot log(phi) vs log(epsilon) for QLT, Fereisl, and TW results.

This script is designed to work with CSV output from saturday.cpp / 4_5pn_comparison.
It supports:
- a single CSV file with epsilon and one or more phi columns
- separate CSV files for QLT, Fereisl, and TW

If current CSV data only contains one phi series, the script still plots it and
prints a warning for any missing solutions.
"""

import argparse
import csv
import math
import os
import sys

import matplotlib.pyplot as plt
import numpy as np


def load_csv(filename):
    with open(filename, newline='') as fp:
        reader = csv.reader(fp)
        header = next(reader)
        rows = [row for row in reader if row and any(field.strip() for field in row)]
    return header, rows


def find_column(header, candidates):
    normalized = [h.strip().lower() for h in header]
    for candidate in candidates:
        candidate_lower = candidate.lower()
        if candidate_lower in normalized:
            return normalized.index(candidate_lower)
    for candidate in candidates:
        candidate_lower = candidate.lower()
        for idx, h in enumerate(normalized):
            if h.endswith(candidate_lower):
                return idx
            if candidate_lower in h:
                return idx
    return None


def parse_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def extract_series(header, rows, eps_col, phi_col):
    eps = []
    phi = []
    for row in rows:
        if len(row) <= max(eps_col, phi_col):
            continue
        eps_val = parse_float(row[eps_col])
        phi_val = parse_float(row[phi_col])
        if eps_val > 0 and phi_val > 0 and math.isfinite(eps_val) and math.isfinite(phi_val):
            eps.append(eps_val)
            phi.append(phi_val)
    return np.array(eps), np.array(phi)


def build_series_from_file(filename, eps_col=None, phi_col=None, phi_name=None):
    header, rows = load_csv(filename)
    normalized = [h.strip().lower() for h in header]

    if eps_col is not None:
        eps_col_idx = normalized.index(eps_col.lower()) if eps_col.lower() in normalized else None
    else:
        eps_col_idx = find_column(header, ['eps', 'epsilon', 'e'])

    if eps_col_idx is None:
        raise ValueError(f"Could not find an epsilon column in {filename}. Header: {header}")

    if phi_col is not None:
        phi_col_idx = normalized.index(phi_col.lower()) if phi_col.lower() in normalized else None
    else:
        possible_phi_cols = [phi_name or 'phi', 'phi_total', 'phase', 'phi_qlt', 'phi_fereisl', 'phi_tw']
        phi_col_idx = find_column(header, possible_phi_cols)

    if phi_col_idx is None:
        raise ValueError(f"Could not find a phi column in {filename}. Header: {header}")

    eps, phi = extract_series(header, rows, eps_col_idx, phi_col_idx)
    if eps.size == 0:
        raise ValueError(f"No valid positive epsilon/phi pairs found in {filename}.")

    label = phi_name or header[phi_col_idx].strip() or os.path.basename(filename)
    return label, eps, phi


def build_plot(series_list, output_file, title):
    plt.figure(figsize=(10, 7))
    for label, eps, phi in series_list:
        log_eps = np.log10(eps)
        log_phi = np.log10(phi)
        sorted_idx = np.argsort(log_eps)
        plt.plot(log_eps[sorted_idx], log_phi[sorted_idx], 'o-', linewidth=2, markersize=6, label=label)

    plt.xlabel('log10(\u03B5)', fontsize=14)
    plt.ylabel('log10(\u03C6)', fontsize=14)
    plt.title(title, fontsize=15)
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=11)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved log(phi) vs log(eps) plot to {output_file}")
    plt.show()


def main():
    parser = argparse.ArgumentParser(
        description='Plot log(phi) vs log(epsilon) for QLT, Fereisl, and TW solutions.'
    )
    parser.add_argument('--csv', default='fereisl_tw_eps_phi_convergence.csv',
                        help='Single CSV file containing eps and phi columns')
    parser.add_argument('--qlt-file', help='CSV file containing QLT eps/phi data')
    parser.add_argument('--fereisl-file', help='CSV file containing Fereisl eps/phi data')
    parser.add_argument('--tw-file', help='CSV file containing TW eps/phi data')
    parser.add_argument('--eps-col', help='epsilon column name in the CSV file(s)')
    parser.add_argument('--phi-col', help='phi column name in the CSV file(s)')
    parser.add_argument('--output', default='qtl_fereisl_tw_logphi.png',
                        help='Output image filename')
    parser.add_argument('--title', default='log(phi) vs log(epsilon) for QLT, Fereisl, and TW',
                        help='Plot title')
    args = parser.parse_args()

    series_list = []
    if args.qlt_file:
        if not os.path.isfile(args.qlt_file):
            parser.error(f"QLT file not found: {args.qlt_file}")
        series_list.append(build_series_from_file(args.qlt_file, eps_col=args.eps_col, phi_col=args.phi_col, phi_name='QLT'))

    if args.fereisl_file:
        if not os.path.isfile(args.fereisl_file):
            parser.error(f"Fereisl file not found: {args.fereisl_file}")
        series_list.append(build_series_from_file(args.fereisl_file, eps_col=args.eps_col, phi_col=args.phi_col, phi_name='Fereisl'))

    if args.tw_file:
        if not os.path.isfile(args.tw_file):
            parser.error(f"TW file not found: {args.tw_file}")
        series_list.append(build_series_from_file(args.tw_file, eps_col=args.eps_col, phi_col=args.phi_col, phi_name='TW'))

    if not series_list:
        if not os.path.isfile(args.csv):
            parser.error(f"CSV file not found: {args.csv}")
        header, _ = load_csv(args.csv)
        normalized = [h.strip().lower() for h in header]
        eps_col = args.eps_col if args.eps_col else find_column(header, ['eps', 'epsilon', 'e'])
        if eps_col is None:
            parser.error(f"Could not find an epsilon column in {args.csv}. Header: {header}")

        possible_phi_cols = ['phi_qlt', 'phi_qtl', 'qlt', 'phi_fereisl', 'fereisl', 'phi_tw', 'tw', 'phi', 'phase']
        found = set()
        for name in possible_phi_cols:
            idx = find_column(header, [name])
            if idx is not None and idx not in found:
                found.add(idx)
                label = header[idx].strip() or f'phi_{len(found)}'
                eps, phi = extract_series(header, load_csv(args.csv)[1], normalized.index(eps_col.lower()) if isinstance(eps_col, str) else eps_col, idx)
                if eps.size > 0:
                    series_list.append((label, eps, phi))

        if not series_list:
            parser.error(
                f"Could not find any phi columns in {args.csv}. Header: {header}\n"
                "Try passing --phi-col or separate --qlt-file/--fereisl-file/--tw-file."
            )

    if len(series_list) == 0:
        parser.error('No data series found to plot.')

    if len(series_list) < 3:
        print("Warning: fewer than three solution series were found."
              " Use --qlt-file, --fereisl-file, and --tw-file if you have separate CSVs.")

    build_plot(series_list, args.output, args.title)


if __name__ == '__main__':
    main()
