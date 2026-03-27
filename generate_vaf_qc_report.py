#!/usr/bin/env python3
"""
VAF Quality Control Report Generator for Repun Pipeline

Generates a grouped QC report with AF distributions organized by gene/mutation,
comparing across platforms, assays, coverage depths, labs, and batches.
Samples are CRISPR-edited HG002 clones with known somatic variants.

This project was developed with the assistance of Claude (Anthropic).
"""

import os
import sys
import gzip
import json
import re
from pathlib import Path
from collections import defaultdict
from statistics import mean, median, stdev

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Gene patterns sorted longest-first so greedy match works
KNOWN_GENES = [
    'EGFR_A763_Y764insFQEA',
    'PDGFRA_I843del',
    'ERBB2-V659E', 'ERBB2_V659E',
    'FGFR3-S249C', 'FGFR3_S249C',
    'BRAF-V600E', 'BRAF_V600E',
    'CCDC6-RET',
    'TPM3-NTRK1',
    'PDGFRA',
    'ERBB2',
    'EGFR',
]

GENE_CANONICAL = {
    'BRAF-V600E': 'BRAF-V600E', 'BRAF_V600E': 'BRAF-V600E',
    'CCDC6-RET': 'CCDC6-RET',
    'EGFR_A763_Y764insFQEA': 'EGFR', 'EGFR': 'EGFR',
    'ERBB2-V659E': 'ERBB2-V659E', 'ERBB2_V659E': 'ERBB2-V659E',
    'ERBB2': 'ERBB2-V659E',
    'FGFR3-S249C': 'FGFR3-S249C', 'FGFR3_S249C': 'FGFR3-S249C',
    'PDGFRA_I843del': 'PDGFRA', 'PDGFRA': 'PDGFRA',
    'TPM3-NTRK1': 'TPM3-NTRK1',
}

# Chart.js color palette (10 distinct colors with alpha)
PALETTE = [
    ('31,119,180',),   # blue
    ('255,127,14',),   # orange
    ('44,160,44',),    # green
    ('214,39,40',),    # red
    ('148,103,189',),  # purple
    ('140,86,75',),    # brown
    ('227,119,194',),  # pink
    ('127,127,127',),  # gray
    ('188,189,34',),   # olive
    ('23,190,207',),   # cyan
]

AF_BIN_EDGES = [0, 1, 3, 5, 8, 10, 15, 20, 25, 50, 100]
AF_BIN_LABELS = [f'{AF_BIN_EDGES[i]}-{AF_BIN_EDGES[i+1]}%'
                 for i in range(len(AF_BIN_EDGES) - 1)]

# ---------------------------------------------------------------------------
# VCF parsing
# ---------------------------------------------------------------------------

def extract_vaf_from_vcf(vcf_path):
    """
    Extract VAF values and FILTER status from a VCF file.
    Returns (pass_vafs, lowvaf_vafs, total_variants) where VAFs are in percent.
    """
    pass_vafs = []
    lowvaf_vafs = []
    total = 0

    try:
        opener = gzip.open if vcf_path.endswith('.gz') else open
        with opener(vcf_path, 'rt') as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                fields = line.rstrip('\n').split('\t')
                if len(fields) < 10:
                    continue
                total += 1

                filt = fields[6]
                fmt = fields[8].split(':')
                sample = fields[9].split(':')

                try:
                    af_idx = fmt.index('AF')
                except ValueError:
                    continue
                if af_idx >= len(sample):
                    continue
                af_str = sample[af_idx]
                if not af_str or af_str == '.':
                    continue
                try:
                    af_pct = float(af_str) * 100
                except ValueError:
                    continue

                if filt == 'PASS':
                    pass_vafs.append(af_pct)
                else:
                    lowvaf_vafs.append(af_pct)

    except Exception as e:
        print(f"Error reading {vcf_path}: {e}", file=sys.stderr)
        return [], [], 0

    return pass_vafs, lowvaf_vafs, total


def calculate_stats(values):
    if not values:
        return None
    s = sorted(values)
    n = len(s)
    return {
        'count': n,
        'min': min(s),
        'max': max(s),
        'mean': mean(s),
        'median': median(s),
        'stdev': stdev(s) if n > 1 else 0,
        'q25': s[int(n * 0.25)],
        'q75': s[int(n * 0.75)],
    }


def bin_vafs(vafs):
    """Bin VAF values into AF_BIN_EDGES histogram."""
    hist = [0] * (len(AF_BIN_EDGES) - 1)
    for v in vafs:
        for i in range(len(AF_BIN_EDGES) - 1):
            if AF_BIN_EDGES[i] <= v < AF_BIN_EDGES[i + 1]:
                hist[i] += 1
                break
        else:
            if v >= AF_BIN_EDGES[-1]:
                hist[-1] += 1
    return hist

# ---------------------------------------------------------------------------
# Metadata parsing
# ---------------------------------------------------------------------------

def parse_metadata(sample_name):
    """Parse experimental metadata from the results subdirectory name."""
    meta = {
        'sample_name': sample_name,
        'gene': None,
        'platform': 'Illumina',
        'assay': None,
        'coverage': None,
        'lab': 'Broad',
        'batch': None,
        'replicate': None,
        'sample_type': 'edited',
        'prep': 'Fresh',
        'passage': None,
        'label': '',
    }

    name = sample_name

    # --- Controls ---
    if name.startswith('HG002_'):
        meta['sample_type'] = 'control'
        meta['gene'] = 'HG002'
        rest = name[6:]  # after HG002_
        if 'FFPE' in rest:
            meta['prep'] = 'FFPE'
        if 'RNA' in rest:
            meta['assay'] = 'RNA'
        elif 'WGS' in rest:
            meta['assay'] = 'WGS'
        m = re.search(r'(\d+)x', rest)
        if m:
            meta['coverage'] = m.group(0)
        meta['label'] = f"HG002 {meta['assay'] or ''} {meta['coverage'] or ''} {meta['prep']}".strip()
        return meta

    if name.startswith('PAR_DNA_'):
        meta['sample_type'] = 'control'
        meta['gene'] = 'PARENT-DNA'
        rest = name[8:]
        if 'WGS' in rest:
            meta['assay'] = 'WGS'
        m = re.search(r'(\d+)x', rest)
        if m:
            meta['coverage'] = m.group(0)
        if 'Azenta' in rest:
            meta['lab'] = 'Azenta'
        meta['label'] = f"Parent DNA {meta['assay'] or ''} {meta['coverage'] or ''} {meta['lab']}"
        return meta

    # --- FFPE Mixtures (P5_, P15_) ---
    pm = re.match(r'P(\d+)_FFPE_(.+)', name)
    if pm:
        meta['sample_type'] = 'mixture'
        meta['prep'] = 'FFPE'
        meta['passage'] = f'P{pm.group(1)}'
        meta['gene'] = f'SRSmix-{meta["passage"]}'
        rest = pm.group(2)
        if 'WGS' in rest:
            meta['assay'] = 'WGS'
        elif 'WES' in rest:
            meta['assay'] = 'WES'
        elif 'RNA' in rest:
            meta['assay'] = 'RNA'
        m = re.search(r'(\d+)x', rest)
        if m:
            meta['coverage'] = m.group(0)
        # Check for replicate number at end
        rm = re.search(r'_(\d+)$', rest)
        if rm:
            meta['replicate'] = rm.group(1)
        meta['label'] = f"{meta['passage']} FFPE {meta['assay'] or ''} {meta['coverage'] or ''}"
        if meta['replicate']:
            meta['label'] += f" Rep{meta['replicate']}"
        return meta

    # --- Strip numeric prefix (Azenta RNA samples like 64184_CCDC6-RET_RNA_Azenta) ---
    if name[0].isdigit():
        m = re.match(r'\d+_(.*)', name)
        if m:
            name = m.group(1)

    # --- Gene-specific samples ---
    # Try to match a known gene name at the start
    matched_gene = None
    remainder = name
    for gene_pat in KNOWN_GENES:
        if name.startswith(gene_pat):
            matched_gene = gene_pat
            remainder = name[len(gene_pat):]
            if remainder.startswith('_'):
                remainder = remainder[1:]
            break

    if matched_gene:
        meta['gene'] = GENE_CANONICAL[matched_gene]
        meta['sample_type'] = 'edited'
    else:
        # Fallback: treat whole name as unknown
        meta['gene'] = name
        meta['label'] = name
        return meta

    # Parse remainder for assay, coverage, lab, batch, platform
    # Check for PacBio
    if 'broad_pacbio' in sample_name or 'pacbio' in remainder.lower():
        meta['platform'] = 'PacBio'
        meta['lab'] = 'Broad'

    parts = remainder.split('_')
    for p in parts:
        p_up = p.upper()
        if p_up in ('WGS', 'WES'):
            meta['assay'] = p_up
        elif p_up == 'RNA':
            meta['assay'] = 'RNA'
        elif re.match(r'^\d+x$', p, re.IGNORECASE):
            meta['coverage'] = p.lower()
        elif p_up == 'AZENTA':
            meta['lab'] = 'Azenta'
        elif p_up == 'BROAD':
            meta['lab'] = 'Broad'
        elif p_up == 'PACBIO':
            meta['platform'] = 'PacBio'
        elif re.match(r'^B\d+$', p, re.IGNORECASE):
            meta['batch'] = p.upper()
        elif re.match(r'^\d+$', p) and meta['batch']:
            meta['replicate'] = p

    # Build human-readable label
    parts_label = []
    parts_label.append(meta['platform'])
    if meta['coverage']:
        parts_label.append(meta['coverage'])
    if meta['assay']:
        parts_label.append(meta['assay'])
    parts_label.append(meta['lab'])
    if meta['batch']:
        parts_label.append(meta['batch'])
    if meta['replicate']:
        parts_label.append(f'r{meta["replicate"]}')
    meta['label'] = ' '.join(parts_label)

    return meta

# ---------------------------------------------------------------------------
# Sample processing
# ---------------------------------------------------------------------------

def process_all_samples(results_dir):
    """Process all sample result directories."""
    samples = {}
    results_path = Path(results_dir)

    for sample_dir in sorted(results_path.iterdir()):
        if not sample_dir.is_dir() or sample_dir.name in ('pipeline_info', 'qc_reports', 's3_cache'):
            continue

        name = sample_dir.name
        repun_dirs = list(sample_dir.glob('repun_*'))
        if not repun_dirs:
            continue
        repun_dir = repun_dirs[0]

        meta = parse_metadata(name)

        sample_data = {
            'meta': meta,
            'repun_dir': str(repun_dir),
            'vcfs': {},
        }

        for vcf_type, vcf_name in [('truths', 'unified_truths.vcf.gz'),
                                    ('candidates', 'unified_candidates.vcf.gz')]:
            vcf_path = repun_dir / vcf_name
            if vcf_path.exists():
                pass_vafs, lowvaf_vafs, total = extract_vaf_from_vcf(str(vcf_path))
                sample_data['vcfs'][vcf_type] = {
                    'exists': True,
                    'total_variants': total,
                    'pass_count': len(pass_vafs),
                    'lowvaf_count': len(lowvaf_vafs),
                    'pass_vafs': pass_vafs,
                    'lowvaf_vafs': lowvaf_vafs,
                    'pass_stats': calculate_stats(pass_vafs),
                    'pass_hist': bin_vafs(pass_vafs),
                }
            else:
                sample_data['vcfs'][vcf_type] = {'exists': False}

        samples[name] = sample_data

    return samples


def group_samples(samples):
    """Group samples into gene groups, mixtures, and controls."""
    gene_groups = defaultdict(list)
    mixture_groups = defaultdict(list)
    control_list = []

    for name, data in sorted(samples.items()):
        st = data['meta']['sample_type']
        gene = data['meta']['gene']
        if st == 'edited':
            gene_groups[gene].append(data)
        elif st == 'mixture':
            mixture_groups[gene].append(data)
        elif st == 'control':
            control_list.append(data)

    return dict(sorted(gene_groups.items())), dict(sorted(mixture_groups.items())), control_list

# ---------------------------------------------------------------------------
# HTML report generation
# ---------------------------------------------------------------------------

def _stat_cell(vcf_info, field='truths'):
    """Return formatted stats for a table cell."""
    vi = vcf_info.get(field, {})
    if not vi.get('exists'):
        return '-', '-', '-', 0, 0, 0
    stats = vi.get('pass_stats')
    if not stats:
        return '0', '-', '-', vi.get('total_variants', 0), 0, vi.get('lowvaf_count', 0)
    return (
        str(stats['count']),
        f"{stats['mean']:.1f}%",
        f"{stats['median']:.1f}%",
        vi.get('total_variants', 0),
        stats['count'],
        vi.get('lowvaf_count', 0),
    )


def generate_html_report(samples, output_path):
    gene_groups, mixture_groups, controls = group_samples(samples)

    html = ["""<!DOCTYPE html>
<html><head><meta charset="utf-8">
<title>Repun VAF QC Report</title>
<script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
<style>
body { font-family: 'Segoe UI', Arial, sans-serif; margin: 0; background: #f0f2f5; color: #1a1a2e; }
.container { max-width: 1500px; margin: 0 auto; padding: 24px; }
h1 { font-size: 1.6em; border-bottom: 3px solid #2563eb; padding-bottom: 8px; margin-bottom: 4px; }
.subtitle { color: #6b7280; font-size: 0.9em; margin-bottom: 24px; }
h2 { font-size: 1.25em; margin-top: 36px; color: #1e40af; border-left: 4px solid #2563eb; padding-left: 10px; }
h3 { font-size: 1.05em; color: #374151; margin-top: 20px; }
table { width: 100%; border-collapse: collapse; margin: 12px 0 24px 0; font-size: 0.85em; }
th { background: #1e40af; color: white; padding: 8px 10px; text-align: left; position: sticky; top: 0; }
td { padding: 6px 10px; border-bottom: 1px solid #e5e7eb; }
tr:hover td { background: #f0f4ff; }
.gene-header td { background: #dbeafe; font-weight: 700; font-size: 0.95em; border-bottom: 2px solid #93c5fd; }
.num { text-align: right; font-variant-numeric: tabular-nums; }
.chart-section { display: flex; flex-wrap: wrap; gap: 20px; margin: 16px 0 32px 0; }
.chart-box { flex: 1 1 620px; min-width: 500px; max-width: 760px; background: white;
  border-radius: 8px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); padding: 16px; }
.chart-box h4 { margin: 0 0 8px 0; font-size: 0.95em; color: #374151; }
.chart-box canvas { width: 100% !important; height: 320px !important; }
.section { background: white; border-radius: 8px; box-shadow: 0 1px 3px rgba(0,0,0,0.1);
  padding: 20px; margin-bottom: 24px; }
.footer { color: #9ca3af; font-size: 0.75em; margin-top: 32px; text-align: center; }
.badge { display: inline-block; padding: 2px 7px; border-radius: 3px; font-size: 0.78em; font-weight: 600; }
.badge-wgs { background: #dbeafe; color: #1e40af; }
.badge-wes { background: #d1fae5; color: #065f46; }
.badge-rna { background: #fce7f3; color: #9d174d; }
.badge-pacbio { background: #fef3c7; color: #92400e; }
.badge-ilmn { background: #e0e7ff; color: #3730a3; }
.tabs { display: flex; gap: 4px; margin-bottom: 8px; }
.tab { padding: 4px 12px; border: 1px solid #d1d5db; border-radius: 4px 4px 0 0;
  cursor: pointer; font-size: 0.82em; background: #f9fafb; }
.tab.active { background: #2563eb; color: white; border-color: #2563eb; }
</style>
</head><body><div class="container">
<h1>Repun Pipeline VAF QC Report</h1>
<p class="subtitle">AF distributions for CRISPR-edited HG002 clones &mdash; PASS variants only (LowVAF excluded from distributions)</p>
"""]

    # --- Summary table ---
    html.append('<div class="section"><h2>Summary by Gene / Mutation</h2>')
    html.append("""<table>
<thead><tr>
  <th>Gene</th><th>Condition</th><th>Platform</th><th>Assay</th><th>Coverage</th>
  <th>Lab</th><th>Batch</th>
  <th class="num">Total</th><th class="num">PASS</th><th class="num">LowVAF</th>
  <th class="num">Mean AF</th><th class="num">Median AF</th>
</tr></thead><tbody>""")

    for gene, group in gene_groups.items():
        html.append(f'<tr class="gene-header"><td colspan="12">{gene}</td></tr>')
        for s in group:
            m = s['meta']
            t = s['vcfs'].get('truths', {})
            ps = t.get('pass_stats')
            assay_badge = {'WGS': 'wgs', 'WES': 'wes', 'RNA': 'rna'}.get(m['assay'], 'wgs')
            plat_badge = 'pacbio' if m['platform'] == 'PacBio' else 'ilmn'
            html.append(f'''<tr>
  <td></td>
  <td>{m["label"]}</td>
  <td><span class="badge badge-{plat_badge}">{m["platform"]}</span></td>
  <td><span class="badge badge-{assay_badge}">{m["assay"] or "-"}</span></td>
  <td class="num">{m["coverage"] or "-"}</td>
  <td>{m["lab"]}</td>
  <td>{m["batch"] or "-"}</td>
  <td class="num">{t.get("total_variants", "-")}</td>
  <td class="num">{t.get("pass_count", "-")}</td>
  <td class="num">{t.get("lowvaf_count", "-")}</td>
  <td class="num">{f'{ps["mean"]:.1f}%' if ps else "-"}</td>
  <td class="num">{f'{ps["median"]:.1f}%' if ps else "-"}</td>
</tr>''')

    # Mixture rows
    if mixture_groups:
        html.append('<tr class="gene-header"><td colspan="12">SRSmix FFPE Mixtures</td></tr>')
        for gene, group in mixture_groups.items():
            for s in group:
                m = s['meta']
                t = s['vcfs'].get('truths', {})
                ps = t.get('pass_stats')
                assay_cls = {'WGS': 'wgs', 'WES': 'wes', 'RNA': 'rna'}.get(m['assay'], 'wgs')
                mean_s = f'{ps["mean"]:.1f}%' if ps else '-'
                med_s = f'{ps["median"]:.1f}%' if ps else '-'
                html.append(f'''<tr>
  <td>{m["passage"]}</td>
  <td>{m["label"]}</td>
  <td><span class="badge badge-ilmn">Illumina</span></td>
  <td><span class="badge badge-{assay_cls}">{m["assay"] or "-"}</span></td>
  <td class="num">{m["coverage"] or "-"}</td>
  <td>{m["lab"]}</td><td>-</td>
  <td class="num">{t.get("total_variants", "-")}</td>
  <td class="num">{t.get("pass_count", "-")}</td>
  <td class="num">{t.get("lowvaf_count", "-")}</td>
  <td class="num">{mean_s}</td>
  <td class="num">{med_s}</td>
</tr>''')

    # Control rows
    if controls:
        html.append('<tr class="gene-header"><td colspan="12">Controls</td></tr>')
        for s in controls:
            m = s['meta']
            t = s['vcfs'].get('truths', {})
            ps = t.get('pass_stats')
            assay_cls = {'WGS': 'wgs', 'WES': 'wes', 'RNA': 'rna'}.get(m.get('assay'), 'wgs')
            mean_s = f'{ps["mean"]:.1f}%' if ps else '-'
            med_s = f'{ps["median"]:.1f}%' if ps else '-'
            html.append(f'''<tr>
  <td>{m["gene"]}</td>
  <td>{m["label"]}</td>
  <td><span class="badge badge-ilmn">Illumina</span></td>
  <td><span class="badge badge-{assay_cls}">{m.get("assay") or "-"}</span></td>
  <td class="num">{m.get("coverage") or "-"}</td>
  <td>{m["lab"]}</td><td>-</td>
  <td class="num">{t.get("total_variants", "-")}</td>
  <td class="num">{t.get("pass_count", "-")}</td>
  <td class="num">{t.get("lowvaf_count", "-")}</td>
  <td class="num">{mean_s}</td>
  <td class="num">{med_s}</td>
</tr>''')

    html.append('</tbody></table></div>')

    # --- Per-gene AF distribution charts ---
    html.append('<h2>AF Distributions by Gene (PASS Variants Only)</h2>')
    chart_id = 0
    chart_scripts = []

    for gene, group in gene_groups.items():
        html.append(f'<div class="section"><h3>{gene}</h3>')
        html.append('<div class="chart-section">')

        for vcf_type, vcf_label in [('truths', 'Truths'), ('candidates', 'Candidates')]:
            cid = f'chart_{chart_id}'
            chart_id += 1
            html.append(f'<div class="chart-box"><h4>{vcf_label} — AF Distribution</h4>')
            html.append(f'<canvas id="{cid}"></canvas></div>')

            datasets = []
            for i, s in enumerate(group):
                vi = s['vcfs'].get(vcf_type, {})
                if not vi.get('exists') or not vi.get('pass_vafs'):
                    continue
                hist = bin_vafs(vi['pass_vafs'])
                ci = i % len(PALETTE)
                rgb = PALETTE[ci][0]
                datasets.append({
                    'label': s['meta']['label'],
                    'data': hist,
                    'backgroundColor': f'rgba({rgb},0.35)',
                    'borderColor': f'rgba({rgb},1)',
                    'borderWidth': 1,
                })

            chart_scripts.append(f"""
new Chart(document.getElementById('{cid}'), {{
  type: 'bar',
  data: {{ labels: {json.dumps(AF_BIN_LABELS)}, datasets: {json.dumps(datasets)} }},
  options: {{
    responsive: true, maintainAspectRatio: false,
    plugins: {{ legend: {{ position: 'bottom', labels: {{ font: {{ size: 11 }} }} }} }},
    scales: {{
      x: {{ title: {{ display: true, text: 'VAF Range' }} }},
      y: {{ beginAtZero: true, title: {{ display: true, text: 'Variant Count' }} }}
    }}
  }}
}});""")

        html.append('</div></div>')  # chart-section, section

    # --- Mixture section ---
    if mixture_groups:
        html.append('<h2>SRSmix FFPE Mixture AF Distributions</h2>')
        # Group mixtures by assay for comparison
        all_mix = []
        for g, group in mixture_groups.items():
            all_mix.extend(group)

        html.append('<div class="section">')
        html.append('<div class="chart-section">')
        for vcf_type, vcf_label in [('truths', 'Truths'), ('candidates', 'Candidates')]:
            cid = f'chart_{chart_id}'
            chart_id += 1
            html.append(f'<div class="chart-box"><h4>{vcf_label} — AF Distribution</h4>')
            html.append(f'<canvas id="{cid}"></canvas></div>')

            datasets = []
            for i, s in enumerate(all_mix):
                vi = s['vcfs'].get(vcf_type, {})
                if not vi.get('exists') or not vi.get('pass_vafs'):
                    continue
                hist = bin_vafs(vi['pass_vafs'])
                ci = i % len(PALETTE)
                rgb = PALETTE[ci][0]
                datasets.append({
                    'label': s['meta']['label'],
                    'data': hist,
                    'backgroundColor': f'rgba({rgb},0.35)',
                    'borderColor': f'rgba({rgb},1)',
                    'borderWidth': 1,
                })

            chart_scripts.append(f"""
new Chart(document.getElementById('{cid}'), {{
  type: 'bar',
  data: {{ labels: {json.dumps(AF_BIN_LABELS)}, datasets: {json.dumps(datasets)} }},
  options: {{
    responsive: true, maintainAspectRatio: false,
    plugins: {{ legend: {{ position: 'bottom', labels: {{ font: {{ size: 11 }} }} }} }},
    scales: {{
      x: {{ title: {{ display: true, text: 'VAF Range' }} }},
      y: {{ beginAtZero: true, title: {{ display: true, text: 'Variant Count' }} }}
    }}
  }}
}});""")
        html.append('</div></div>')

    # --- Control section ---
    if controls:
        html.append('<h2>Control Sample AF Distributions</h2>')
        html.append('<div class="section"><div class="chart-section">')
        for vcf_type, vcf_label in [('truths', 'Truths'), ('candidates', 'Candidates')]:
            cid = f'chart_{chart_id}'
            chart_id += 1
            html.append(f'<div class="chart-box"><h4>{vcf_label} — AF Distribution</h4>')
            html.append(f'<canvas id="{cid}"></canvas></div>')

            datasets = []
            for i, s in enumerate(controls):
                vi = s['vcfs'].get(vcf_type, {})
                if not vi.get('exists') or not vi.get('pass_vafs'):
                    continue
                hist = bin_vafs(vi['pass_vafs'])
                ci = i % len(PALETTE)
                rgb = PALETTE[ci][0]
                datasets.append({
                    'label': s['meta']['label'],
                    'data': hist,
                    'backgroundColor': f'rgba({rgb},0.35)',
                    'borderColor': f'rgba({rgb},1)',
                    'borderWidth': 1,
                })

            chart_scripts.append(f"""
new Chart(document.getElementById('{cid}'), {{
  type: 'bar',
  data: {{ labels: {json.dumps(AF_BIN_LABELS)}, datasets: {json.dumps(datasets)} }},
  options: {{
    responsive: true, maintainAspectRatio: false,
    plugins: {{ legend: {{ position: 'bottom', labels: {{ font: {{ size: 11 }} }} }} }},
    scales: {{
      x: {{ title: {{ display: true, text: 'VAF Range' }} }},
      y: {{ beginAtZero: true, title: {{ display: true, text: 'Variant Count' }} }}
    }}
  }}
}});""")
        html.append('</div></div>')

    # --- Footer and scripts ---
    html.append('<div class="footer">Generated by generate_vaf_qc_report.py</div>')
    html.append('<script>')
    html.append('\n'.join(chart_scripts))
    html.append('</script>')
    html.append('</div></body></html>')

    with open(output_path, 'w') as f:
        f.write('\n'.join(html))


# ---------------------------------------------------------------------------
# JSON report (with metadata, for playground consumption)
# ---------------------------------------------------------------------------

def generate_json_report(samples, output_path):
    gene_groups, mixture_groups, controls = group_samples(samples)
    out = {'gene_groups': {}, 'mixtures': {}, 'controls': [], 'af_bin_labels': AF_BIN_LABELS}

    def _sample_json(s):
        m = s['meta']
        result = {
            'sample_name': m['sample_name'],
            'label': m['label'],
            'gene': m['gene'],
            'platform': m['platform'],
            'assay': m['assay'],
            'coverage': m['coverage'],
            'lab': m['lab'],
            'batch': m['batch'],
            'sample_type': m['sample_type'],
            'prep': m['prep'],
            'vcfs': {},
        }
        for vt in ('truths', 'candidates'):
            vi = s['vcfs'].get(vt, {})
            if not vi.get('exists'):
                result['vcfs'][vt] = {'exists': False}
                continue
            entry = {
                'exists': True,
                'total_variants': vi['total_variants'],
                'pass_count': vi['pass_count'],
                'lowvaf_count': vi['lowvaf_count'],
                'pass_hist': vi['pass_hist'],
            }
            if vi.get('pass_stats'):
                st = vi['pass_stats']
                entry['stats'] = {k: round(float(v), 4) for k, v in st.items()}
            result['vcfs'][vt] = entry
        return result

    for gene, group in gene_groups.items():
        out['gene_groups'][gene] = [_sample_json(s) for s in group]
    for gene, group in mixture_groups.items():
        out['mixtures'][gene] = [_sample_json(s) for s in group]
    out['controls'] = [_sample_json(s) for s in controls]

    with open(output_path, 'w') as f:
        json.dump(out, f, indent=2)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    results_dir = '/wrk/mdic_repun/results'
    output_dir = Path(results_dir) / 'qc_reports'
    output_dir.mkdir(exist_ok=True)

    print("Processing samples...")
    samples = process_all_samples(results_dir)
    print(f"Found {len(samples)} samples with Repun output")

    gene_groups, mixture_groups, controls = group_samples(samples)
    print(f"  {len(gene_groups)} gene groups, {sum(len(v) for v in mixture_groups.values())} mixture samples, {len(controls)} controls")

    html_out = output_dir / 'vaf_qc_report.html'
    json_out = output_dir / 'vaf_qc_report.json'

    print(f"Generating HTML report: {html_out}")
    generate_html_report(samples, str(html_out))

    print(f"Generating JSON report: {json_out}")
    generate_json_report(samples, str(json_out))

    print(f"\nReports generated:")
    print(f"  HTML: {html_out}")
    print(f"  JSON: {json_out}")

    # Console summary
    print(f"\n{'='*72}")
    print("SUMMARY (PASS variants, Truths VCF)")
    print(f"{'='*72}")
    for gene, group in gene_groups.items():
        print(f"\n  {gene}:")
        for s in group:
            m = s['meta']
            t = s['vcfs'].get('truths', {})
            ps = t.get('pass_stats')
            if ps:
                print(f"    {m['label']:40s}  PASS: {ps['count']:>6d}  Mean AF: {ps['mean']:5.1f}%  Median: {ps['median']:5.1f}%")
            else:
                print(f"    {m['label']:40s}  No PASS variants")


if __name__ == '__main__':
    main()
