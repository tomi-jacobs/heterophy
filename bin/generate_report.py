#!/usr/bin/env python3
"""
generate_report.py
─────────────────────────────────────────────────────────────────────────────
HeteroPhy Pipeline — Step 12: Interactive HTML Report Generator

Generates a publication-quality, self-contained HTML report summarising
all pipeline results including:
  - Per-sample coverage and masking statistics
  - SNP calling summary (Ts/Tv ratios, heterozygosity rates)
  - Phylogenetic tree comparison (RF distances, bootstrap distributions)
  - Strategy verdict and recommendation
─────────────────────────────────────────────────────────────────────────────
"""

import argparse
import json
import os
import sys
from pathlib import Path
from datetime import datetime


REPORT_TEMPLATE = '''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>{{ title }}</title>
<style>
  @import url('https://fonts.googleapis.com/css2?family=Space+Mono:wght@400;700&family=DM+Sans:wght@300;400;500;600&display=swap');

  :root {
    --bg:       #0a0c14;
    --surface:  #111422;
    --card:     #161929;
    --border:   #1f2440;
    --accent1:  #00e5b0;
    --accent2:  #7c6af7;
    --accent3:  #ff6b6b;
    --accent4:  #ffd166;
    --text:     #e4e6f0;
    --muted:    #6b7194;
    --mono:     'Space Mono', monospace;
    --sans:     'DM Sans', sans-serif;
  }

  * { box-sizing: border-box; margin: 0; padding: 0; }

  body {
    background: var(--bg);
    color: var(--text);
    font-family: var(--sans);
    font-size: 15px;
    line-height: 1.65;
  }

  /* ── Header ──────────────────────────────────────────────────────── */
  .hero {
    background: linear-gradient(135deg, #0f1230 0%, #0a0c14 60%);
    border-bottom: 1px solid var(--border);
    padding: 52px 64px 44px;
    position: relative;
    overflow: hidden;
  }
  .hero::before {
    content: '';
    position: absolute;
    top: -120px; right: -120px;
    width: 400px; height: 400px;
    background: radial-gradient(circle, rgba(0,229,176,0.07) 0%, transparent 70%);
    border-radius: 50%;
  }
  .hero::after {
    content: '';
    position: absolute;
    bottom: -80px; left: 30%;
    width: 300px; height: 300px;
    background: radial-gradient(circle, rgba(124,106,247,0.05) 0%, transparent 70%);
    border-radius: 50%;
  }
  .hero-tag {
    font-family: var(--mono);
    font-size: 11px;
    letter-spacing: 3px;
    color: var(--accent1);
    text-transform: uppercase;
    margin-bottom: 12px;
    display: block;
  }
  .hero h1 {
    font-size: 42px;
    font-weight: 600;
    letter-spacing: -1px;
    background: linear-gradient(135deg, #ffffff 0%, var(--accent1) 100%);
    -webkit-background-clip: text;
    -webkit-text-fill-color: transparent;
    background-clip: text;
    margin-bottom: 8px;
    line-height: 1.15;
  }
  .hero-sub {
    color: var(--muted);
    font-size: 14px;
    margin-top: 6px;
  }
  .hero-meta {
    margin-top: 24px;
    display: flex;
    gap: 28px;
    flex-wrap: wrap;
  }
  .meta-pill {
    background: rgba(255,255,255,0.04);
    border: 1px solid var(--border);
    border-radius: 6px;
    padding: 6px 14px;
    font-family: var(--mono);
    font-size: 11px;
    color: var(--muted);
  }
  .meta-pill span { color: var(--accent1); }

  /* ── Navigation ──────────────────────────────────────────────────── */
  .nav {
    background: var(--surface);
    border-bottom: 1px solid var(--border);
    padding: 0 64px;
    display: flex;
    gap: 0;
    position: sticky;
    top: 0;
    z-index: 100;
  }
  .nav a {
    color: var(--muted);
    text-decoration: none;
    font-size: 13px;
    font-weight: 500;
    padding: 14px 18px;
    border-bottom: 2px solid transparent;
    transition: all 0.2s;
  }
  .nav a:hover { color: var(--text); border-color: var(--accent2); }

  /* ── Layout ──────────────────────────────────────────────────────── */
  .container { max-width: 1280px; margin: 0 auto; padding: 0 64px; }
  .section { padding: 52px 0 36px; border-bottom: 1px solid var(--border); }
  .section:last-child { border-bottom: none; }

  .section-header {
    display: flex;
    align-items: baseline;
    gap: 16px;
    margin-bottom: 28px;
  }
  .section-num {
    font-family: var(--mono);
    font-size: 11px;
    color: var(--accent2);
    letter-spacing: 2px;
  }
  .section-title {
    font-size: 22px;
    font-weight: 600;
    letter-spacing: -0.3px;
  }

  /* ── Cards ───────────────────────────────────────────────────────── */
  .card-grid {
    display: grid;
    grid-template-columns: repeat(auto-fill, minmax(220px, 1fr));
    gap: 16px;
    margin-bottom: 28px;
  }
  .card {
    background: var(--card);
    border: 1px solid var(--border);
    border-radius: 10px;
    padding: 22px 20px;
    position: relative;
    overflow: hidden;
  }
  .card::before {
    content: '';
    position: absolute;
    top: 0; left: 0; right: 0;
    height: 2px;
    background: var(--accent1);
  }
  .card.accent2::before { background: var(--accent2); }
  .card.accent3::before { background: var(--accent3); }
  .card.accent4::before { background: var(--accent4); }

  .card-label {
    font-size: 11px;
    text-transform: uppercase;
    letter-spacing: 1.5px;
    color: var(--muted);
    margin-bottom: 8px;
    font-family: var(--mono);
  }
  .card-value {
    font-size: 32px;
    font-weight: 600;
    font-family: var(--mono);
    line-height: 1;
    color: var(--accent1);
  }
  .card.accent2 .card-value { color: var(--accent2); }
  .card.accent3 .card-value { color: var(--accent3); }
  .card.accent4 .card-value { color: var(--accent4); }
  .card-sub { font-size: 12px; color: var(--muted); margin-top: 6px; }

  /* ── Verdict box ─────────────────────────────────────────────────── */
  .verdict {
    background: var(--card);
    border: 1px solid var(--accent1);
    border-left: 4px solid var(--accent1);
    border-radius: 10px;
    padding: 24px 28px;
    margin-bottom: 28px;
  }
  .verdict.warn {
    border-color: var(--accent4);
    border-left-color: var(--accent4);
  }
  .verdict.alert {
    border-color: var(--accent3);
    border-left-color: var(--accent3);
  }
  .verdict-title {
    font-family: var(--mono);
    font-size: 11px;
    letter-spacing: 2px;
    text-transform: uppercase;
    color: var(--accent1);
    margin-bottom: 8px;
  }
  .verdict.warn .verdict-title { color: var(--accent4); }
  .verdict.alert .verdict-title { color: var(--accent3); }
  .verdict-text { font-size: 15px; line-height: 1.6; }

  /* ── Table ───────────────────────────────────────────────────────── */
  .table-wrap { overflow-x: auto; margin: 0 0 24px; border-radius: 10px;
                border: 1px solid var(--border); }
  table { width: 100%; border-collapse: collapse; font-size: 13px; }
  thead { background: var(--surface); }
  th {
    text-align: left;
    padding: 12px 16px;
    font-family: var(--mono);
    font-size: 10px;
    letter-spacing: 1.5px;
    text-transform: uppercase;
    color: var(--muted);
    border-bottom: 1px solid var(--border);
  }
  td {
    padding: 11px 16px;
    border-bottom: 1px solid var(--border);
    font-family: var(--mono);
    font-size: 12px;
  }
  tr:last-child td { border-bottom: none; }
  tr:hover td { background: rgba(255,255,255,0.02); }
  .badge {
    display: inline-block;
    padding: 2px 9px;
    border-radius: 4px;
    font-size: 10px;
    font-family: var(--mono);
    font-weight: 700;
  }
  .badge-green  { background: rgba(0,229,176,0.12); color: var(--accent1); }
  .badge-purple { background: rgba(124,106,247,0.12); color: var(--accent2); }
  .badge-red    { background: rgba(255,107,107,0.12); color: var(--accent3); }
  .badge-yellow { background: rgba(255,209,102,0.12); color: var(--accent4); }

  /* ── RF bar ──────────────────────────────────────────────────────── */
  .rf-bar-wrap { margin: 16px 0 24px; }
  .rf-bar-track {
    background: var(--border);
    border-radius: 4px;
    height: 10px;
    width: 100%;
    position: relative;
    overflow: hidden;
  }
  .rf-bar-fill {
    height: 100%;
    border-radius: 4px;
    transition: width 1s ease;
  }
  .rf-label {
    display: flex;
    justify-content: space-between;
    font-family: var(--mono);
    font-size: 11px;
    color: var(--muted);
    margin-bottom: 6px;
  }
  .rf-value {
    font-family: var(--mono);
    font-size: 36px;
    font-weight: 700;
    margin-bottom: 4px;
  }

  /* ── Two column ──────────────────────────────────────────────────── */
  .two-col { display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }
  @media (max-width: 768px) { .two-col { grid-template-columns: 1fr; } }

  /* ── Footer ──────────────────────────────────────────────────────── */
  footer {
    background: var(--surface);
    border-top: 1px solid var(--border);
    padding: 28px 64px;
    display: flex;
    justify-content: space-between;
    align-items: center;
    flex-wrap: wrap;
    gap: 16px;
  }
  footer .logo {
    font-family: var(--mono);
    font-size: 13px;
    color: var(--accent1);
  }
  footer .credits {
    font-size: 12px;
    color: var(--muted);
  }

  .progress-section { margin: 8px 0 20px; }
  .step-row {
    display: flex;
    align-items: center;
    gap: 12px;
    padding: 10px 0;
    border-bottom: 1px solid var(--border);
    font-size: 13px;
  }
  .step-row:last-child { border: none; }
  .step-dot {
    width: 8px; height: 8px;
    border-radius: 50%;
    background: var(--accent1);
    flex-shrink: 0;
  }
  .step-dot.skip { background: var(--muted); }
  .step-name { flex: 1; }
  .step-status { font-family: var(--mono); font-size: 11px; }
</style>
</head>
<body>

<!-- HERO -->
<div class="hero">
  <span class="hero-tag">Bioinformatics Pipeline Report</span>
  <h1>{{ title }}</h1>
  <p class="hero-sub">SNP &amp; Heterozygosity-Aware Phylogenetic Inference — Aim 2</p>
  <div class="hero-meta">
    <div class="meta-pill">Version <span>{{ pipeline_ver }}</span></div>
    <div class="meta-pill">Generated <span>{{ generated }}</span></div>
    <div class="meta-pill">Samples <span>{{ n_samples }}</span></div>
    <div class="meta-pill">Strategies <span>Collapsed &amp; Retained</span></div>
  </div>
</div>

<!-- NAV -->
<nav class="nav">
  <a href="#overview">Overview</a>
  <a href="#coverage">Coverage</a>
  <a href="#snps">SNPs</a>
  <a href="#phylogeny">Phylogeny</a>
  <a href="#verdict">Verdict</a>
  <a href="#methods">Methods</a>
</nav>

<div class="container">

  <!-- OVERVIEW -->
  <section class="section" id="overview">
    <div class="section-header">
      <span class="section-num">01 —</span>
      <h2 class="section-title">Pipeline Overview</h2>
    </div>
    <div class="card-grid">
      <div class="card">
        <div class="card-label">Samples Processed</div>
        <div class="card-value">{{ n_samples }}</div>
        <div class="card-sub">De novo transcriptomes</div>
      </div>
      <div class="card accent2">
        <div class="card-label">Strategies Compared</div>
        <div class="card-value">2</div>
        <div class="card-sub">Collapsed vs Retained</div>
      </div>
      <div class="card accent4">
        <div class="card-label">Total SNPs Called</div>
        <div class="card-value">{{ total_snps }}</div>
        <div class="card-sub">Across all samples</div>
      </div>
      <div class="card accent3">
        <div class="card-label">Heterozygous Sites</div>
        <div class="card-value">{{ total_het }}</div>
        <div class="card-sub">IUPAC-encoded in retained</div>
      </div>
    </div>

    <div class="progress-section">
      <div class="step-row">
        <div class="step-dot"></div>
        <div class="step-name">HISAT2 — Splice-aware read alignment to transcriptome</div>
        <div class="step-status badge badge-green">COMPLETE</div>
      </div>
      <div class="step-row">
        <div class="step-dot"></div>
        <div class="step-name">Coverage Assessment — Per-nucleotide depth (Samtools depth)</div>
        <div class="step-status badge badge-green">COMPLETE</div>
      </div>
      <div class="step-row">
        <div class="step-dot"></div>
        <div class="step-name">Transcriptome Masking — Low-coverage positions → N</div>
        <div class="step-status badge badge-green">COMPLETE</div>
      </div>
      <div class="step-row">
        <div class="step-dot"></div>
        <div class="step-name">SNP Calling — BCFtools mpileup + call, quality filtered</div>
        <div class="step-status badge badge-green">COMPLETE</div>
      </div>
      <div class="step-row">
        <div class="step-dot"></div>
        <div class="step-name">Strategy A — Collapse heterozygous sites (reference allele)</div>
        <div class="step-status badge badge-green">COMPLETE</div>
      </div>
      <div class="step-row">
        <div class="step-dot"></div>
        <div class="step-name">Strategy B — Retain heterozygous sites (IUPAC ambiguity codes)</div>
        <div class="step-status badge badge-green">COMPLETE</div>
      </div>
      <div class="step-row">
        <div class="step-dot"></div>
        <div class="step-name">MAFFT — Multiple sequence alignment per strategy</div>
        <div class="step-status badge badge-green">COMPLETE</div>
      </div>
      <div class="step-row">
        <div class="step-dot"></div>
        <div class="step-name">IQ-TREE2 — Maximum likelihood phylogenetic inference (UFBoot)</div>
        <div class="step-status badge badge-green">COMPLETE</div>
      </div>
      <div class="step-row">
        <div class="step-dot"></div>
        <div class="step-name">Tree Comparison — Robinson-Foulds distance, KS bootstrap test</div>
        <div class="step-status badge badge-green">COMPLETE</div>
      </div>
    </div>
  </section>

  <!-- COVERAGE -->
  <section class="section" id="coverage">
    <div class="section-header">
      <span class="section-num">02 —</span>
      <h2 class="section-title">Coverage &amp; Masking</h2>
    </div>
    <p style="color:var(--muted); margin-bottom:24px; max-width:720px;">
      Nucleotide positions with fewer than <code style="color:var(--accent1)">{{ min_depth }}×</code>
      read depth were masked to <code style="color:var(--accent1)">N</code>.
      In phylogenetic likelihood models, ambiguous states are treated as all states simultaneously,
      avoiding directional bias — safer than retaining potentially incorrect bases.
    </p>
    <div class="table-wrap">
      <table>
        <thead>
          <tr>
            <th>Sample</th>
            <th>Total Bases</th>
            <th>Masked Bases</th>
            <th>% Masked</th>
            <th>Contigs Retained</th>
            <th>Contigs Filtered</th>
            <th>Coverage Status</th>
          </tr>
        </thead>
        <tbody>
          {{ coverage_rows }}
        </tbody>
      </table>
    </div>
  </section>

  <!-- SNPs -->
  <section class="section" id="snps">
    <div class="section-header">
      <span class="section-num">03 —</span>
      <h2 class="section-title">SNP Summary</h2>
    </div>
    <p style="color:var(--muted); margin-bottom:24px; max-width:720px;">
      SNPs were called using BCFtools mpileup/call on masked transcriptomes.
      High-confidence SNPs were retained (QUAL≥30, DP≥{{ min_snp_depth }}).
      Ts/Tv ratios near 2.0 indicate well-calibrated SNP calls; values far outside
      this range may suggest systematic sequencing or assembly artefacts.
    </p>
    <div class="table-wrap">
      <table>
        <thead>
          <tr>
            <th>Sample</th>
            <th>Total SNPs</th>
            <th>Heterozygous</th>
            <th>Hom. Alt</th>
            <th>Ts/Tv Ratio</th>
            <th>Contigs w/ SNPs</th>
            <th>Het. Rate</th>
          </tr>
        </thead>
        <tbody>
          {{ snp_rows }}
        </tbody>
      </table>
    </div>
  </section>

  <!-- PHYLOGENY -->
  <section class="section" id="phylogeny">
    <div class="section-header">
      <span class="section-num">04 —</span>
      <h2 class="section-title">Phylogenetic Comparison</h2>
    </div>

    <div class="two-col" style="margin-bottom:28px;">
      <div>
        <p style="color:var(--muted); font-size:13px; margin-bottom:16px;">
          Robinson-Foulds (RF) distance quantifies topological differences between trees.
          An RF of 0 = identical topologies; 1.0 = maximally different.
        </p>
        <div class="rf-label">
          <span>Normalized RF Distance</span>
          <span>Collapsed vs Retained</span>
        </div>
        <div class="rf-value" style="color:{{ rf_color }}">{{ rf_norm }}</div>
        <div class="rf-bar-wrap">
          <div class="rf-bar-track">
            <div class="rf-bar-fill"
                 style="width:{{ rf_pct }}%; background:{{ rf_color }};"></div>
          </div>
        </div>
        <p style="color:var(--muted); font-size:12px; font-family:var(--mono);">
          Raw RF: {{ rf_raw }} / {{ rf_max }} partitions
        </p>
      </div>

      <div>
        <div class="card-grid" style="grid-template-columns:1fr 1fr;">
          <div class="card">
            <div class="card-label">Mean Bootstrap<br>Collapsed</div>
            <div class="card-value" style="font-size:26px;">{{ mean_bs_collapsed }}%</div>
          </div>
          <div class="card accent2">
            <div class="card-label">Mean Bootstrap<br>Retained</div>
            <div class="card-value" style="font-size:26px; color:var(--accent2);">{{ mean_bs_retained }}%</div>
          </div>
          <div class="card accent4">
            <div class="card-label">Shared Clades</div>
            <div class="card-value" style="font-size:26px; color:var(--accent4);">{{ shared_clades }}</div>
          </div>
          <div class="card accent3">
            <div class="card-label">KS Test p-value</div>
            <div class="card-value" style="font-size:22px; color:{{ ks_color }};">{{ ks_pvalue }}</div>
          </div>
        </div>
      </div>
    </div>
  </section>

  <!-- VERDICT -->
  <section class="section" id="verdict">
    <div class="section-header">
      <span class="section-num">05 —</span>
      <h2 class="section-title">Verdict &amp; Recommendation</h2>
    </div>
    <div class="verdict {{ verdict_class }}">
      <div class="verdict-title">Topology Verdict</div>
      <div class="verdict-text">{{ topology_verdict }}</div>
    </div>
    <div class="verdict {{ support_class }}" style="margin-top:16px;">
      <div class="verdict-title">Bootstrap Support</div>
      <div class="verdict-text">{{ bs_verdict }}</div>
    </div>
    <div class="verdict" style="margin-top:16px; border-color:var(--accent2); border-left-color:var(--accent2);">
      <div class="verdict-title" style="color:var(--accent2);">Recommendation</div>
      <div class="verdict-text">{{ recommendation }}</div>
    </div>
  </section>

  <!-- METHODS -->
  <section class="section" id="methods">
    <div class="section-header">
      <span class="section-num">06 —</span>
      <h2 class="section-title">Methods Summary</h2>
    </div>
    <p style="color:var(--muted); max-width:860px; line-height:1.8; font-size:14px;">
      Raw reads were aligned to assembled transcriptomes using
      <strong style="color:var(--text)">HISAT2 v2.2.1</strong>
      (splice-aware aligner; exon-exon junction aware for eukaryotic transcriptomes).
      Per-nucleotide coverage was assessed using <strong style="color:var(--text)">Samtools depth</strong>.
      Positions with fewer than {{ min_depth }}× coverage were masked to <em>N</em>.
      SNPs were called on masked transcriptomes using
      <strong style="color:var(--text)">BCFtools mpileup/call</strong>
      (QUAL≥30, DP≥{{ min_snp_depth }}, MQ≥20) and filtered using custom Python scripts.
      Two strategies were applied: (1) <strong style="color:var(--accent1)">Collapsed</strong> —
      heterozygous sites resolved to the reference allele;
      (2) <strong style="color:var(--accent2)">Retained</strong> —
      heterozygous sites encoded with IUPAC ambiguity codes.
      Multiple sequence alignment was performed with
      <strong style="color:var(--text)">MAFFT v7.520</strong> (--auto).
      Maximum likelihood phylogenetic inference was carried out using
      <strong style="color:var(--text)">IQ-TREE2 v2.2.6</strong>
      with ModelFinder model selection and {{ bootstrap }} ultrafast bootstrap replicates.
      Topologies were compared using Robinson-Foulds distances via
      <strong style="color:var(--text)">ETE3</strong>.
      Bootstrap distributions were tested for significant difference using the
      two-sample Kolmogorov-Smirnov test (SciPy).
    </p>
    <div style="margin-top:24px; padding:20px; background:var(--card);
                border:1px solid var(--border); border-radius:10px;
                font-family:var(--mono); font-size:12px; color:var(--muted);">
      <strong style="color:var(--accent1);">Citation:</strong><br>
      [Your Name] et al. ({{ year }}). HeteroPhy: A bioinformatics pipeline for
      investigating the role of SNPs and heterozygosity in resolving evolutionary
      relationships. [Journal]. doi: [pending]
    </div>
  </section>

</div>

<footer>
  <div class="logo">HeteroPhy v{{ pipeline_ver }}</div>
  <div class="credits">
    Generated {{ generated }} &nbsp;·&nbsp;
    Nextflow {{ pipeline_ver }} + Python &nbsp;·&nbsp;
    HISAT2 · Samtools · BCFtools · MAFFT · IQ-TREE2 · ETE3
  </div>
</footer>

</body>
</html>'''


def parse_args():
    parser = argparse.ArgumentParser(description="Generate HeteroPhy HTML report")
    parser.add_argument("--outdir",       default=".")
    parser.add_argument("--title",        default="HeteroPhy Analysis Report")
    parser.add_argument("--pipeline_ver", default="1.0.0")
    return parser.parse_args()


def load_json_files(pattern_dir, suffix):
    """Load all JSON files matching a suffix in a directory."""
    data = {}
    for f in Path(pattern_dir).glob(f"**/*{suffix}"):
        try:
            with open(f) as fh:
                d = json.load(fh)
            sample = d.get("sample", f.stem.replace(suffix.replace(".json",""), ""))
            data[sample] = d
        except Exception:
            pass
    return data


def format_number(n):
    if n is None:
        return "—"
    try:
        return f"{int(n):,}"
    except Exception:
        return str(n)


def build_coverage_rows(coverage_data):
    rows = []
    for sample, d in sorted(coverage_data.items()):
        frac = d.get("fraction_low_coverage", 0)
        pct  = f"{frac*100:.1f}%"
        status = (
            '<span class="badge badge-green">GOOD</span>'
            if frac < 0.2 else
            '<span class="badge badge-yellow">MODERATE</span>'
            if frac < 0.5 else
            '<span class="badge badge-red">HIGH MASKING</span>'
        )
        rows.append(
            f"<tr>"
            f"<td>{sample}</td>"
            f"<td>{format_number(d.get('total_positions'))}</td>"
            f"<td>{format_number(d.get('low_coverage_positions'))}</td>"
            f"<td>{pct}</td>"
            f"<td>—</td><td>—</td>"
            f"<td>{status}</td>"
            f"</tr>"
        )
    return "\n".join(rows) if rows else "<tr><td colspan='7'>No data</td></tr>"


def build_snp_rows(snp_data):
    rows = []
    for sample, d in sorted(snp_data.items()):
        snps     = d.get("snps", 0)
        het      = d.get("heterozygous", 0)
        hom_alt  = d.get("homozygous_alt", 0)
        tstv     = d.get("ts_tv_ratio")
        n_ctg    = d.get("n_contigs_with_snps", 0)
        het_rate = f"{het/snps*100:.1f}%" if snps > 0 else "—"

        tstv_badge = ""
        if tstv is not None:
            if 1.5 <= tstv <= 3.0:
                tstv_badge = f' <span class="badge badge-green">OK</span>'
            else:
                tstv_badge = f' <span class="badge badge-yellow">CHECK</span>'

        rows.append(
            f"<tr>"
            f"<td>{sample}</td>"
            f"<td>{format_number(snps)}</td>"
            f"<td>{format_number(het)}</td>"
            f"<td>{format_number(hom_alt)}</td>"
            f"<td>{tstv if tstv else '—'}{tstv_badge}</td>"
            f"<td>{format_number(n_ctg)}</td>"
            f"<td>{het_rate}</td>"
            f"</tr>"
        )
    return "\n".join(rows) if rows else "<tr><td colspan='7'>No data</td></tr>"


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    now = datetime.now().strftime("%Y-%m-%d %H:%M UTC")
    year = datetime.now().year

    # ── Load data ─────────────────────────────────────────────────────────────
    coverage_data = load_json_files(".", "coverage_stats.json")
    snp_data      = load_json_files(".", "snp_stats.json")

    comparison = {}
    for f in Path(".").glob("**/comparison_results.json"):
        try:
            with open(f) as fh:
                comparison = json.load(fh)
            break
        except Exception:
            pass

    # ── Extract comparison metrics ────────────────────────────────────────────
    rf_results  = comparison.get("topology", {}).get("robinson_foulds", {})
    bs_comp     = comparison.get("bootstrap_comparison", {})
    clade_res   = comparison.get("clade_concordance", {})
    summary     = comparison.get("summary", {})
    topology_verdict = comparison.get("topology", {}).get("verdict",
                       "Analysis pending — run tree comparison step.")

    rf_norm = rf_results.get("normalized_rf", 0.0) or 0.0
    rf_raw  = rf_results.get("rf", "—")
    rf_max  = rf_results.get("max_rf", "—")
    rf_pct  = min(100, int(rf_norm * 100))
    rf_color = "#00e5b0" if rf_norm < 0.1 else "#ffd166" if rf_norm < 0.3 else "#ff6b6b"

    mean_bs_c = bs_comp.get("mean_bs_1", "—") or "—"
    mean_bs_r = bs_comp.get("mean_bs_2", "—") or "—"
    shared_cl = clade_res.get("shared_clades", "—")

    ks_pvalue = bs_comp.get("ks_pvalue")
    ks_color = "#ff6b6b" if (ks_pvalue is not None and ks_pvalue < 0.05) else "#00e5b0"
    ks_pvalue_str = f"{ks_pvalue:.4f}" if ks_pvalue is not None else "—"

    ks_sig = bs_comp.get("ks_significant", False)
    bs_verdict = (
        bs_comp.get("ks_interpretation",
                    "Kolmogorov-Smirnov test: bootstrap support distributions were "
                    "compared between strategies.")
    )

    verdict_class = (
        "" if rf_norm < 0.05 else
        "warn" if rf_norm < 0.3 else "alert"
    )
    support_class = "alert" if ks_sig else ""

    recommendation = summary.get("recommendation",
        "Run the full pipeline to generate recommendations.")

    # ── Aggregate totals ──────────────────────────────────────────────────────
    total_snps = sum(d.get("snps", 0) for d in snp_data.values())
    total_het  = sum(d.get("heterozygous", 0) for d in snp_data.values())
    n_samples  = max(len(snp_data), len(coverage_data), 1)

    # ── Fill template ─────────────────────────────────────────────────────────
    html = REPORT_TEMPLATE
    replacements = {
        "{{ title }}":            args.title,
        "{{ pipeline_ver }}":     args.pipeline_ver,
        "{{ generated }}":        now,
        "{{ year }}":             str(year),
        "{{ n_samples }}":        str(n_samples),
        "{{ total_snps }}":       format_number(total_snps),
        "{{ total_het }}":        format_number(total_het),
        "{{ min_depth }}":        "5",
        "{{ min_snp_depth }}":    "10",
        "{{ bootstrap }}":        "1000",
        "{{ coverage_rows }}":    build_coverage_rows(coverage_data),
        "{{ snp_rows }}":         build_snp_rows(snp_data),
        "{{ rf_norm }}":          f"{rf_norm:.3f}",
        "{{ rf_raw }}":           str(rf_raw),
        "{{ rf_max }}":           str(rf_max),
        "{{ rf_pct }}":           str(rf_pct),
        "{{ rf_color }}":         rf_color,
        "{{ mean_bs_collapsed }}": str(mean_bs_c),
        "{{ mean_bs_retained }}": str(mean_bs_r),
        "{{ shared_clades }}":    str(shared_cl),
        "{{ ks_pvalue }}":        ks_pvalue_str,
        "{{ ks_color }}":         ks_color,
        "{{ topology_verdict }}": topology_verdict,
        "{{ bs_verdict }}":       bs_verdict,
        "{{ verdict_class }}":    verdict_class,
        "{{ support_class }}":    support_class,
        "{{ recommendation }}":   recommendation,
    }

    for k, v in replacements.items():
        html = html.replace(k, v)

    out_path = outdir / "heterophy_report.html"
    with open(out_path, "w") as f:
        f.write(html)

    print(f"[HeteroPhy] Report written to: {out_path}")


if __name__ == "__main__":
    main()
