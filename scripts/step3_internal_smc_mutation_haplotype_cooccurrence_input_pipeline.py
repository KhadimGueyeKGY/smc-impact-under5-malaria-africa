#!/usr/bin/env python3

import argparse
import csv
import json
import os
import re
import shlex
import time
from pathlib import Path

import pandas as pd
import pysam
from rich.console import Console
from rich.panel import Panel
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TimeElapsedColumn, TimeRemainingColumn
from rich.table import Table

console = Console()

AA1_TO_AA3 = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
    "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His", "I": "Ile",
    "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
    "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
    "*": "Ter", "X": "Xaa"
}

AA3_TO_AA1 = {v: k for k, v in AA1_TO_AA3.items()}

DEFAULT_VCF_PATTERNS = [
    "vcf_good_qual/{barcode}.ensemble.good.vcf.gz",
    "vcf_annotated/{barcode}.ensemble.anno.corrected.vcf.gz",
    "vcf_annotated/{barcode}.ensemble.anno.vcf.gz",
]

PREFERRED_FOLDERS = {
    "Jigawa": ["amplicon_Jigawa", "Jigawa_plate_1_panel_2"],
    "Katsina": ["Katsina_final", "Katsina_plate_1_panel_1"],
}

GENE_FIX = {
    "pfdhfr": "dhfr",
    "dhfr": "dhfr",
    "pfdhps": "dhps",
    "dhps": "dhps",
    "pfmdr1": "mdr1",
    "mdr1": "mdr1",
    "pfcrt": "crt",
    "crt": "crt",
    "pfk13": "k13",
    "k13": "k13",
    "mdr2": "mdr2",
    "ferredoxin": "ferredoxin",
    "exonuclease": "exonuclease",
    "coronin": "coronin",
    "msp1": "msp1",
    "msp2": "msp2",
    "polya": "polyA",
    "polya1": "polyA",
    "polyA": "polyA",
    "ta1": "ta1",
}


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--base-dir", required=True)
    p.add_argument("--metadata-tsv", required=True)
    p.add_argument("--fragment-manifest", required=True)
    p.add_argument("--bed-file", required=True)
    p.add_argument("--haplotype-file", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--min-reads", type=int, default=100)
    p.add_argument("--vcf-patterns", nargs="+", default=DEFAULT_VCF_PATTERNS)
    p.add_argument("--bam-pattern", default="bam/{barcode}.bam")
    p.add_argument("--state-filter", nargs="*", default=None)
    return p.parse_args()


def ensure_dir(path):
    Path(path).mkdir(parents=True, exist_ok=True)


def print_header():
    console.print(
        Panel.fit(
            "[bold cyan]STEP 3 - INTERNAL SMC CO-OCCURRENCE INPUT PIPELINE[/bold cyan]\n"
            "[bold white]High-quality callable samples, global mutation matrix, global haplotype matrix[/bold white]",
            border_style="bright_blue",
        )
    )


def print_step(title):
    console.print(Panel.fit(f"[bold white]{title}[/bold white]", border_style="green"))


def print_success(msg):
    console.print(f"[bold green]✔ {msg}[/bold green]")


def print_info(msg):
    console.print(f"[bold cyan]{msg}[/bold cyan]")


def norm_col(s):
    return str(s).strip().lower().replace(" ", "_")


def clean_string_series(s):
    s = s.astype(str).str.strip()
    s = s.replace({"nan": pd.NA, "NaN": pd.NA, "None": pd.NA, "": pd.NA, ".": pd.NA})
    return s


def normalize_gene_name(x):
    if x is None or pd.isna(x):
        return None
    x = str(x).strip()
    return GENE_FIX.get(x.lower(), x)


def mut_1to3(m):
    m = str(m).strip()
    g = re.match(r"^([A-Z\*])(\d+)([A-Z\*])$", m)
    if not g:
        return None
    return f"{AA1_TO_AA3[g.group(1)]}{g.group(2)}{AA1_TO_AA3[g.group(3)]}"


def load_metadata(path, state_filter=None):
    meta = pd.read_csv(path, sep="\t", dtype=str, engine="python", quoting=csv.QUOTE_NONE, on_bad_lines="skip")
    meta.columns = [norm_col(c) for c in meta.columns]
    required = {"sample", "folder", "barcode", "sample_id", "state"}
    if not required.issubset(set(meta.columns)):
        raise ValueError("metadata missing required columns")
    meta["sample"] = clean_string_series(meta["sample"])
    meta["folder"] = clean_string_series(meta["folder"])
    meta["barcode"] = clean_string_series(meta["barcode"])
    meta["sample_id"] = clean_string_series(meta["sample_id"])
    meta["state"] = clean_string_series(meta["state"])
    meta = meta.dropna(subset=["sample", "folder", "barcode", "sample_id", "state"]).copy()
    meta["year"] = 2021
    meta["age_group"] = "<5"
    if state_filter:
        meta = meta[meta["state"].isin(set(state_filter))].copy()
    return meta


def choose_one_record_per_sample(st):
    st = st.copy()
    pref = PREFERRED_FOLDERS.get(str(st["state"].iloc[0]) if len(st) else "", [])
    if pref:
        st["_rank"] = st["folder"].apply(lambda x: pref.index(x) if x in pref else len(pref))
    else:
        st["_rank"] = 0
    st = st.sort_values(["sample_id", "_rank", "folder", "barcode"]).drop_duplicates(subset=["sample_id"], keep="first")
    st = st.drop(columns=["_rank"], errors="ignore")
    return st


def load_fragment_manifest(path):
    df = pd.read_csv(path, sep=None, engine="python", dtype=str)
    df.columns = [norm_col(c) for c in df.columns]
    rename = {}
    if "chr" in df.columns and "chromosome" not in df.columns:
        rename["chr"] = "chromosome"
    if "reference" in df.columns and "reference_id" not in df.columns:
        rename["reference"] = "reference_id"
    if "gene_symbol" in df.columns and "reference_id" not in df.columns:
        rename["gene_symbol"] = "reference_id"
    if "description" not in df.columns and "marker_type" in df.columns:
        rename["marker_type"] = "description"
    df = df.rename(columns=rename)

    if {"chromosome", "start", "end", "reference_id"}.issubset(df.columns):
        df["start"] = pd.to_numeric(df["start"], errors="coerce").astype("Int64")
        df["end"] = pd.to_numeric(df["end"], errors="coerce").astype("Int64")
    elif {"chromosome", "start_position", "end_position", "reference_id"}.issubset(df.columns):
        df["start"] = pd.to_numeric(df["start_position"], errors="coerce").astype("Int64")
        df["end"] = pd.to_numeric(df["end_position"], errors="coerce").astype("Int64")
    else:
        raise ValueError("fragment manifest missing required columns")

    if "gene_id" not in df.columns:
        df["gene_id"] = None
    if "description" not in df.columns:
        df["description"] = None

    df["chromosome"] = clean_string_series(df["chromosome"])
    df["reference_id"] = df["reference_id"].apply(normalize_gene_name)
    df["gene_id"] = clean_string_series(df["gene_id"].astype(str)).astype(object)
    df["description"] = clean_string_series(df["description"].astype(str)).astype(object)

    df = df.dropna(subset=["chromosome", "start", "end", "reference_id"]).copy()
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df["fragment_id"] = (
        df["reference_id"].astype(str)
        + "__"
        + df["chromosome"].astype(str)
        + "__"
        + df["start"].astype(str)
        + "__"
        + df["end"].astype(str)
    )
    return df[["fragment_id", "chromosome", "start", "end", "reference_id", "gene_id", "description"]].drop_duplicates()


def load_bed(path):
    rows = []
    with open(path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.lower().startswith("#chrom"):
                continue
            parts = re.split(r"\s+", line)
            if len(parts) < 6:
                continue
            chrom, start, end, mutation, gene, drug = parts[:6]
            rows.append([chrom, start, end, mutation, gene, drug])

    if not rows:
        raise ValueError("bed file empty or malformed")

    df = pd.DataFrame(rows, columns=["chromosome", "start", "end", "mutation", "gene", "drug"])
    df["start"] = pd.to_numeric(df["start"], errors="coerce").astype("Int64")
    df["end"] = pd.to_numeric(df["end"], errors="coerce").astype("Int64")
    df["chromosome"] = clean_string_series(df["chromosome"])
    df["mutation"] = clean_string_series(df["mutation"])
    df["gene"] = df["gene"].apply(normalize_gene_name)
    df["drug"] = clean_string_series(df["drug"])
    df = df.dropna(subset=["chromosome", "start", "end", "mutation", "gene", "drug"]).copy()
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df["mutation_aa3"] = df["mutation"].apply(mut_1to3)
    return df


def load_haplotypes(path):
    hap = pd.read_csv(path, sep="\t", engine="python", dtype=str)
    hap.columns = [norm_col(c) for c in hap.columns]
    required = {"drug", "gene", "mutations", "haplotype"}
    if not required.issubset(set(hap.columns)):
        raise ValueError("haplotype file missing required columns")
    hap["drug"] = clean_string_series(hap["drug"])
    hap["gene"] = clean_string_series(hap["gene"])
    hap["mutations"] = clean_string_series(hap["mutations"])
    hap["haplotype"] = clean_string_series(hap["haplotype"])
    hap = hap.dropna(subset=["drug", "gene", "mutations", "haplotype"]).copy()
    hap["gene_list"] = hap["gene"].apply(lambda x: [normalize_gene_name(t.strip()) for t in str(x).split("/") if t.strip()])
    hap["mutation_list"] = hap["mutations"].apply(lambda x: [t.strip() for t in str(x).split("+") if t.strip()])
    hap["mutation_list_aa3"] = hap["mutation_list"].apply(lambda xs: [mut_1to3(x) for x in xs if mut_1to3(x) is not None])
    return hap


def resolve_bam_path(base_dir, folder, barcode, bam_pattern):
    p = Path(base_dir) / str(folder) / bam_pattern.format(barcode=barcode)
    return p if p.exists() else None


def resolve_vcf_path(base_dir, folder, barcode, patterns):
    for pat in patterns:
        p = Path(base_dir) / str(folder) / pat.format(barcode=barcode)
        if p.exists():
            return p
    return None


def bam_count_reads(bam_path, chrom, start, end):
    if bam_path is None or not Path(bam_path).exists():
        return None
    try:
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            return int(bam.count(reference=str(chrom), start=int(start), end=int(end)))
    except Exception:
        return None


def map_markers_to_fragments(markers, fragments):
    out = []
    for r in markers.itertuples(index=False):
        g = fragments[
            (fragments["chromosome"] == r.chromosome)
            & (fragments["reference_id"] == r.gene)
            & (fragments["start"] <= int(r.start))
            & (fragments["end"] >= int(r.end))
        ]
        if g.empty:
            g = fragments[
                (fragments["chromosome"] == r.chromosome)
                & (fragments["reference_id"] == r.gene)
                & (fragments["start"] <= int(r.start))
                & (fragments["end"] >= int(r.start))
            ]
        if g.empty:
            fragment_id = None
        else:
            best = g.sort_values(["start", "end"]).iloc[0]
            fragment_id = best["fragment_id"]

        out.append(
            {
                "chromosome": r.chromosome,
                "marker_start": int(r.start),
                "marker_end": int(r.end),
                "mutation": r.mutation,
                "mutation_aa3": r.mutation_aa3,
                "gene": r.gene,
                "drug": r.drug,
                "fragment_id": fragment_id,
            }
        )
    return pd.DataFrame(out)


def grep_first_line(vcf_path, token):
    if vcf_path is None or not Path(vcf_path).exists():
        return ""
    q_vcf = shlex.quote(str(vcf_path))
    q_tok = shlex.quote(str(token))
    cmd = f"zcat {q_vcf} | grep -F {q_tok} | head -1"
    try:
        return os.popen(cmd).read().strip()
    except Exception:
        return ""


def summarize_state_sample_table(df):
    tbl = Table(title="Samples selected", show_lines=False)
    tbl.add_column("State", style="cyan")
    tbl.add_column("Samples", justify="right", style="green")
    for row in df.itertuples(index=False):
        tbl.add_row(str(row.state), str(row.n_samples))
    console.print(tbl)


def main():
    args = parse_args()
    start_time = time.time()

    print_header()
    ensure_dir(args.output_dir)

    print_step("STEP 1/6 Loading metadata and references")
    meta = load_metadata(args.metadata_tsv, args.state_filter)
    fragments = load_fragment_manifest(args.fragment_manifest)
    bed = load_bed(args.bed_file)
    haplotypes = load_haplotypes(args.haplotype_file)
    marker_map = map_markers_to_fragments(bed, fragments)
    counts_state = meta.groupby("state", as_index=False).agg(n_samples=("sample_id", "nunique")).sort_values("state")
    summarize_state_sample_table(counts_state)
    print_info(f"Metadata samples: {meta['sample_id'].nunique()}")
    print_info(f"Fragments loaded: {fragments['fragment_id'].nunique()}")
    print_info(f"Markers loaded: {marker_map['mutation'].nunique()}")
    print_info(f"Haplotypes loaded: {haplotypes['haplotype'].nunique()}")
    print_success("References loaded")

    print_step("STEP 2/6 Resolving sample records and file paths")
    selected_rows = []
    for state in sorted(meta["state"].dropna().unique()):
        st = meta[meta["state"] == state].copy()
        st = choose_one_record_per_sample(st)
        selected_rows.append(st)
    sample_df = pd.concat(selected_rows, ignore_index=True) if selected_rows else meta.iloc[0:0].copy()
    sample_df = sample_df.dropna(subset=["sample_id", "state", "folder", "barcode", "sample"]).copy()
    sample_df = sample_df.drop_duplicates(subset=["sample_id"], keep="first").copy()
    sample_df["bam_path"] = sample_df.apply(lambda r: resolve_bam_path(args.base_dir, r["folder"], r["barcode"], args.bam_pattern), axis=1)
    sample_df["vcf_path"] = sample_df.apply(lambda r: resolve_vcf_path(args.base_dir, r["folder"], r["barcode"], args.vcf_patterns), axis=1)
    sample_df["bam_exists"] = sample_df["bam_path"].apply(lambda x: x is not None)
    sample_df["vcf_exists"] = sample_df["vcf_path"].apply(lambda x: x is not None)

    sample_manifest = sample_df.copy()
    sample_manifest["bam_path"] = sample_manifest["bam_path"].apply(lambda x: str(x) if x is not None else "")
    sample_manifest["vcf_path"] = sample_manifest["vcf_path"].apply(lambda x: str(x) if x is not None else "")
    sample_manifest_tsv = Path(args.output_dir) / "step3_sample_manifest.tsv"
    sample_manifest.to_csv(sample_manifest_tsv, sep="\t", index=False)

    print_info(f"Samples with BAM: {int(sample_df['bam_exists'].sum())}")
    print_info(f"Samples with VCF: {int(sample_df['vcf_exists'].sum())}")
    print_success("Sample paths resolved")

    print_step("STEP 3/6 Computing fragment coverage and selecting high-quality samples")
    coverage_rows = []
    jobs = [(idx, row) for idx, row in sample_df.iterrows()]

    with Progress(
        SpinnerColumn(),
        TextColumn("[bold cyan]{task.description}"),
        BarColumn(bar_width=40),
        TextColumn("[bold white]{task.completed}/{task.total}"),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
        console=console,
    ) as progress:
        task = progress.add_task("Counting fragment reads", total=len(jobs))
        for _, srow in jobs:
            bam_path = srow["bam_path"]
            for frag in fragments.itertuples(index=False):
                n_reads = bam_count_reads(bam_path, frag.chromosome, frag.start, frag.end)
                coverage_rows.append(
                    {
                        "sample": srow["sample"],
                        "sample_id": srow["sample_id"],
                        "folder": srow["folder"],
                        "barcode": srow["barcode"],
                        "state": srow["state"],
                        "year": srow["year"],
                        "age_group": srow["age_group"],
                        "bam_exists": int(bool(srow["bam_exists"])),
                        "vcf_exists": int(bool(srow["vcf_exists"])),
                        "fragment_id": frag.fragment_id,
                        "chromosome": frag.chromosome,
                        "fragment_start": int(frag.start),
                        "fragment_end": int(frag.end),
                        "reference_id": frag.reference_id,
                        "gene_id": frag.gene_id,
                        "description": frag.description,
                        "n_reads": n_reads,
                        "fragment_pass_100reads": int(n_reads is not None and n_reads >= args.min_reads),
                    }
                )
            progress.advance(task)

    coverage_df = pd.DataFrame(coverage_rows)
    coverage_tsv = Path(args.output_dir) / "step3_fragment_coverage_by_sample.tsv"
    coverage_df.to_csv(coverage_tsv, sep="\t", index=False)

    coverage_summary_sample = (
        coverage_df.groupby(["sample_id", "sample", "state", "year", "age_group"], as_index=False)
        .agg(
            n_fragments=("fragment_id", "nunique"),
            n_callable_fragments=("fragment_pass_100reads", "sum"),
            mean_reads=("n_reads", "mean"),
            median_reads=("n_reads", "median"),
            min_reads=("n_reads", "min"),
            max_reads=("n_reads", "max"),
        )
        .sort_values(["state", "sample_id"])
    )
    coverage_summary_sample["all_fragments_callable"] = (coverage_summary_sample["n_fragments"] == coverage_summary_sample["n_callable_fragments"]).astype(int)
    coverage_summary_sample["callable_fraction"] = coverage_summary_sample["n_callable_fragments"] / coverage_summary_sample["n_fragments"]
    coverage_summary_sample_tsv = Path(args.output_dir) / "step3_fragment_coverage_summary_by_sample.tsv"
    coverage_summary_sample.to_csv(coverage_summary_sample_tsv, sep="\t", index=False)

    high_quality_samples = coverage_summary_sample[coverage_summary_sample["all_fragments_callable"] == 1]["sample_id"].astype(str).tolist()
    high_quality_df = sample_df[sample_df["sample_id"].astype(str).isin(high_quality_samples)].copy()
    high_quality_manifest_tsv = Path(args.output_dir) / "step3_high_quality_samples_all_fragments_callable.tsv"
    high_quality_df.to_csv(high_quality_manifest_tsv, sep="\t", index=False)

    print_info(f"Selected samples after all-fragment filter: {len(high_quality_samples)}")
    print_success("High-quality samples selected")

    print_step("STEP 4/6 Building global mutation matrix")
    marker_rows = []
    sample_lookup = high_quality_df.set_index("sample_id").to_dict(orient="index")
    selected_ids = high_quality_df["sample_id"].astype(str).drop_duplicates().tolist()

    with Progress(
        SpinnerColumn(),
        TextColumn("[bold yellow]{task.description}"),
        BarColumn(bar_width=40),
        TextColumn("[bold white]{task.completed}/{task.total}"),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
        console=console,
    ) as progress:
        task = progress.add_task("Searching mutations", total=len(selected_ids))
        for sid in selected_ids:
            meta_row = sample_lookup[sid]
            vcf_path = meta_row["vcf_path"]

            for m in marker_map.itertuples(index=False):
                hit_line = ""
                present = 0
                if m.mutation_aa3:
                    hit_line = grep_first_line(vcf_path, m.mutation_aa3)
                    present = int(bool(hit_line))

                marker_rows.append(
                    {
                        "sample_id": sid,
                        "sample": meta_row["sample"],
                        "state": meta_row["state"],
                        "year": meta_row["year"],
                        "age_group": meta_row["age_group"],
                        "drug": m.drug,
                        "gene": m.gene,
                        "mutation": m.mutation,
                        "mutation_aa3": m.mutation_aa3,
                        "mutation_present": present,
                        "hit_line": hit_line,
                    }
                )
            progress.advance(task)

    marker_long_df = pd.DataFrame(marker_rows).sort_values(["state", "sample_id", "drug", "gene", "mutation"])
    marker_long_tsv = Path(args.output_dir) / "step3_mutation_presence_long.tsv"
    marker_long_df.to_csv(marker_long_tsv, sep="\t", index=False)

    marker_matrix = marker_long_df.pivot_table(
        index=["sample_id", "sample", "state", "year", "age_group"],
        columns="mutation",
        values="mutation_present",
        aggfunc="first",
        fill_value=0,
    ).reset_index()
    marker_matrix.columns.name = None
    marker_matrix_tsv = Path(args.output_dir) / "step3_mutation_matrix_global_high_quality_samples.tsv"
    marker_matrix.to_csv(marker_matrix_tsv, sep="\t", index=False)

    mutation_prevalence = (
        marker_long_df.groupby(["drug", "gene", "mutation", "mutation_aa3"], as_index=False)
        .agg(
            numerator=("mutation_present", "sum"),
            denominator=("mutation_present", "size"),
        )
        .sort_values(["drug", "gene", "mutation"])
    )
    mutation_prevalence["prevalence"] = mutation_prevalence["numerator"] / mutation_prevalence["denominator"]
    mutation_prevalence_tsv = Path(args.output_dir) / "step3_mutation_prevalence_overall_high_quality_samples.tsv"
    mutation_prevalence.to_csv(mutation_prevalence_tsv, sep="\t", index=False)

    print_success("Global mutation matrix created")

    print_step("STEP 5/6 Building global haplotype matrix")
    hap_rows = []

    with Progress(
        SpinnerColumn(),
        TextColumn("[bold magenta]{task.description}"),
        BarColumn(bar_width=40),
        TextColumn("[bold white]{task.completed}/{task.total}"),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
        console=console,
    ) as progress:
        task = progress.add_task("Searching haplotypes", total=len(selected_ids))
        for sid in selected_ids:
            meta_row = sample_lookup[sid]
            vcf_path = meta_row["vcf_path"]

            for h in haplotypes.itertuples(index=False):
                component_statuses = []
                any_hit = ""
                present = 1

                for mut1, mut3 in zip(h.mutation_list, h.mutation_list_aa3):
                    hit_line = grep_first_line(vcf_path, mut3)
                    hit_flag = int(bool(hit_line))
                    component_statuses.append(f"{mut1}:{'present' if hit_flag == 1 else 'absent'}")
                    if hit_line and not any_hit:
                        any_hit = hit_line
                    if hit_flag == 0:
                        present = 0

                hap_rows.append(
                    {
                        "sample_id": sid,
                        "sample": meta_row["sample"],
                        "state": meta_row["state"],
                        "year": meta_row["year"],
                        "age_group": meta_row["age_group"],
                        "drug": h.drug,
                        "gene": h.gene,
                        "haplotype": h.haplotype,
                        "mutations": h.mutations,
                        "mutations_aa3": " + ".join(h.mutation_list_aa3),
                        "n_required_mutations": len(h.mutation_list),
                        "haplotype_present": present,
                        "component_statuses": ";".join(component_statuses),
                        "hit_line": any_hit if present == 1 else "",
                    }
                )
            progress.advance(task)

    hap_long_df = pd.DataFrame(hap_rows).sort_values(["state", "sample_id", "drug", "haplotype"])
    hap_long_tsv = Path(args.output_dir) / "step3_haplotype_presence_long.tsv"
    hap_long_df.to_csv(hap_long_tsv, sep="\t", index=False)

    hap_matrix = hap_long_df.pivot_table(
        index=["sample_id", "sample", "state", "year", "age_group"],
        columns="haplotype",
        values="haplotype_present",
        aggfunc="first",
        fill_value=0,
    ).reset_index()
    hap_matrix.columns.name = None
    hap_matrix_tsv = Path(args.output_dir) / "step3_haplotype_matrix_global_high_quality_samples.tsv"
    hap_matrix.to_csv(hap_matrix_tsv, sep="\t", index=False)

    hap_prevalence = (
        hap_long_df.groupby(["drug", "gene", "haplotype", "mutations", "mutations_aa3"], as_index=False)
        .agg(
            numerator=("haplotype_present", "sum"),
            denominator=("haplotype_present", "size"),
        )
        .sort_values(["drug", "haplotype"])
    )
    hap_prevalence["prevalence"] = hap_prevalence["numerator"] / hap_prevalence["denominator"]
    hap_prevalence_tsv = Path(args.output_dir) / "step3_haplotype_prevalence_overall_high_quality_samples.tsv"
    hap_prevalence.to_csv(hap_prevalence_tsv, sep="\t", index=False)

    print_success("Global haplotype matrix created")

    print_step("STEP 6/6 Writing manifest and QC outputs")
    qc_summary = pd.DataFrame(
        [
            {
                "n_total_samples_input": int(sample_df["sample_id"].nunique()),
                "n_samples_with_bam": int(sample_df[sample_df["bam_exists"]]["sample_id"].nunique()),
                "n_samples_with_vcf": int(sample_df[sample_df["vcf_exists"]]["sample_id"].nunique()),
                "n_high_quality_samples_all_fragments_callable": int(len(high_quality_samples)),
                "n_fragments_required": int(fragments["fragment_id"].nunique()),
                "n_mutations": int(marker_map["mutation"].nunique()),
                "n_haplotypes": int(haplotypes["haplotype"].nunique()),
            }
        ]
    )
    qc_summary_tsv = Path(args.output_dir) / "step3_qc_summary.tsv"
    qc_summary.to_csv(qc_summary_tsv, sep="\t", index=False)

    outputs = {
        "sample_manifest_tsv": str(sample_manifest_tsv),
        "fragment_coverage_by_sample_tsv": str(coverage_tsv),
        "fragment_coverage_summary_by_sample_tsv": str(coverage_summary_sample_tsv),
        "high_quality_samples_manifest_tsv": str(high_quality_manifest_tsv),
        "mutation_presence_long_tsv": str(marker_long_tsv),
        "mutation_matrix_global_tsv": str(marker_matrix_tsv),
        "mutation_prevalence_overall_tsv": str(mutation_prevalence_tsv),
        "haplotype_presence_long_tsv": str(hap_long_tsv),
        "haplotype_matrix_global_tsv": str(hap_matrix_tsv),
        "haplotype_prevalence_overall_tsv": str(hap_prevalence_tsv),
        "qc_summary_tsv": str(qc_summary_tsv),
    }

    manifest = {
        "base_dir": str(args.base_dir),
        "metadata_tsv": str(args.metadata_tsv),
        "fragment_manifest": str(args.fragment_manifest),
        "bed_file": str(args.bed_file),
        "haplotype_file": str(args.haplotype_file),
        "output_dir": str(args.output_dir),
        "min_reads": args.min_reads,
        "n_input_samples": int(sample_df["sample_id"].nunique()),
        "n_high_quality_samples": int(len(high_quality_samples)),
        "n_fragments": int(fragments["fragment_id"].nunique()),
        "n_mutations": int(marker_map["mutation"].nunique()),
        "n_haplotypes": int(haplotypes["haplotype"].nunique()),
        "outputs": outputs,
    }
    manifest_path = Path(args.output_dir) / "step3_pipeline_manifest.json"
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2)

    summary_table = Table(title="Step 3 summary", show_lines=False)
    summary_table.add_column("Metric", style="cyan")
    summary_table.add_column("Value", style="green")
    summary_table.add_row("Input samples", str(sample_df["sample_id"].nunique()))
    summary_table.add_row("High-quality samples", str(len(high_quality_samples)))
    summary_table.add_row("Fragments required", str(fragments["fragment_id"].nunique()))
    summary_table.add_row("Mutations", str(marker_map["mutation"].nunique()))
    summary_table.add_row("Haplotypes", str(haplotypes["haplotype"].nunique()))
    summary_table.add_row("Mutation matrix rows", str(len(marker_matrix)))
    summary_table.add_row("Haplotype matrix rows", str(len(hap_matrix)))
    console.print(summary_table)

    elapsed = time.time() - start_time
    console.print(
        Panel.fit(
            f"[bold green]STEP 3 PIPELINE COMPLETED SUCCESSFULLY[/bold green]\n"
            f"[bold white]Manifest:[/bold white] {manifest_path}\n"
            f"[bold white]Elapsed:[/bold white] {elapsed:.2f} seconds",
            border_style="bright_green",
        )
    )


if __name__ == "__main__":
    main()