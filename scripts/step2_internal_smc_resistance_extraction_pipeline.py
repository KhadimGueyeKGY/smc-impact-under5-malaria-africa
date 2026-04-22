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
            "[bold cyan]STEP 2 - INTERNAL SMC SAMPLE EXTRACTION PIPELINE[/bold cyan]\n"
            "[bold white]Coverage, markers, haplotypes, prevalence, QC, plot-ready tables[/bold white]",
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


def mut_3to1(m):
    m = str(m).strip().replace("p.", "")
    g = re.match(r"^([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|Ter)$", m)
    if not g:
        return None
    return f"{AA3_TO_AA1[g.group(1)]}{g.group(2)}{AA3_TO_AA1[g.group(3)]}"


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
            fragment_start = None
            fragment_end = None
        else:
            best = g.sort_values(["start", "end"]).iloc[0]
            fragment_id = best["fragment_id"]
            fragment_start = int(best["start"])
            fragment_end = int(best["end"])

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
                "fragment_start": fragment_start,
                "fragment_end": fragment_end,
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


def grep_first_any_line(vcf_path, tokens):
    if vcf_path is None or not Path(vcf_path).exists():
        return ""
    tokens = [str(t).strip() for t in tokens if str(t).strip()]
    if not tokens:
        return ""
    pattern = "|".join(re.escape(t) for t in tokens)
    q_vcf = shlex.quote(str(vcf_path))
    q_pat = shlex.quote(pattern)
    cmd = f"zcat {q_vcf} | grep -E {q_pat} | head -1"
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

    print_step("STEP 1/8 Loading metadata and references")
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

    print_step("STEP 2/8 Resolving sample records and file paths")
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
    sample_manifest_tsv = Path(args.output_dir) / "step2_sample_manifest.tsv"
    sample_manifest.to_csv(sample_manifest_tsv, sep="\t", index=False)
    print_info(f"Samples with BAM: {int(sample_df['bam_exists'].sum())}")
    print_info(f"Samples with VCF: {int(sample_df['vcf_exists'].sum())}")
    print_success("Sample paths resolved")

    print_step("STEP 3/8 Computing fragment coverage for all samples")
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
    coverage_long_tsv = Path(args.output_dir) / "step2_fragment_coverage_by_sample.tsv"
    coverage_df.to_csv(coverage_long_tsv, sep="\t", index=False)

    coverage_summary_state = (
        coverage_df.groupby(["state", "reference_id"], as_index=False)
        .agg(
            n_samples=("sample_id", "nunique"),
            n_callable_samples=("fragment_pass_100reads", "sum"),
            mean_reads=("n_reads", "mean"),
            median_reads=("n_reads", "median"),
            min_reads=("n_reads", "min"),
            max_reads=("n_reads", "max"),
        )
        .sort_values(["state", "reference_id"])
    )
    coverage_summary_state["callable_fraction"] = coverage_summary_state["n_callable_samples"] / coverage_summary_state["n_samples"]
    coverage_summary_state_tsv = Path(args.output_dir) / "step2_fragment_coverage_summary_by_state.tsv"
    coverage_summary_state.to_csv(coverage_summary_state_tsv, sep="\t", index=False)

    coverage_summary_sample = (
        coverage_df.groupby(["sample_id", "state"], as_index=False)
        .agg(
            n_fragments=("fragment_id", "nunique"),
            n_callable_fragments=("fragment_pass_100reads", "sum"),
            mean_reads=("n_reads", "mean"),
            median_reads=("n_reads", "median"),
        )
        .sort_values(["state", "sample_id"])
    )
    coverage_summary_sample["callable_fraction"] = coverage_summary_sample["n_callable_fragments"] / coverage_summary_sample["n_fragments"]
    coverage_summary_sample_tsv = Path(args.output_dir) / "step2_fragment_coverage_summary_by_sample.tsv"
    coverage_summary_sample.to_csv(coverage_summary_sample_tsv, sep="\t", index=False)
    print_success("Fragment coverage computed")

    print_step("STEP 4/8 Searching target mutations in VCFs")
    coverage_lookup = coverage_df[["sample_id", "fragment_id", "reference_id", "n_reads", "fragment_pass_100reads"]].drop_duplicates()
    sample_meta_df = sample_df.drop_duplicates(subset=["sample_id"], keep="first").copy()
    sample_meta_lookup = sample_meta_df.set_index("sample_id").to_dict(orient="index")
    sample_ids = sample_df["sample_id"].astype(str).drop_duplicates().tolist()

    observed_rows = []
    marker_rows = []

    with Progress(
        SpinnerColumn(),
        TextColumn("[bold yellow]{task.description}"),
        BarColumn(bar_width=40),
        TextColumn("[bold white]{task.completed}/{task.total}"),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
        console=console,
    ) as progress:
        task = progress.add_task("Searching markers", total=len(sample_ids))
        for sid in sample_ids:
            meta_row = sample_meta_lookup[sid]
            vcf_path = meta_row["vcf_path"]
            vcf_exists = int(bool(meta_row["vcf_exists"]))

            for m in marker_map.itertuples(index=False):
                cov = coverage_lookup[
                    (coverage_lookup["sample_id"] == sid)
                    & (coverage_lookup["fragment_id"] == m.fragment_id)
                ]
                if cov.empty:
                    n_reads = None
                    pass100 = 0
                else:
                    n_reads = cov["n_reads"].iloc[0]
                    pass100 = int(cov["fragment_pass_100reads"].iloc[0])

                callable_flag = int(pass100 == 1 and vcf_exists == 1 and pd.notna(m.fragment_id))
                hit_line = ""
                present = None
                if callable_flag == 1 and m.mutation_aa3:
                    hit_line = grep_first_line(vcf_path, m.mutation_aa3)
                    present = int(bool(hit_line))

                if callable_flag == 0:
                    call_status = "not_callable"
                else:
                    call_status = "present" if present == 1 else "absent"

                marker_rows.append(
                    {
                        "sample": meta_row["sample"],
                        "sample_id": sid,
                        "folder": meta_row["folder"],
                        "barcode": meta_row["barcode"],
                        "state": meta_row["state"],
                        "year": meta_row["year"],
                        "age_group": meta_row["age_group"],
                        "drug": m.drug,
                        "gene": m.gene,
                        "mutation": m.mutation,
                        "mutation_aa3": m.mutation_aa3,
                        "chromosome": m.chromosome,
                        "marker_start": m.marker_start,
                        "marker_end": m.marker_end,
                        "fragment_id": m.fragment_id,
                        "fragment_start": m.fragment_start,
                        "fragment_end": m.fragment_end,
                        "fragment_reads": n_reads,
                        "fragment_pass_100reads": pass100,
                        "vcf_exists": vcf_exists,
                        "callable": callable_flag,
                        "mutation_present": present,
                        "call_status": call_status,
                        "hit_line": hit_line,
                    }
                )

                if callable_flag == 1 and present == 1:
                    observed_rows.append(
                        {
                            "sample": meta_row["sample"],
                            "sample_id": sid,
                            "folder": meta_row["folder"],
                            "barcode": meta_row["barcode"],
                            "state": meta_row["state"],
                            "year": meta_row["year"],
                            "age_group": meta_row["age_group"],
                            "vcf_path": str(vcf_path),
                            "chromosome": m.chromosome,
                            "pos": m.marker_start,
                            "ref": "",
                            "alt": "",
                            "gene": m.gene,
                            "mutation": m.mutation,
                            "aa_change": m.mutation_aa3,
                            "impact": "",
                            "effect": "",
                            "qual": "",
                            "filter": "",
                            "dp": None,
                            "ad": "",
                            "af": None,
                            "hit_line": hit_line,
                        }
                    )
            progress.advance(task)

    observed_df = pd.DataFrame(observed_rows)
    observed_long_tsv = Path(args.output_dir) / "step2_observed_target_gene_mutations_by_sample_long.tsv"
    observed_df.to_csv(observed_long_tsv, sep="\t", index=False)

    marker_df = pd.DataFrame(marker_rows).sort_values(["state", "sample_id", "drug", "gene", "mutation"])
    marker_long_tsv = Path(args.output_dir) / "step2_marker_calls_by_sample_long.tsv"
    marker_df.to_csv(marker_long_tsv, sep="\t", index=False)
    print_success("Target mutation search completed")

    print_step("STEP 5/8 Writing marker prevalence outputs")
    marker_prevalence_state = (
        marker_df[marker_df["callable"] == 1]
        .groupby(["state", "drug", "gene", "mutation", "mutation_aa3"], as_index=False)
        .agg(
            numerator=("mutation_present", "sum"),
            denominator=("mutation_present", "size"),
        )
        .sort_values(["state", "drug", "gene", "mutation"])
    )
    marker_prevalence_state["prevalence"] = marker_prevalence_state["numerator"] / marker_prevalence_state["denominator"]
    marker_prevalence_state_tsv = Path(args.output_dir) / "step2_marker_prevalence_by_state.tsv"
    marker_prevalence_state.to_csv(marker_prevalence_state_tsv, sep="\t", index=False)

    marker_prevalence_overall = (
        marker_df[marker_df["callable"] == 1]
        .groupby(["drug", "gene", "mutation", "mutation_aa3"], as_index=False)
        .agg(
            numerator=("mutation_present", "sum"),
            denominator=("mutation_present", "size"),
        )
        .sort_values(["drug", "gene", "mutation"])
    )
    marker_prevalence_overall["prevalence"] = marker_prevalence_overall["numerator"] / marker_prevalence_overall["denominator"]
    marker_prevalence_overall_tsv = Path(args.output_dir) / "step2_marker_prevalence_overall.tsv"
    marker_prevalence_overall.to_csv(marker_prevalence_overall_tsv, sep="\t", index=False)

    marker_wide = marker_df.pivot_table(
        index=["sample_id", "sample", "state", "year", "age_group"],
        columns="mutation",
        values="mutation_present",
        aggfunc="first",
    ).reset_index()
    marker_wide.columns.name = None
    marker_wide_tsv = Path(args.output_dir) / "step2_marker_matrix_by_sample.tsv"
    marker_wide.to_csv(marker_wide_tsv, sep="\t", index=False)

    marker_heatmap_state = marker_prevalence_state.pivot_table(
        index="state",
        columns="mutation",
        values="prevalence",
        aggfunc="first",
    ).reset_index()
    marker_heatmap_state.columns.name = None
    marker_heatmap_state_tsv = Path(args.output_dir) / "step2_marker_heatmap_by_state.tsv"
    marker_heatmap_state.to_csv(marker_heatmap_state_tsv, sep="\t", index=False)

    marker_callability_state = (
        marker_df.groupby(["state", "mutation"], as_index=False)
        .agg(
            n_samples=("sample_id", "nunique"),
            n_callable=("callable", "sum"),
            n_present=("mutation_present", lambda x: int(pd.Series(x).fillna(0).sum())),
        )
        .sort_values(["state", "mutation"])
    )
    marker_callability_state["callable_fraction"] = marker_callability_state["n_callable"] / marker_callability_state["n_samples"]
    marker_callability_state_tsv = Path(args.output_dir) / "step2_marker_callability_by_state.tsv"
    marker_callability_state.to_csv(marker_callability_state_tsv, sep="\t", index=False)
    print_success("Marker tables created")

    print_step("STEP 6/8 Searching haplotypes in VCFs")
    fragment_gene_lookup = fragments.groupby("reference_id", as_index=False).first()[["reference_id", "chromosome", "start", "end", "fragment_id"]]
    gene_fragment_map = {
        row.reference_id: {
            "chromosome": row.chromosome,
            "start": int(row.start),
            "end": int(row.end),
            "fragment_id": row.fragment_id,
        }
        for row in fragment_gene_lookup.itertuples(index=False)
    }

    hap_rows = []
    with Progress(
        SpinnerColumn(),
        TextColumn("[bold cyan]{task.description}"),
        BarColumn(bar_width=40),
        TextColumn("[bold white]{task.completed}/{task.total}"),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
        console=console,
    ) as progress:
        task = progress.add_task("Searching haplotypes", total=len(sample_ids))
        for sid in sample_ids:
            meta_row = sample_meta_lookup[sid]
            vcf_path = meta_row["vcf_path"]
            vcf_exists = int(bool(meta_row["vcf_exists"]))

            for h in haplotypes.itertuples(index=False):
                component_mutations = [m.strip() for m in h.mutation_list if m.strip()]
                component_mutations_aa3 = [m for m in h.mutation_list_aa3 if m]

                gene_pass_rows = []
                all_genes_callable = True
                for gene_name in h.gene_list:
                    frag_meta = gene_fragment_map.get(gene_name)
                    if frag_meta is None:
                        all_genes_callable = False
                        gene_pass_rows.append((gene_name, None, 0))
                        continue
                    cov = coverage_lookup[
                        (coverage_lookup["sample_id"] == sid)
                        & (coverage_lookup["fragment_id"] == frag_meta["fragment_id"])
                    ]
                    if cov.empty:
                        gene_reads = None
                        gene_pass = 0
                    else:
                        gene_reads = cov["n_reads"].iloc[0]
                        gene_pass = int(cov["fragment_pass_100reads"].iloc[0])
                    gene_pass_rows.append((gene_name, gene_reads, gene_pass))
                    if gene_pass != 1:
                        all_genes_callable = False

                callable_flag = int(vcf_exists == 1 and all_genes_callable and len(component_mutations_aa3) > 0)
                component_hits = []
                component_statuses = []
                present = None
                any_hit_line = ""

                if callable_flag == 1:
                    for mut1, mut3 in zip(component_mutations, component_mutations_aa3):
                        hit_line = grep_first_line(vcf_path, mut3)
                        hit_flag = int(bool(hit_line))
                        component_hits.append((mut1, mut3, hit_flag, hit_line))
                        component_statuses.append(f"{mut1}:{'present' if hit_flag == 1 else 'absent'}")
                        if hit_line and not any_hit_line:
                            any_hit_line = hit_line
                    present = int(all(x[2] == 1 for x in component_hits))
                if callable_flag == 0:
                    call_status = "not_callable"
                else:
                    call_status = "present" if present == 1 else "absent"

                hap_rows.append(
                    {
                        "sample": meta_row["sample"],
                        "sample_id": sid,
                        "folder": meta_row["folder"],
                        "barcode": meta_row["barcode"],
                        "state": meta_row["state"],
                        "year": meta_row["year"],
                        "age_group": meta_row["age_group"],
                        "drug": h.drug,
                        "gene": h.gene,
                        "haplotype": h.haplotype,
                        "mutations": h.mutations,
                        "mutations_aa3": " + ".join(component_mutations_aa3),
                        "n_required_mutations": len(component_mutations),
                        "vcf_exists": vcf_exists,
                        "callable": callable_flag,
                        "haplotype_present": present,
                        "call_status": call_status,
                        "component_statuses": ";".join(component_statuses),
                        "hit_line": any_hit_line if present == 1 else "",
                        "gene_callable_status": ";".join([f"{g}:{r}:{p}" for g, r, p in gene_pass_rows]),
                    }
                )
            progress.advance(task)

    hap_df = pd.DataFrame(hap_rows).sort_values(["state", "sample_id", "drug", "haplotype"])
    hap_long_tsv = Path(args.output_dir) / "step2_haplotype_calls_by_sample_long.tsv"
    hap_df.to_csv(hap_long_tsv, sep="\t", index=False)

    hap_prevalence_state = (
        hap_df[hap_df["callable"] == 1]
        .groupby(["state", "drug", "gene", "haplotype", "mutations", "mutations_aa3"], as_index=False)
        .agg(
            numerator=("haplotype_present", "sum"),
            denominator=("haplotype_present", "size"),
        )
        .sort_values(["state", "drug", "haplotype"])
    )
    hap_prevalence_state["prevalence"] = hap_prevalence_state["numerator"] / hap_prevalence_state["denominator"]
    hap_prevalence_state_tsv = Path(args.output_dir) / "step2_haplotype_prevalence_by_state.tsv"
    hap_prevalence_state.to_csv(hap_prevalence_state_tsv, sep="\t", index=False)

    hap_prevalence_overall = (
        hap_df[hap_df["callable"] == 1]
        .groupby(["drug", "gene", "haplotype", "mutations", "mutations_aa3"], as_index=False)
        .agg(
            numerator=("haplotype_present", "sum"),
            denominator=("haplotype_present", "size"),
        )
        .sort_values(["drug", "haplotype"])
    )
    hap_prevalence_overall["prevalence"] = hap_prevalence_overall["numerator"] / hap_prevalence_overall["denominator"]
    hap_prevalence_overall_tsv = Path(args.output_dir) / "step2_haplotype_prevalence_overall.tsv"
    hap_prevalence_overall.to_csv(hap_prevalence_overall_tsv, sep="\t", index=False)

    hap_wide = hap_df.pivot_table(
        index=["sample_id", "sample", "state", "year", "age_group"],
        columns="haplotype",
        values="haplotype_present",
        aggfunc="first",
    ).reset_index()
    hap_wide.columns.name = None
    hap_wide_tsv = Path(args.output_dir) / "step2_haplotype_matrix_by_sample.tsv"
    hap_wide.to_csv(hap_wide_tsv, sep="\t", index=False)

    hap_heatmap_state = hap_prevalence_state.pivot_table(
        index="state",
        columns="haplotype",
        values="prevalence",
        aggfunc="first",
    ).reset_index()
    hap_heatmap_state.columns.name = None
    hap_heatmap_state_tsv = Path(args.output_dir) / "step2_haplotype_heatmap_by_state.tsv"
    hap_heatmap_state.to_csv(hap_heatmap_state_tsv, sep="\t", index=False)
    print_success("Haplotype tables created")

    print_step("STEP 7/8 Building target mutation summary tables and workbooks")
    if observed_df.empty:
        exploratory_sample_level = pd.DataFrame(
            columns=["state", "gene", "mutation", "n_positive_samples", "n_samples", "prevalence_sample_level"]
        )
        exploratory_overall = pd.DataFrame(
            columns=["gene", "mutation", "n_positive_samples", "n_samples", "prevalence_sample_level"]
        )
    else:
        exploratory_sample_level = (
            observed_df.groupby(["state", "gene", "mutation"], as_index=False)
            .agg(n_positive_samples=("sample_id", "nunique"))
            .sort_values(["state", "gene", "mutation"])
        )
        state_sample_denoms = sample_df.groupby("state", as_index=False).agg(n_samples=("sample_id", "nunique"))
        exploratory_sample_level = exploratory_sample_level.merge(state_sample_denoms, on="state", how="left")
        exploratory_sample_level["prevalence_sample_level"] = exploratory_sample_level["n_positive_samples"] / exploratory_sample_level["n_samples"]

        exploratory_overall = (
            observed_df.groupby(["gene", "mutation"], as_index=False)
            .agg(n_positive_samples=("sample_id", "nunique"))
            .sort_values(["gene", "mutation"])
        )
        exploratory_overall["n_samples"] = sample_df["sample_id"].nunique()
        exploratory_overall["prevalence_sample_level"] = exploratory_overall["n_positive_samples"] / exploratory_overall["n_samples"]

    exploratory_sample_tsv = Path(args.output_dir) / "step2_exploratory_observed_mutation_prevalence_by_state.tsv"
    exploratory_sample_level.to_csv(exploratory_sample_tsv, sep="\t", index=False)

    exploratory_overall_tsv = Path(args.output_dir) / "step2_exploratory_observed_mutation_prevalence_overall.tsv"
    exploratory_overall.to_csv(exploratory_overall_tsv, sep="\t", index=False)

    excel_path = Path(args.output_dir) / "step2_state_resistance_summary.xlsx"
    with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
        for state in sorted(sample_df["state"].dropna().unique()):
            state_marker = marker_prevalence_state[marker_prevalence_state["state"] == state].copy()
            state_hap = hap_prevalence_state[hap_prevalence_state["state"] == state].copy()
            state_cov = coverage_summary_state[coverage_summary_state["state"] == state].copy()
            state_exp = exploratory_sample_level[exploratory_sample_level["state"] == state].copy()
            if not state_marker.empty:
                state_marker.to_excel(writer, sheet_name=f"{state[:20]}_marker", index=False)
            if not state_hap.empty:
                state_hap.to_excel(writer, sheet_name=f"{state[:20]}_hap", index=False)
            if not state_cov.empty:
                state_cov.to_excel(writer, sheet_name=f"{state[:20]}_cov", index=False)
            if not state_exp.empty:
                state_exp.to_excel(writer, sheet_name=f"{state[:20]}_obs", index=False)
    print_success("Summary tables created")

    print_step("STEP 8/8 Writing QC summaries, plot manifests, and final outputs")
    qc_summary = (
        sample_df.groupby("state", as_index=False)
        .agg(
            n_selected_samples=("sample_id", "nunique"),
            n_bam_found=("bam_exists", "sum"),
            n_vcf_found=("vcf_exists", "sum"),
        )
        .sort_values("state")
    )
    qc_summary_tsv = Path(args.output_dir) / "step2_qc_summary_by_state.tsv"
    qc_summary.to_csv(qc_summary_tsv, sep="\t", index=False)

    inclusion_rows = []
    for state in sorted(sample_df["state"].dropna().unique()):
        st_samples = sample_df[sample_df["state"] == state]["sample_id"].nunique()
        st_bam = sample_df[(sample_df["state"] == state) & (sample_df["bam_exists"])]["sample_id"].nunique()
        st_vcf = sample_df[(sample_df["state"] == state) & (sample_df["vcf_exists"])]["sample_id"].nunique()
        st_cov = coverage_df[(coverage_df["state"] == state) & (coverage_df["fragment_pass_100reads"] == 1)]["sample_id"].nunique()
        st_marker = marker_df[(marker_df["state"] == state) & (marker_df["callable"] == 1)]["sample_id"].nunique()
        st_hap = hap_df[(hap_df["state"] == state) & (hap_df["callable"] == 1)]["sample_id"].nunique()
        inclusion_rows.append(
            {
                "state": state,
                "n_selected_samples": st_samples,
                "n_with_bam": st_bam,
                "n_with_vcf": st_vcf,
                "n_with_any_callable_fragment": st_cov,
                "n_with_any_callable_marker": st_marker,
                "n_with_any_callable_haplotype": st_hap,
            }
        )
    inclusion_df = pd.DataFrame(inclusion_rows).sort_values("state")
    inclusion_tsv = Path(args.output_dir) / "step2_inclusion_flow_by_state.tsv"
    inclusion_df.to_csv(inclusion_tsv, sep="\t", index=False)

    outputs = {
        "sample_manifest_tsv": str(sample_manifest_tsv),
        "fragment_coverage_by_sample_tsv": str(coverage_long_tsv),
        "fragment_coverage_summary_by_state_tsv": str(coverage_summary_state_tsv),
        "fragment_coverage_summary_by_sample_tsv": str(coverage_summary_sample_tsv),
        "observed_target_gene_mutations_by_sample_long_tsv": str(observed_long_tsv),
        "marker_calls_by_sample_long_tsv": str(marker_long_tsv),
        "marker_prevalence_by_state_tsv": str(marker_prevalence_state_tsv),
        "marker_prevalence_overall_tsv": str(marker_prevalence_overall_tsv),
        "marker_matrix_by_sample_tsv": str(marker_wide_tsv),
        "marker_heatmap_by_state_tsv": str(marker_heatmap_state_tsv),
        "marker_callability_by_state_tsv": str(marker_callability_state_tsv),
        "haplotype_calls_by_sample_long_tsv": str(hap_long_tsv),
        "haplotype_prevalence_by_state_tsv": str(hap_prevalence_state_tsv),
        "haplotype_prevalence_overall_tsv": str(hap_prevalence_overall_tsv),
        "haplotype_matrix_by_sample_tsv": str(hap_wide_tsv),
        "haplotype_heatmap_by_state_tsv": str(hap_heatmap_state_tsv),
        "exploratory_observed_mutation_prevalence_by_state_tsv": str(exploratory_sample_tsv),
        "exploratory_observed_mutation_prevalence_overall_tsv": str(exploratory_overall_tsv),
        "state_resistance_summary_xlsx": str(excel_path),
        "qc_summary_by_state_tsv": str(qc_summary_tsv),
        "inclusion_flow_by_state_tsv": str(inclusion_tsv),
    }

    manifest = {
        "base_dir": str(args.base_dir),
        "metadata_tsv": str(args.metadata_tsv),
        "fragment_manifest": str(args.fragment_manifest),
        "bed_file": str(args.bed_file),
        "haplotype_file": str(args.haplotype_file),
        "output_dir": str(args.output_dir),
        "min_reads": args.min_reads,
        "n_samples_selected": int(sample_df["sample_id"].nunique()),
        "n_states": int(sample_df["state"].nunique()),
        "n_fragments": int(fragments["fragment_id"].nunique()),
        "n_markers": int(marker_map["mutation"].nunique()),
        "n_haplotypes": int(haplotypes["haplotype"].nunique()),
        "outputs": outputs,
    }
    manifest_path = Path(args.output_dir) / "step2_pipeline_manifest.json"
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2)

    summary_table = Table(title="Pipeline summary", show_lines=False)
    summary_table.add_column("Metric", style="cyan")
    summary_table.add_column("Value", style="green")
    summary_table.add_row("Selected samples", str(sample_df["sample_id"].nunique()))
    summary_table.add_row("States", str(sample_df["state"].nunique()))
    summary_table.add_row("Fragments", str(fragments["fragment_id"].nunique()))
    summary_table.add_row("Markers", str(marker_map["mutation"].nunique()))
    summary_table.add_row("Haplotypes", str(haplotypes["haplotype"].nunique()))
    summary_table.add_row("Observed mutation rows", str(len(observed_df)))
    summary_table.add_row("Marker call rows", str(len(marker_df)))
    summary_table.add_row("Haplotype call rows", str(len(hap_df)))
    console.print(summary_table)

    elapsed = time.time() - start_time
    console.print(
        Panel.fit(
            f"[bold green]STEP 2 PIPELINE COMPLETED SUCCESSFULLY[/bold green]\n"
            f"[bold white]Manifest:[/bold white] {manifest_path}\n"
            f"[bold white]Elapsed:[/bold white] {elapsed:.2f} seconds",
            border_style="bright_green",
        )
    )


if __name__ == "__main__":
    main()