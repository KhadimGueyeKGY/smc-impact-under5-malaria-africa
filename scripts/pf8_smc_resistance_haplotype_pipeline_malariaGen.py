#!/usr/bin/env python3

import argparse
import gzip
import json
import math
import re
import subprocess
import time
from collections import defaultdict
from pathlib import Path

import pandas as pd
from rich.console import Console
from rich.panel import Panel
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TimeElapsedColumn, TimeRemainingColumn
from rich.table import Table


DEFAULT_COUNTRIES = [
    "Burkina Faso",
    "Cameroon",
    "Chad",
    "Gambia",
    "Ghana",
    "Guinea",
    "Guinea-Bissau",
    "Mali",
    "Mauritania",
    "Niger",
    "Nigeria",
    "Senegal",
    "Togo",
    "Mozambique",
    "Uganda",
]

THREE_TO_ONE = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Cys": "C",
    "Gln": "Q",
    "Glu": "E",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Phe": "F",
    "Pro": "P",
    "Ser": "S",
    "Thr": "T",
    "Trp": "W",
    "Tyr": "Y",
    "Val": "V",
    "Ter": "*",
    "Stop": "*",
}

ONE_TO_THREE = {v: k for k, v in THREE_TO_ONE.items()}
ONE_TO_THREE["*"] = "Ter"

console = Console()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-dir", required=True)
    parser.add_argument("--malariagen-dir", required=True)
    parser.add_argument("--samples-metadata", required=True)
    parser.add_argument("--bed-file", required=True)
    parser.add_argument("--haplotype-file", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--countries", nargs="*", default=DEFAULT_COUNTRIES)
    parser.add_argument("--year-min", type=int, default=2000)
    parser.add_argument("--year-max", type=int, default=2100)
    parser.add_argument("--bcftools", default="bcftools")
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def ensure_dir(path):
    Path(path).mkdir(parents=True, exist_ok=True)


def file_exists_and_nonempty(path):
    p = Path(path)
    return p.exists() and p.stat().st_size > 0


def normalize_country_name(x):
    if pd.isna(x):
        return x
    x = str(x).strip()
    fixes = {
        "Burkina_Faso": "Burkina Faso",
        "Guinea_Bissau": "Guinea-Bissau",
    }
    return fixes.get(x, x)


def extract_gt(sample_blob):
    if sample_blob is None:
        return None
    sample_blob = str(sample_blob).strip()
    if sample_blob == "":
        return None
    return sample_blob.split(":", 1)[0]


def gt_has_mutation(gt):
    if gt is None:
        return None
    gt = str(gt).strip()
    if gt in {"./.", ".|.", ".", "././.", ".|.|."}:
        return None
    alleles = re.split(r"[\/|]", gt)
    alleles = [a.strip() for a in alleles if a.strip() != ""]
    if not alleles:
        return None
    if any(a == "." for a in alleles):
        return None
    if all(a == "0" for a in alleles):
        return 0
    return 1


def mutation_one_to_three(mut):
    m = re.fullmatch(r"([A-Z\*])(\d+)([A-Z\*])", str(mut).strip())
    if not m:
        return str(mut).strip()
    ref, pos, alt = m.groups()
    ref3 = ONE_TO_THREE.get(ref, ref)
    alt3 = ONE_TO_THREE.get(alt, alt)
    return f"{ref3}{pos}{alt3}"


def extract_ann_aa_change(info_field):
    if "ANN=" not in info_field:
        return set()
    ann_blob = info_field.split("ANN=", 1)[1].split(";", 1)[0]
    aa_changes = set()
    for ann in ann_blob.split(","):
        fields = ann.split("|")
        if len(fields) < 11:
            continue
        aa_change = fields[10].strip()
        if aa_change:
            aa_changes.add(aa_change)
    return aa_changes


def parse_p_change_to_three_letter(aa_change):
    aa_change = aa_change.strip()
    aa_change = aa_change.replace("p.", "")
    m = re.fullmatch(r"([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|=|Ter)", aa_change)
    if m:
        ref, pos, alt = m.groups()
        if alt == "=":
            return None
        return f"{ref}{pos}{alt}"
    return None


def sample_metadata_column_map(df):
    mapping = {}
    for c in df.columns:
        mapping[c.strip().lower()] = c
    return mapping


def print_header():
    console.print(
        Panel.fit(
            "[bold cyan]PF8 SMC Resistance Haplotype Pipeline[/bold cyan]\n"
            "[bold white]Targeted extraction, mutation prevalence, haplotype prevalence[/bold white]",
            border_style="bright_blue",
        )
    )


def print_step(title):
    console.print(Panel.fit(f"[bold white]{title}[/bold white]", border_style="green"))


def print_success(msg):
    console.print(f"[bold green]✔ {msg}[/bold green]")


def print_info(msg):
    console.print(f"[bold cyan]{msg}[/bold cyan]")


def print_command(cmd):
    console.print(f"[bold yellow]$ {' '.join(cmd)}[/bold yellow]")


def run(cmd):
    print_command(cmd)
    subprocess.run(cmd, check=True)


def load_samples_metadata(path, countries, year_min, year_max):
    samples = pd.read_csv(path, sep="\t", dtype=str)
    colmap = sample_metadata_column_map(samples)
    sample_col = colmap.get("sample")
    country_col = colmap.get("country")
    year_col = colmap.get("year")
    qc_col = colmap.get("qc pass")

    if sample_col is None or country_col is None or year_col is None:
        raise ValueError("Missing Sample, Country, or Year column in samples metadata")

    samples[country_col] = samples[country_col].apply(normalize_country_name)
    samples[year_col] = pd.to_numeric(samples[year_col], errors="coerce").astype("Int64")

    keep = samples[country_col].isin(countries) & samples[year_col].notna()
    keep &= samples[year_col].between(year_min, year_max)

    if qc_col is not None:
        samples[qc_col] = samples[qc_col].astype(str)
        keep &= samples[qc_col].isin(["True", "TRUE", "true"])

    out = samples.loc[keep].copy()
    out = out.rename(
        columns={
            sample_col: "sample",
            country_col: "country",
            year_col: "year",
        }
    )

    extra_cols = []
    for wanted in ["Study", "Admin level 1", "ENA", "Population", "Sample type", "Sample was in Pf7"]:
        if wanted in out.columns:
            extra_cols.append(wanted)

    out = out[["sample", "country", "year"] + extra_cols].drop_duplicates()
    out["year"] = out["year"].astype(int)
    return out


def load_bed(path):
    bed = pd.read_csv(
        path,
        sep="\t",
        header=None,
        comment="#",
        names=["chrom", "start", "end", "mutation", "gene", "drug"],
        dtype={"chrom": str, "start": int, "end": int, "mutation": str, "gene": str, "drug": str},
    )
    bed["mutation_1letter"] = bed["mutation"].astype(str).str.strip()
    bed["mutation_3letter"] = bed["mutation_1letter"].apply(mutation_one_to_three)
    bed["pos1"] = bed["start"].astype(int)
    return bed


def load_haplotypes(path):
    hap = pd.read_csv(path, sep="\t", dtype=str)
    hap.columns = [c.strip() for c in hap.columns]
    required = {"drug", "gene", "mutations", "haplotype"}
    if not required.issubset(set(hap.columns)):
        raise ValueError("Haplotype file missing required columns")
    hap["mutation_list_1letter"] = hap["mutations"].apply(
        lambda x: [m.strip() for m in re.split(r"\s*\+\s*", str(x).strip()) if m.strip()]
    )
    hap["mutation_list_3letter"] = hap["mutation_list_1letter"].apply(
        lambda muts: [mutation_one_to_three(m) for m in muts]
    )
    return hap


def write_samples_outputs(samples_df, output_dir):
    samples_tsv = Path(output_dir) / "pf8_smc_15countries_samples_2000_2024.tsv"
    samples_df.sort_values(["country", "year", "sample"]).to_csv(samples_tsv, sep="\t", index=False)

    sample_ids_txt = Path(output_dir) / "pf8_smc_15countries_sample_ids_2000_2024.txt"
    samples_df["sample"].drop_duplicates().to_csv(sample_ids_txt, sep="\t", index=False, header=False)

    counts_country_year = (
        samples_df.groupby(["country", "year"], as_index=False)
        .agg(n_samples=("sample", "nunique"))
        .sort_values(["country", "year"])
    )
    counts_country_year_tsv = Path(output_dir) / "pf8_smc_15countries_sample_counts_by_country_year.tsv"
    counts_country_year.to_csv(counts_country_year_tsv, sep="\t", index=False)

    counts_country = (
        samples_df.groupby(["country"], as_index=False)
        .agg(n_samples=("sample", "nunique"), year_min=("year", "min"), year_max=("year", "max"))
        .sort_values(["country"])
    )
    counts_country_tsv = Path(output_dir) / "pf8_smc_15countries_sample_counts_by_country.tsv"
    counts_country.to_csv(counts_country_tsv, sep="\t", index=False)

    return {
        "samples_tsv": str(samples_tsv),
        "sample_ids_txt": str(sample_ids_txt),
        "counts_country_year_tsv": str(counts_country_year_tsv),
        "counts_country_tsv": str(counts_country_tsv),
        "counts_country_year_df": counts_country_year,
        "counts_country_df": counts_country,
    }


def build_site_lists(bed_df, output_dir):
    regions_file = Path(output_dir) / "pf3d7_sp_aq_resistance_markers_regions_v1.tsv"
    bed_regions_file = Path(output_dir) / "pf3d7_sp_aq_resistance_markers_regions_v1.bed"
    per_chrom_json = Path(output_dir) / "pf3d7_sp_aq_resistance_markers_by_chromosome_v1.json"

    regions = bed_df[["chrom", "pos1"]].drop_duplicates().sort_values(["chrom", "pos1"]).copy()
    regions["region"] = regions["chrom"] + ":" + regions["pos1"].astype(str) + "-" + regions["pos1"].astype(str)
    regions[["region"]].to_csv(regions_file, sep="\t", index=False, header=False)

    bed_df[["chrom", "start", "end", "mutation", "gene", "drug"]].drop_duplicates().sort_values(
        ["chrom", "start", "mutation"]
    ).to_csv(bed_regions_file, sep="\t", index=False, header=False)

    chrom_map = defaultdict(list)
    for _, row in regions.iterrows():
        chrom_map[row["chrom"]].append(int(row["pos1"]))
    with open(per_chrom_json, "w") as f:
        json.dump({k: sorted(v) for k, v in sorted(chrom_map.items())}, f, indent=2)

    return {
        "regions_file": str(regions_file),
        "bed_regions_file": str(bed_regions_file),
        "per_chrom_json": str(per_chrom_json),
        "n_markers": int(bed_df[["chrom", "start"]].drop_duplicates().shape[0]),
        "n_chromosomes": int(bed_df["chrom"].nunique()),
    }


def extract_targeted_vcfs_per_chrom(args, bed_df, sample_ids_file, output_dir):
    targeted_dir = Path(output_dir) / "targeted_vcfs"
    ensure_dir(targeted_dir)

    metadata_rows = []
    chroms = sorted(bed_df["chrom"].unique())

    with Progress(
        SpinnerColumn(),
        TextColumn("[bold cyan]{task.description}"),
        BarColumn(bar_width=40),
        TextColumn("[bold white]{task.completed}/{task.total}"),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
        console=console,
    ) as progress:
        task = progress.add_task("Extracting targeted VCFs", total=len(chroms))

        for chrom in chroms:
            src_vcf = Path(args.malariagen_dir) / f"{chrom}.clean.snps.vcf.gz"
            if not src_vcf.exists():
                progress.console.print(f"[bold red]Missing source VCF:[/bold red] {src_vcf}")
                progress.advance(task)
                continue

            subset_vcf = targeted_dir / f"{chrom}.smc15.targeted.markers.vcf.gz"
            subset_index_tbi = targeted_dir / f"{chrom}.smc15.targeted.markers.vcf.gz.tbi"
            subset_index_csi = targeted_dir / f"{chrom}.smc15.targeted.markers.vcf.gz.csi"
            chrom_bed = targeted_dir / f"{chrom}.markers.bed"

            bed_df.loc[
                bed_df["chrom"] == chrom, ["chrom", "start", "end", "mutation", "gene", "drug"]
            ].drop_duplicates().sort_values(["start", "mutation"]).to_csv(
                chrom_bed, sep="\t", index=False, header=False
            )

            if args.force or not file_exists_and_nonempty(subset_vcf):
                cmd = [
                    args.bcftools,
                    "view",
                    "-R",
                    str(chrom_bed),
                    "-S",
                    str(sample_ids_file),
                    "-Oz",
                    "-o",
                    str(subset_vcf),
                    str(src_vcf),
                ]
                run(cmd)

            if args.force or (not subset_index_tbi.exists() and not subset_index_csi.exists()):
                run([args.bcftools, "index", "-f", str(subset_vcf)])

            metadata_rows.append(
                {
                    "chromosome": chrom,
                    "source_vcf": str(src_vcf),
                    "targeted_vcf": str(subset_vcf),
                    "markers_bed": str(chrom_bed),
                }
            )

            progress.advance(task)

    meta_df = pd.DataFrame(metadata_rows).sort_values(["chromosome"])
    meta_tsv = Path(output_dir) / "pf8_smc_15countries_targeted_vcf_manifest.tsv"
    meta_df.to_csv(meta_tsv, sep="\t", index=False)
    return meta_df, str(meta_tsv)


def parse_targeted_vcf(vcf_path, bed_df, samples_df):
    sample_info = samples_df[["sample", "country", "year"]].drop_duplicates().copy()
    sample_info = sample_info.set_index("sample")

    bed_lookup = {
        (str(r.chrom), int(r.pos1)): {
            "mutation_1letter": r.mutation_1letter,
            "mutation_3letter": r.mutation_3letter,
            "gene": r.gene,
            "drug": r.drug,
        }
        for r in bed_df.itertuples(index=False)
    }

    rows = []
    with gzip.open(vcf_path, "rt") as fh:
        sample_names = []
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                sample_names = header[9:]
                continue

            parts = line.rstrip("\n").split("\t")
            chrom = parts[0]
            pos = int(parts[1])
            info = parts[7]
            sample_values = parts[9:]

            key = (chrom, pos)
            if key not in bed_lookup:
                continue

            lookup = bed_lookup[key]
            ann_three_set = set()
            ann_one_set = set()

            for aa_change in extract_ann_aa_change(info):
                p3 = parse_p_change_to_three_letter(aa_change)
                if p3:
                    ann_three_set.add(p3)
                    m = re.fullmatch(r"([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|Ter)", p3)
                    if m:
                        ref3, aa_pos, alt3 = m.groups()
                        ann_one_set.add(f"{THREE_TO_ONE.get(ref3, ref3)}{aa_pos}{THREE_TO_ONE.get(alt3, alt3)}")

            for sample_name, sample_blob in zip(sample_names, sample_values):
                gt = extract_gt(sample_blob)
                mut_flag = gt_has_mutation(gt)
                if mut_flag is None:
                    continue
                if sample_name not in sample_info.index:
                    continue

                call_status = "mutant" if mut_flag == 1 else "wildtype"
                rows.append(
                    {
                        "sample": sample_name,
                        "country": sample_info.loc[sample_name, "country"],
                        "year": int(sample_info.loc[sample_name, "year"]),
                        "chrom": chrom,
                        "pos": pos,
                        "mutation_1letter": lookup["mutation_1letter"],
                        "mutation_3letter_expected": lookup["mutation_3letter"],
                        "mutation_3letter_ann": ";".join(sorted(ann_three_set)) if ann_three_set else "",
                        "mutation_1letter_ann": ";".join(sorted(ann_one_set)) if ann_one_set else "",
                        "gene": lookup["gene"],
                        "drug": lookup["drug"],
                        "gt": gt,
                        "mutation_present": int(mut_flag),
                        "call_status": call_status,
                    }
                )

    return pd.DataFrame(rows)


def write_long_calls(all_calls_df, output_dir):
    long_tsv = Path(output_dir) / "pf8_smc_15countries_targeted_marker_calls_long.tsv"
    all_calls_df.sort_values(
        ["country", "year", "sample", "chrom", "pos", "mutation_1letter"]
    ).to_csv(long_tsv, sep="\t", index=False)

    per_mut_sample = (
        all_calls_df.groupby(
            ["country", "year", "sample", "drug", "gene", "mutation_1letter", "mutation_3letter_expected"],
            as_index=False,
        )
        .agg(
            n_called_sites=("mutation_present", "size"),
            mutation_present=("mutation_present", "max"),
        )
        .sort_values(["country", "year", "sample", "drug", "gene", "mutation_1letter"])
    )
    per_mut_sample_tsv = Path(output_dir) / "pf8_smc_15countries_marker_presence_by_sample.tsv"
    per_mut_sample.to_csv(per_mut_sample_tsv, sep="\t", index=False)

    per_mut_country_year = (
        per_mut_sample.groupby(
            ["country", "year", "drug", "gene", "mutation_1letter", "mutation_3letter_expected"],
            as_index=False,
        )
        .agg(
            n_called_samples=("mutation_present", "size"),
            n_mutant_samples=("mutation_present", "sum"),
        )
        .sort_values(["country", "year", "drug", "gene", "mutation_1letter"])
    )
    per_mut_country_year["prevalence"] = (
        per_mut_country_year["n_mutant_samples"] / per_mut_country_year["n_called_samples"]
    )
    per_mut_country_year_tsv = Path(output_dir) / "pf8_smc_15countries_marker_prevalence_by_country_year.tsv"
    per_mut_country_year.to_csv(per_mut_country_year_tsv, sep="\t", index=False)

    return {
        "long_tsv": str(long_tsv),
        "per_mut_sample_tsv": str(per_mut_sample_tsv),
        "per_mut_country_year_tsv": str(per_mut_country_year_tsv),
        "per_mut_sample_df": per_mut_sample,
        "per_mut_country_year_df": per_mut_country_year,
    }


def compute_haplotype_prevalence(per_mut_sample_df, hap_df):
    records = []
    per_mut_sample_df = per_mut_sample_df.copy()
    mut_map = defaultdict(dict)

    for row in per_mut_sample_df.itertuples(index=False):
        key = (row.country, int(row.year), row.sample)
        mut_map[key][str(row.mutation_1letter)] = int(row.mutation_present)

    sample_meta = per_mut_sample_df[["country", "year", "sample"]].drop_duplicates()

    grouped_country_year = list(sample_meta.groupby(["country", "year"]))
    with Progress(
        SpinnerColumn(),
        TextColumn("[bold magenta]{task.description}"),
        BarColumn(bar_width=40),
        TextColumn("[bold white]{task.completed}/{task.total}"),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
        console=console,
    ) as progress:
        task = progress.add_task("Computing haplotype prevalence", total=len(grouped_country_year))
        for (country, year), grp in grouped_country_year:
            samples = grp["sample"].tolist()
            for hr in hap_df.itertuples(index=False):
                muts = [str(m).strip() for m in hr.mutation_list_1letter]
                n_called_samples = 0
                n_haplotype_positive = 0
                n_missing_any_required = 0

                for sample in samples:
                    state_map = mut_map.get((country, int(year), sample), {})
                    values = [state_map.get(m, None) for m in muts]
                    if any(v is None for v in values):
                        n_missing_any_required += 1
                        continue
                    n_called_samples += 1
                    if all(v == 1 for v in values):
                        n_haplotype_positive += 1

                prevalence = (
                    n_haplotype_positive / n_called_samples
                    if n_called_samples > 0
                    else math.nan
                )

                records.append(
                    {
                        "country": country,
                        "year": int(year),
                        "drug": hr.drug,
                        "gene": hr.gene,
                        "haplotype": hr.haplotype,
                        "mutations": hr.mutations,
                        "n_required_mutations": len(muts),
                        "n_total_samples": len(samples),
                        "n_called_samples": n_called_samples,
                        "n_missing_any_required": n_missing_any_required,
                        "n_haplotype_positive": n_haplotype_positive,
                        "prevalence": prevalence,
                    }
                )
            progress.advance(task)

    return pd.DataFrame(records).sort_values(["country", "year", "drug", "haplotype"])


def write_haplotype_outputs(haplotype_df, output_dir):
    hap_tsv = Path(output_dir) / "pf8_smc_15countries_haplotype_prevalence_by_country_year.tsv"
    haplotype_df.to_csv(hap_tsv, sep="\t", index=False)

    plot_tsv = Path(output_dir) / "pf8_smc_15countries_haplotype_prevalence_for_plot.tsv"
    haplotype_df.to_csv(plot_tsv, sep="\t", index=False)

    summary = (
        haplotype_df.groupby(["drug", "haplotype"], as_index=False)
        .agg(
            countries_with_data=("country", "nunique"),
            years_with_data=("year", "nunique"),
            mean_prevalence=("prevalence", "mean"),
        )
        .sort_values(["drug", "haplotype"])
    )
    summary_tsv = Path(output_dir) / "pf8_smc_15countries_haplotype_prevalence_summary.tsv"
    summary.to_csv(summary_tsv, sep="\t", index=False)

    return {
        "hap_tsv": str(hap_tsv),
        "plot_tsv": str(plot_tsv),
        "summary_tsv": str(summary_tsv),
        "summary_df": summary,
    }


def write_mutation_plot_outputs(per_mut_country_year_df, output_dir):
    mutation_plot_tsv = Path(output_dir) / "pf8_smc_15countries_mutation_prevalence_for_plot.tsv"
    per_mut_country_year_df.to_csv(mutation_plot_tsv, sep="\t", index=False)

    wide = per_mut_country_year_df.pivot_table(
        index=["country", "year"],
        columns="mutation_1letter",
        values="prevalence",
        aggfunc="first",
    ).reset_index()
    wide.columns.name = None
    wide_tsv = Path(output_dir) / "pf8_smc_15countries_mutation_prevalence_wide.tsv"
    wide.to_csv(wide_tsv, sep="\t", index=False)

    return {
        "mutation_plot_tsv": str(mutation_plot_tsv),
        "wide_tsv": str(wide_tsv),
    }


def write_manifest(output_dir, args, outputs):
    manifest = {
        "project_dir": args.project_dir,
        "malariagen_dir": args.malariagen_dir,
        "samples_metadata": args.samples_metadata,
        "bed_file": args.bed_file,
        "haplotype_file": args.haplotype_file,
        "output_dir": args.output_dir,
        "countries": args.countries,
        "year_min": args.year_min,
        "year_max": args.year_max,
        "outputs": outputs,
    }
    manifest_path = Path(output_dir) / "pf8_smc_15countries_pipeline_manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)
    return str(manifest_path)


def show_samples_summary(samples_df, counts_country_df):
    table = Table(title="Samples retained", show_lines=False)
    table.add_column("Country", style="cyan")
    table.add_column("Samples", justify="right", style="green")
    table.add_column("Min year", justify="right", style="white")
    table.add_column("Max year", justify="right", style="white")

    for row in counts_country_df.itertuples(index=False):
        table.add_row(str(row.country), str(row.n_samples), str(row.year_min), str(row.year_max))

    console.print(table)
    console.print(f"[bold white]Total retained samples:[/bold white] [bold green]{samples_df['sample'].nunique()}[/bold green]")


def show_run_summary(site_outputs, targeted_manifest_df, all_calls_df, hap_summary_df):
    table = Table(title="Pipeline summary", show_lines=False)
    table.add_column("Metric", style="cyan")
    table.add_column("Value", style="green")

    table.add_row("Markers", str(site_outputs["n_markers"]))
    table.add_row("Chromosomes", str(site_outputs["n_chromosomes"]))
    table.add_row("Targeted VCF files", str(len(targeted_manifest_df)))
    table.add_row("Marker calls", str(len(all_calls_df)))
    table.add_row("Haplotype summaries", str(len(hap_summary_df)))

    console.print(table)


def main():
    args = parse_args()
    start_time = time.time()

    print_header()
    ensure_dir(args.output_dir)

    print_step("STEP 1/6 Filtering samples")
    samples_df = load_samples_metadata(
        args.samples_metadata,
        [normalize_country_name(c) for c in args.countries],
        args.year_min,
        args.year_max,
    )
    sample_outputs = write_samples_outputs(samples_df, args.output_dir)
    show_samples_summary(samples_df, sample_outputs["counts_country_df"])
    print_success("Sample filtering completed")

    print_step("STEP 2/6 Loading marker and haplotype panels")
    bed_df = load_bed(args.bed_file)
    hap_df = load_haplotypes(args.haplotype_file)
    site_outputs = build_site_lists(bed_df, args.output_dir)
    print_info(f"Markers loaded: {site_outputs['n_markers']}")
    print_info(f"Chromosomes loaded: {site_outputs['n_chromosomes']}")
    print_success("Panel loading completed")

    print_step("STEP 3/6 Extracting targeted VCFs")
    targeted_manifest_df, targeted_manifest_tsv = extract_targeted_vcfs_per_chrom(
        args,
        bed_df,
        sample_outputs["sample_ids_txt"],
        args.output_dir,
    )
    print_success("Targeted VCF extraction completed")

    print_step("STEP 4/6 Parsing targeted VCFs")
    call_tables = []
    with Progress(
        SpinnerColumn(),
        TextColumn("[bold yellow]{task.description}"),
        BarColumn(bar_width=40),
        TextColumn("[bold white]{task.completed}/{task.total}"),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
        console=console,
    ) as progress:
        task = progress.add_task("Parsing chromosomes", total=len(targeted_manifest_df))
        for row in targeted_manifest_df.itertuples(index=False):
            call_df = parse_targeted_vcf(row.targeted_vcf, bed_df, samples_df)
            if not call_df.empty:
                call_tables.append(call_df)
            progress.advance(task)

    if call_tables:
        all_calls_df = pd.concat(call_tables, ignore_index=True)
    else:
        all_calls_df = pd.DataFrame(
            columns=[
                "sample",
                "country",
                "year",
                "chrom",
                "pos",
                "mutation_1letter",
                "mutation_3letter_expected",
                "mutation_3letter_ann",
                "mutation_1letter_ann",
                "gene",
                "drug",
                "gt",
                "mutation_present",
                "call_status",
            ]
        )
    print_success("Targeted VCF parsing completed")

    print_step("STEP 5/6 Writing mutation prevalence outputs")
    call_outputs = write_long_calls(all_calls_df, args.output_dir)
    mut_plot_outputs = write_mutation_plot_outputs(
        call_outputs["per_mut_country_year_df"],
        args.output_dir,
    )
    print_success("Mutation prevalence outputs completed")

    print_step("STEP 6/6 Computing haplotype prevalence")
    haplotype_df = compute_haplotype_prevalence(call_outputs["per_mut_sample_df"], hap_df)
    hap_outputs = write_haplotype_outputs(haplotype_df, args.output_dir)
    print_success("Haplotype prevalence outputs completed")

    outputs = {}
    outputs.update({k: v for k, v in sample_outputs.items() if not k.endswith("_df")})
    outputs.update(site_outputs)
    outputs["targeted_manifest_tsv"] = targeted_manifest_tsv
    outputs.update({k: v for k, v in call_outputs.items() if not k.endswith("_df")})
    outputs.update({k: v for k, v in hap_outputs.items() if not k.endswith("_df")})
    outputs.update(mut_plot_outputs)

    manifest_path = write_manifest(args.output_dir, args, outputs)

    show_run_summary(site_outputs, targeted_manifest_df, all_calls_df, hap_outputs["summary_df"])

    elapsed = time.time() - start_time
    console.print(
        Panel.fit(
            f"[bold green]PIPELINE COMPLETED SUCCESSFULLY[/bold green]\n"
            f"[bold white]Manifest:[/bold white] {manifest_path}\n"
            f"[bold white]Elapsed:[/bold white] {elapsed:.2f} seconds",
            border_style="bright_green",
        )
    )


if __name__ == "__main__":
    main()