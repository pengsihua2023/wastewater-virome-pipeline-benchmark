#!/usr/bin/env python3
"""
Simulate PacBio HiFi reads for a metagenomic community using PBSIM3.

Uses the same genome references and abundance distribution as the CAMISIM
Illumina/Nanopore simulations. Each genome is simulated independently with
depth proportional to its abundance, then all reads are merged.

Usage:
    python simulate_pacbio.py \
        --genome_locations config/genome_locations.tsv \
        --distribution config/distribution_0.txt \
        --output_dir out/wastewater_sim_pacbio \
        --total_bases 1000000000 \
        --seed 42
"""

import argparse
import csv
import gzip
import logging
import os
import subprocess
import sys
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
log = logging.getLogger(__name__)


def get_genome_size(fasta_path):
    """Calculate total bases in a FASTA file."""
    total = 0
    with open(fasta_path) as f:
        for line in f:
            if not line.startswith(">"):
                total += len(line.strip())
    return total


def main():
    parser = argparse.ArgumentParser(
        description="Simulate PacBio HiFi reads for a metagenomic community using PBSIM3"
    )
    parser.add_argument(
        "--genome_locations", required=True,
        help="TSV file: genome_ID<tab>path_to_fasta (no header)",
    )
    parser.add_argument(
        "--distribution", required=True,
        help="TSV file: genome_ID<tab>abundance (no header, normalized to 1.0)",
    )
    parser.add_argument(
        "--output_dir", default="out/wastewater_sim_pacbio",
        help="Output directory (default: out/wastewater_sim_pacbio)",
    )
    parser.add_argument(
        "--total_bases", type=float, default=1e9,
        help="Total bases to simulate in bp (default: 1e9 = 1 Gbp)",
    )
    parser.add_argument(
        "--model", default="",
        help="Path to PBSIM3 quality model file. If empty, auto-detect from pbsim3 installation.",
    )
    parser.add_argument(
        "--method", default="errhmm",
        help="PBSIM3 simulation method: errhmm or qshmm (default: errhmm)",
    )
    parser.add_argument(
        "--seed", type=int, default=42,
        help="Random seed (default: 42)",
    )
    parser.add_argument(
        "--threads", type=int, default=8,
        help="Number of parallel pbsim3 jobs (default: 8)",
    )
    parser.add_argument(
        "--mean_length", type=int, default=15000,
        help="Mean read length in bp (default: 15000 for HiFi)",
    )
    parser.add_argument(
        "--sd_length", type=int, default=2000,
        help="Standard deviation of read length (default: 2000)",
    )
    args = parser.parse_args()

    # Check pbsim3
    try:
        result = subprocess.run(["pbsim", "--version"], capture_output=True, text=True)
        log.info(f"pbsim3 found")
    except FileNotFoundError:
        log.error("pbsim3 not found. Install with: conda install -c bioconda pbsim3")
        sys.exit(1)

    # Find model file if not specified
    model_path = args.model
    if not model_path:
        # Try to find the model from pbsim3 installation
        pbsim_bin = subprocess.run(["which", "pbsim"], capture_output=True, text=True).stdout.strip()
        if pbsim_bin:
            # Models are typically in ../share/pbsim3/ relative to bin
            share_dir = Path(pbsim_bin).parent.parent / "share" / "pbsim3"
            if share_dir.exists():
                # Look for HiFi model (errhmm)
                if args.method == "errhmm":
                    candidates = sorted(share_dir.glob("ERRHMM-*.model"))
                    if not candidates:
                        candidates = sorted(share_dir.glob("*.model"))
                else:
                    candidates = sorted(share_dir.glob("QSHMM-*.model"))
                    if not candidates:
                        candidates = sorted(share_dir.glob("*.model"))

                if candidates:
                    model_path = str(candidates[0])
                    log.info(f"Auto-detected model: {model_path}")

        if not model_path:
            log.error(
                "Could not find PBSIM3 model file. Specify with --model. "
                "Models are usually in $(conda info --base)/share/pbsim3/"
            )
            sys.exit(1)

    # Read genome locations
    genomes = {}
    with open(args.genome_locations) as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if len(row) >= 2:
                genomes[row[0]] = row[1]

    # Read abundance distribution
    abundances = {}
    with open(args.distribution) as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if len(row) >= 2:
                abundances[row[0]] = float(row[1])

    # Normalize abundances
    total_ab = sum(abundances.values())
    abundances = {k: v / total_ab for k, v in abundances.items()}

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir = output_dir / "tmp_pbsim"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    log.info(f"Total bases to simulate: {args.total_bases:.0f} bp ({args.total_bases/1e9:.1f} Gbp)")
    log.info(f"Number of genomes: {len(genomes)}")
    log.info(f"Method: {args.method}, Model: {model_path}")
    log.info(f"Mean read length: {args.mean_length} bp")
    log.info("=" * 70)

    # Calculate depth for each genome based on abundance and genome size
    genome_tasks = []
    for genome_id in genomes:
        if genome_id not in abundances:
            log.warning(f"  {genome_id}: no abundance found, skipping")
            continue

        fasta_path = genomes[genome_id]
        if not os.path.isabs(fasta_path):
            fasta_path = os.path.join(os.getcwd(), fasta_path)

        if not os.path.exists(fasta_path):
            log.warning(f"  {genome_id}: FASTA not found at {fasta_path}, skipping")
            continue

        genome_size = get_genome_size(fasta_path)
        if genome_size == 0:
            log.warning(f"  {genome_id}: empty FASTA, skipping")
            continue

        # Bases allocated to this genome
        allocated_bases = args.total_bases * abundances[genome_id]
        # Depth = allocated_bases / genome_size
        depth = allocated_bases / genome_size

        if depth < 0.01:
            log.info(f"  {genome_id}: depth={depth:.4f}x (very low, skipping)")
            continue

        genome_tasks.append({
            "genome_id": genome_id,
            "fasta_path": fasta_path,
            "depth": depth,
            "abundance": abundances[genome_id],
            "genome_size": genome_size,
        })

    log.info(f"Genomes to simulate: {len(genome_tasks)}")
    log.info("=" * 70)

    # Run PBSIM3 for each genome
    fastq_files = []
    for i, task in enumerate(genome_tasks, 1):
        genome_id = task["genome_id"]
        prefix = tmp_dir / genome_id
        log.info(
            f"[{i}/{len(genome_tasks)}] {genome_id}: "
            f"abundance={task['abundance']:.4f}, "
            f"genome_size={task['genome_size']:,}bp, "
            f"depth={task['depth']:.2f}x"
        )

        # Skip if already simulated (resume support)
        existing = list(tmp_dir.glob(f"{genome_id}_*.fq.gz")) + list(tmp_dir.glob(f"{genome_id}_*.fastq.gz"))
        if existing:
            log.info(f"  Already simulated ({len(existing)} files), skipping")
            for fq in sorted(existing):
                fastq_files.append(fq)
            continue

        cmd = [
            "pbsim",
            "--strategy", "wgs",
            "--method", args.method,
            f"--{args.method}", model_path,
            "--depth", f"{task['depth']:.4f}",
            "--genome", task["fasta_path"],
            "--prefix", str(prefix),
            "--seed", str(args.seed + i),
            "--length-mean", str(args.mean_length),
            "--length-sd", str(args.sd_length),
        ]

        try:
            result = subprocess.run(
                cmd, capture_output=True, text=True, timeout=3600
            )
            if result.returncode != 0:
                log.warning(f"  PBSIM3 failed for {genome_id}: {result.stderr[:200]}")
                continue
        except subprocess.TimeoutExpired:
            log.warning(f"  PBSIM3 timed out for {genome_id}")
            continue

        # Collect output FASTQ files (pbsim3 outputs _0001.fq.gz, .fq, .fastq, or .fastq.gz)
        for pattern in [f"{genome_id}_*.fq.gz", f"{genome_id}_*.fastq.gz", f"{genome_id}_*.fq", f"{genome_id}_*.fastq"]:
            for fq in sorted(tmp_dir.glob(pattern)):
                fastq_files.append(fq)
                log.info(f"  Generated: {fq.name}")

    log.info("=" * 70)

    if not fastq_files:
        log.error("No FASTQ files generated!")
        sys.exit(1)

    # Merge all FASTQ files
    merged_fastq = output_dir / "pacbio_hifi_reads.fastq.gz"
    log.info(f"Merging {len(fastq_files)} FASTQ files into {merged_fastq}")

    total_reads = 0
    total_bases_out = 0
    with gzip.open(merged_fastq, "wt") as out:
        for fq in fastq_files:
            opener = gzip.open if str(fq).endswith(".gz") else open
            with opener(fq, "rt") as f:
                line_count = 0
                for line in f:
                    out.write(line)
                    line_count += 1
                    if line_count % 4 == 2:  # sequence line
                        total_reads += 1
                        total_bases_out += len(line.strip())

    log.info(f"Total reads: {total_reads:,}")
    log.info(f"Total bases: {total_bases_out:,} ({total_bases_out/1e9:.2f} Gbp)")
    log.info(f"Output: {merged_fastq}")

    # Write summary
    summary_file = output_dir / "simulation_summary.txt"
    with open(summary_file, "w") as f:
        f.write(f"PacBio HiFi Simulation Summary\n")
        f.write(f"{'='*50}\n")
        f.write(f"Total reads: {total_reads:,}\n")
        f.write(f"Total bases: {total_bases_out:,} ({total_bases_out/1e9:.2f} Gbp)\n")
        f.write(f"Mean read length setting: {args.mean_length} bp\n")
        f.write(f"Method: {args.method}\n")
        f.write(f"Model: {model_path}\n")
        f.write(f"Seed: {args.seed}\n")
        f.write(f"\nPer-genome details:\n")
        f.write(f"{'Genome ID':<15}{'Abundance':>12}{'Genome Size':>15}{'Depth':>10}\n")
        for task in genome_tasks:
            f.write(
                f"{task['genome_id']:<15}"
                f"{task['abundance']:>12.6f}"
                f"{task['genome_size']:>15,}"
                f"{task['depth']:>10.2f}x\n"
            )

    log.info(f"Summary: {summary_file}")
    log.info("Done!")


if __name__ == "__main__":
    main()
