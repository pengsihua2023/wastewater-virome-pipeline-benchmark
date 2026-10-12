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
import shutil
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


def _bam_to_fastq(bam, fq_out, args):
    """Convert a subread BAM to FASTQ using ccs (HiFi) or samtools (CLR)."""
    if fq_out.exists():
        log.info(f"  Already converted: {fq_out.name}")
        return [fq_out]

    use_ccs = (args.method == "qshmm" and args.pass_num > 1)
    if use_ccs:
        min_passes = min(3, args.pass_num)
        cmd = [
            "ccs", str(bam), str(fq_out),
            "--min-passes", str(min_passes),
            "--min-rq", "0.99",
            "--min-length", str(args.min_len),
            "--num-threads", str(args.threads),
        ]
        label = "ccs"
    else:
        cmd = ["samtools", "fastq", "-0", str(fq_out), str(bam)]
        label = "samtools fastq"

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        log.warning(f"  {label} failed for {bam.name}: {result.stderr[:300]}")
        return []

    if fq_out.exists() and fq_out.stat().st_size > 0:
        log.info(f"  {label} → {fq_out.name}")
        return [fq_out]
    else:
        log.warning(f"  {label} produced empty output for {bam.name}")
        return []


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
        "--method", default="qshmm",
        help="PBSIM3 simulation method: errhmm (CLR) or qshmm (HiFi) (default: qshmm)",
    )
    parser.add_argument(
        "--pass_num", type=int, default=10,
        help="Number of CCS passes for HiFi simulation (default: 10, only used with qshmm)",
    )
    parser.add_argument(
        "--min_depth", type=float, default=1.0,
        help="Minimum sequencing depth per genome (default: 1.0x).",
    )
    parser.add_argument(
        "--min_reads", type=int, default=5,
        help="Minimum number of reads per genome (default: 5). "
             "For genomes smaller than mean_length, this overrides min_depth "
             "to ensure every genome produces at least this many reads.",
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
    parser.add_argument(
        "--min_len", type=int, default=1000,
        help="Minimum read length to keep after merging (default: 1000 bp)",
    )
    parser.add_argument(
        "--keep_tmp", action="store_true", default=False,
        help="Keep tmp_pbsim/ directory after merging (default: delete it to save disk space)",
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
    log.info(f"Pass num: {args.pass_num}, Mean read length: {args.mean_length} bp")
    log.info(f"Min depth: {args.min_depth}x, Min reads per genome: {args.min_reads}")
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

        # Ensure minimum depth, and also ensure minimum reads regardless of genome size
        min_depth_for_reads = (args.min_reads * args.mean_length) / genome_size
        effective_min = max(args.min_depth, min_depth_for_reads)
        if depth < effective_min:
            log.info(
                f"  {genome_id}: depth={depth:.4f}x → raised to {effective_min:.4f}x "
                f"(min_depth={args.min_depth}x, min_reads={args.min_reads} × "
                f"{args.mean_length}bp / {genome_size}bp = {min_depth_for_reads:.4f}x)"
            )
            depth = effective_min

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

        # Skip if already simulated (resume support).
        existing_fq = (list(tmp_dir.glob(f"{genome_id}_*.fq.gz")) +
                       list(tmp_dir.glob(f"{genome_id}_*.fastq.gz")))
        existing_bam = list(tmp_dir.glob(f"{genome_id}_*.bam"))
        if existing_fq:
            log.info(f"  Already converted ({len(existing_fq)} FASTQ files), skipping")
            for fq in sorted(existing_fq):
                fastq_files.append(fq)
            continue
        if existing_bam:
            log.info(f"  BAM exists ({len(existing_bam)} files), converting to FASTQ...")
            for bam in sorted(existing_bam):
                fq_out = bam.with_suffix(".fastq.gz")
                fqs = _bam_to_fastq(bam, fq_out, args)
                fastq_files.extend(fqs)
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
            "--length-min", str(args.min_len),
        ]
        # qshmm + pass-num enables CCS consensus → true HiFi quality (Q20+)
        if args.method == "qshmm" and args.pass_num > 1:
            cmd += ["--pass-num", str(args.pass_num)]

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

        # PBSIM3 with --pass-num outputs subread BAM; convert via ccs (HiFi) or samtools (CLR).
        bam_files = sorted(tmp_dir.glob(f"{genome_id}_*.bam"))
        if bam_files:
            for bam in bam_files:
                fq_out = bam.with_suffix(".fastq.gz")
                fqs = _bam_to_fastq(bam, fq_out, args)
                fastq_files.extend(fqs)
        else:
            for pattern in [f"{genome_id}_*.fq.gz", f"{genome_id}_*.fastq.gz",
                            f"{genome_id}_*.fq",    f"{genome_id}_*.fastq"]:
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
    skipped_reads = 0
    log.info(f"Minimum read length filter: {args.min_len} bp")
    with gzip.open(merged_fastq, "wt") as out:
        for fq in fastq_files:
            opener = gzip.open if str(fq).endswith(".gz") else open
            with opener(fq, "rt") as f:
                while True:
                    header = f.readline()
                    if not header:
                        break
                    seq = f.readline()
                    plus = f.readline()
                    qual = f.readline()
                    read_len = len(seq.strip())
                    if read_len < args.min_len:
                        skipped_reads += 1
                        continue
                    out.write(header)
                    out.write(seq)
                    out.write(plus)
                    out.write(qual)
                    total_reads += 1
                    total_bases_out += read_len

    log.info(f"Total reads: {total_reads:,}")
    log.info(f"Skipped reads (< {args.min_len} bp): {skipped_reads:,}")
    log.info(f"Total bases: {total_bases_out:,} ({total_bases_out/1e9:.2f} Gbp)")
    log.info(f"Output: {merged_fastq}")

    # Write summary
    summary_file = output_dir / "simulation_summary.txt"
    with open(summary_file, "w") as f:
        f.write(f"PacBio HiFi Simulation Summary\n")
        f.write(f"{'='*50}\n")
        f.write(f"Total reads: {total_reads:,}\n")
        f.write(f"Skipped reads (< {args.min_len} bp): {skipped_reads:,}\n")
        f.write(f"Total bases: {total_bases_out:,} ({total_bases_out/1e9:.2f} Gbp)\n")
        f.write(f"Mean read length setting: {args.mean_length} bp\n")
        f.write(f"Min read length filter: {args.min_len} bp\n")
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

    # Clean up tmp directory unless --keep_tmp is set
    if args.keep_tmp:
        log.info(f"tmp_pbsim kept at: {tmp_dir}  (--keep_tmp)")
    else:
        log.info(f"Removing tmp_pbsim directory ({tmp_dir}) ...")
        shutil.rmtree(tmp_dir)
        log.info("tmp_pbsim removed.")

    log.info("Done!")


if __name__ == "__main__":
    main()
