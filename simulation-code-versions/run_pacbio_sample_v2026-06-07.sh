#!/bin/bash
#SBATCH --job-name=pacbio_10e9
#SBATCH --partition=bahl_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=72:00:00
#SBATCH --array=0-2
#SBATCH --output=pacbio_10e9_%A_%a.out
#SBATCH --error=pacbio_10e9_%A_%a.err

set -euo pipefail

source /apps/eb/Miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate pacbio_sim
export PATH=~/.conda/envs/samtools_env/bin:$PATH
PROJECT=/scratch/$USER/Meta-genome-data-analysis/CAMISIM_project

SCRIPT=$PROJECT/simulate_pacbio.py
GENOME_LOC=$PROJECT/camisim_input/genome_locations.tsv
MODEL=/home/$USER/.conda/envs/pacbio_sim/data/QSHMM-RSII.model

i=$SLURM_ARRAY_TASK_ID
OUTDIR=$PROJECT/pacbio_hifi_output/sample${i}
mkdir -p "$OUTDIR"

echo "samtools: $(which samtools) -> $(samtools --version | head -1)"
echo "model: $MODEL"

python "$SCRIPT" \
 --genome_locations "$GENOME_LOC" \
 --distribution     "$PROJECT/camisim_input/distribution_${i}.txt" \
 --output_dir       "$OUTDIR" \
 --total_bases 10e9 \
 --method      qshmm \
 --model       "$MODEL" \
 --pass_num    10 \
 --min_depth   0 \
 --min_reads   1 \
 --mean_length 15000 \
 --sd_length   2000 \
 --min_len     1000 \
 --seed        $((42 + i)) \
 --threads     ${SLURM_CPUS_PER_TASK}
