#!/bin/bash
#SBATCH --job-name=Dinucleotide
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --time=2:00:00
#SBATCH --mem=64G
#SBATCH --output=%x_%j.out
#SBATCH --mail-type=END,FAIL


set -euo pipefail


# ============================================================
# MODULES
# ============================================================

module purge
module load Python/3.14.2-GCCcore-15.2.0
module load BEDTools/2.31.1-GCC-13.3.0


# ============================================================
# WORK DIRECTORY
# ============================================================

WORKDIR="/scratch/sj31280/06_projects/NHF1_UV_Target_Damage-seq_workflow2_sh/Final/03.Multi_Analysis/02.Dinu_seq"
GROUP_DIR="/scratch/sj31280/06_projects/NHF1_UV_Target_Damage-seq_workflow2_sh/Final/02.Cluster/02.Target_4_group"
FASTA="/work/wllab/shan/05.scripts/Reference/hg19_male/refined.hg19.male.fa"

mkdir -p "$WORKDIR"
cd "$WORKDIR"


BED_OUT="${WORKDIR}/01.Cluster_bed"
FASTA_OUT="${WORKDIR}/02.Cluster_fasta"
DINU_OUT="${WORKDIR}/03.Dinucleotide_results"

mkdir -p "$BED_OUT"
mkdir -p "$FASTA_OUT"
mkdir -p "$DINU_OUT"


declare -A GROUP_FILES

GROUP_FILES[M_H_C_H]="${GROUP_DIR}/Mutation_High_CPD_High.bed"
GROUP_FILES[M_H_C_L]="${GROUP_DIR}/Mutation_High_CPD_Low.bed"
GROUP_FILES[M_L_C_H]="${GROUP_DIR}/Mutation_Low_CPD_High.bed"
GROUP_FILES[M_L_C_L]="${GROUP_DIR}/Mutation_Low_CPD_Low.bed"

GROUP_ORDER=(
    M_H_C_H
    M_H_C_L
    M_L_C_H
    M_L_C_L
)


if [ ! -f "$FASTA" ]; then
    echo "ERROR: reference FASTA not found:"
    echo "$FASTA"
    exit 1
fi

for GROUP in "${GROUP_ORDER[@]}"; do

    if [ ! -f "${GROUP_FILES[$GROUP]}" ]; then
        echo "ERROR: group BED file not found:"
        echo "${GROUP_FILES[$GROUP]}"
        exit 1
    fi

done



for GROUP in "${GROUP_ORDER[@]}"; do

    INPUT_BED="${GROUP_FILES[$GROUP]}"
    OUTPUT_BED="${BED_OUT}/${GROUP}.bed"


    awk 'BEGIN{OFS="\t"} {
        name=$4 "::" $1 ":" $2 "-" $3
        print $1,$2,$3,name
    }' "$INPUT_BED" \
    > "$OUTPUT_BED"

    echo "Created:"
    echo "$OUTPUT_BED"

done



for GROUP in "${GROUP_ORDER[@]}"; do

    INPUT_BED="${BED_OUT}/${GROUP}.bed"
    OUTPUT_FASTA="${FASTA_OUT}/${GROUP}.fasta"


    bedtools getfasta \
        -fi "$FASTA" \
        -bed "$INPUT_BED" \
        -name \
        > "$OUTPUT_FASTA"

    echo "Created:"
    echo "$OUTPUT_FASTA"

done


python "${WORKDIR}/Each_dinucleotide_percentage.py"



FINAL_RESULT="${DINU_OUT}/Dinucleotide_percentage_4groups.tsv"

if [ -f "$FINAL_RESULT" ]; then


    cat "$FINAL_RESULT"

    awk '
    NR > 1 {
        g1 += $2
        g2 += $3
        g3 += $4
        g4 += $5
    }
    END {
        printf "M_H_C_H\t%.6f\n", g1
        printf "M_H_C_L\t%.6f\n", g2
        printf "M_L_C_H\t%.6f\n", g3
        printf "M_L_C_L\t%.6f\n", g4
    }
    ' "$FINAL_RESULT"

else

    echo "ERROR: final dinucleotide result was not generated:"
    echo "$FINAL_RESULT"
    exit 1

fi
