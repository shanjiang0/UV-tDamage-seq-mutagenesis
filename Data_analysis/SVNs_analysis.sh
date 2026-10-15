#!/bin/bash
set -euo pipefail

WORKDIR="/scratch/sj31280/06_projects/NHF1_UV_Target_Damage-seq_workflow2_sh/Final/03.Multi_Analysis/01.SVNs"

mkdir -p "$WORKDIR"
cd "$WORKDIR"

GROUP_DIR="/scratch/sj31280/06_projects/NHF1_UV_Target_Damage-seq_workflow2_sh/Final/02.Cluster/02.Target_4_group"

MUTATION="/scratch/sj31280/01_database/MELA-AU/Processed/ICGC_snvs_unique_by_coordinate_sample.bed"

CPD="/scratch/sj31280/06_projects/NHF1_UV_Target_Damage-seq_workflow2_sh/Final/02.Cluster/01.Input_data/CPD_vs_CTR_Fold.in_regions.bedgraph"


OUTDIR="${WORKDIR}/Results"
mkdir -p "$OUTDIR"


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

if [ ! -f "$MUTATION" ]; then
    echo "ERROR: Mutation file not found:"
    echo "$MUTATION"
    exit 1
fi

if [ ! -f "$CPD" ]; then
    echo "ERROR: CPD file not found:"
    echo "$CPD"
    exit 1
fi

for GROUP in "${GROUP_ORDER[@]}"; do

    if [ ! -f "${GROUP_FILES[$GROUP]}" ]; then
        echo "ERROR: Group file not found:"
        echo "${GROUP_FILES[$GROUP]}"
        exit 1
    fi

done


for GROUP in "${GROUP_ORDER[@]}"; do

    CLUSTER="${GROUP_FILES[$GROUP]}"


    # --------------------------------------------------------
    # Count number of clusters
    # --------------------------------------------------------

    N_CLUSTER=$(awk 'NF >= 3 {n++} END {print n+0}' "$CLUSTER")

    echo "Number of clusters: $N_CLUSTER"


    # ========================================================
    # PART 1: MAP MUTATIONS TO CLUSTERS
    # ========================================================

    bedtools intersect \
        -a "$CLUSTER" \
        -b "$MUTATION" \
        -wa -wb \
        > "${OUTDIR}/${GROUP}.cluster_mutation_overlap.tsv"


    # ========================================================
    # PART 2: CALCULATE SNV TYPE FREQUENCY
    # ========================================================


    awk -v N="$N_CLUSTER" '
    BEGIN {

        OFS="\t"

        type[1]="A>C"
        type[2]="A>G"
        type[3]="A>T"
        type[4]="C>A"
        type[5]="C>G"
        type[6]="C>T"
        type[7]="G>A"
        type[8]="G>C"
        type[9]="G>T"
        type[10]="T>A"
        type[11]="T>C"
        type[12]="T>G"
    }

    FNR==NR {
        cluster=$1 ":" $2 "-" $3
        cluster_len[cluster]=$3-$2
        next
    }

    {
        cluster=$1 ":" $2 "-" $3
        ref=$10
        alt=$11
        snv=ref ">" alt

        if (ref ~ /^[ACGT]$/ && alt ~ /^[ACGT]$/ && ref != alt) {
            count[cluster,snv]++
        }
    }

    END {

        for (i=1; i<=12; i++) {
            t=type[i]
            sum_norm=0
            for (cluster in cluster_len) {
                sum_norm += count[cluster,t] / cluster_len[cluster]
            }

            if (N > 0) {
                frequency=sum_norm/N
            } else {
                frequency=0
            }
            printf "%s\t%.6f\n", t, frequency
        }
    }
    ' \
    "$CLUSTER" \
    "${OUTDIR}/${GROUP}.cluster_mutation_overlap.tsv" \
    > "${OUTDIR}/${GROUP}.SNV_frequency.tsv"


    # ========================================================
    # PART 3: UNIQUE MUTATION SITES WITHIN EACH GROUP
    # ========================================================

    bedtools intersect \
        -a "$MUTATION" \
        -b "$CLUSTER" \
        -wa -u \
    | cut -f1-3 \
    | sort -k1,1V -k2,2n -k3,3n -u \
    > "${OUTDIR}/${GROUP}.unique_mutation_sites.bed"


    TOTAL_UNIQUE=$(wc -l < "${OUTDIR}/${GROUP}.unique_mutation_sites.bed")


    # ========================================================
    # PART 4: EXACT MUTATION-CPD OVERLAP
    # ========================================================

    bedtools intersect \
        -a "${OUTDIR}/${GROUP}.unique_mutation_sites.bed" \
        -b "$CPD" \
        -u \
        > "${OUTDIR}/${GROUP}.unique_mutation_CPD_overlap.bed"


    OVERLAP_UNIQUE=$(wc -l < "${OUTDIR}/${GROUP}.unique_mutation_CPD_overlap.bed")


    if [ "$TOTAL_UNIQUE" -gt 0 ]; then

        OVERLAP_PERCENT=$(awk \
            -v x="$OVERLAP_UNIQUE" \
            -v n="$TOTAL_UNIQUE" \
            'BEGIN {printf "%.6f", (x/n)*100}')

    else

        OVERLAP_PERCENT="0.000000"

    fi


    # ========================================================
    # PART 5: OUTPUT EXACT OVERLAP WITH CPD SIGNAL
    # ========================================================


    bedtools intersect \
        -a "${OUTDIR}/${GROUP}.unique_mutation_sites.bed" \
        -b "$CPD" \
        -wa -wb \
        > "${OUTDIR}/${GROUP}.mutation_CPD_exact_overlap_with_signal.tsv"


# ============================================================
# COMBINE 4 GROUP SNV FREQUENCIES
# ============================================================

paste \
    "${OUTDIR}/M_H_C_H.SNV_frequency.tsv" \
    "${OUTDIR}/M_H_C_L.SNV_frequency.tsv" \
    "${OUTDIR}/M_L_C_H.SNV_frequency.tsv" \
    "${OUTDIR}/M_L_C_L.SNV_frequency.tsv" \
| awk '
BEGIN {
    OFS="\t"
}
{
    print $1,$2,$4,$6,$8
}
' \
>> "${OUTDIR}/SNV_frequency_4groups.tsv"


