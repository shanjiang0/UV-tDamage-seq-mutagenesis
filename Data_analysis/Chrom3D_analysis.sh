#!/bin/bash
#SBATCH --job-name=Chrom3D
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --mem=32G
#SBATCH --output=%x_%j.out
#SBATCH --mail-type=END,FAIL

set -euo pipefail

module purge
module load BEDTools/2.31.1-GCC-13.3.0


INPUT_DIR="/scratch/sj31280/06_projects/NHF1_UV_Target_Damage-seq_workflow2_sh/Final/03.Multi_Analysis/06.Chrom3D/01.Input_data"
OUTDIR="/scratch/sj31280/06_projects/NHF1_UV_Target_Damage-seq_workflow2_sh/Final/03.Multi_Analysis/06.Chrom3D/02.Result"
GROUP_DIR="/scratch/sj31280/06_projects/NHF1_UV_Target_Damage-seq_workflow2_sh/Final/02.Cluster/02.Target_4_group"

mkdir -p "$OUTDIR"

A_INFO="${INPUT_DIR}/A_info.bed"
B_INFO="${INPUT_DIR}/B_info.bed"


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


if [ ! -f "$A_INFO" ]; then
    echo "ERROR: A_info.bed not found:"
    echo "$A_INFO"
    exit 1
fi

if [ ! -f "$B_INFO" ]; then
    echo "ERROR: B_info.bed not found:"
    echo "$B_INFO"
    exit 1
fi

for GROUP in "${GROUP_ORDER[@]}"; do

    if [ ! -f "${GROUP_FILES[$GROUP]}" ]; then
        echo "ERROR: group file not found:"
        echo "${GROUP_FILES[$GROUP]}"
        exit 1
    fi

done


calculate_distance() {

    CLUSTER_FILE=$1
    DISTANCE_FILE=$2
    OUTPUT_FILE=$3
    LABEL=$4


    INTERSECT_TMP="${OUTPUT_FILE}.intersect.tmp"


    bedtools intersect \
        -a "$CLUSTER_FILE" \
        -b "$DISTANCE_FILE" \
        -wao \
        > "$INTERSECT_TMP"


    awk -v label="$LABEL" '
    BEGIN {
        OFS="\t"

        print "Chr","Start","End","Gene",
              "ClusterLength",
              "Covered_bp",
              "WeightedSum",
              label
    }

    {
        chr=$1
        start=$2
        end=$3
        gene=$4

        key=chr OFS start OFS end OFS gene

        cluster_len=end-start

        distance=$10
        overlap=$11

        clen[key]=cluster_len
        seen[key]=1

        if (overlap > 0 && distance != ".") {

            weighted[key] += overlap * distance
            covered[key] += overlap

        }
    }

    END {

        for (key in seen) {

            if (clen[key] > 0) {
                result=weighted[key]/clen[key]
            } else {
                result=0
            }

            printf "%s\t%d\t%d\t%.10f\t%.10f\n",
                   key,
                   clen[key],
                   covered[key],
                   weighted[key],
                   result
        }
    }
    ' "$INTERSECT_TMP" \
    | {
        read header
        echo "$header"
        sort -k1,1V -k2,2n -k3,3n
    } \
    > "$OUTPUT_FILE"


    rm -f "$INTERSECT_TMP"
}


for GROUP in "${GROUP_ORDER[@]}"; do

    CLUSTER="${GROUP_FILES[$GROUP]}"


    calculate_distance \
        "$CLUSTER" \
        "$A_INFO" \
        "${OUTDIR}/${GROUP}.A_distance.tsv" \
        "A_Distance"


    calculate_distance \
        "$CLUSTER" \
        "$B_INFO" \
        "${OUTDIR}/${GROUP}.B_distance.tsv" \
        "B_Distance"


    paste \
        "${OUTDIR}/${GROUP}.A_distance.tsv" \
        "${OUTDIR}/${GROUP}.B_distance.tsv" \
    | awk '
    BEGIN {
        OFS="\t"
    }

    NR==1 {
        print "Chr","Start","End","Gene",
              "ClusterLength",
              "A_Covered_bp",
              "A_Distance",
              "B_Covered_bp",
              "B_Distance"
        next
    }

    {
        print $1,$2,$3,$4,$5,$6,$8,$14,$16
    }
    ' \
    > "${OUTDIR}/${GROUP}.Chrom3D.tsv"


    echo "Created:"
    echo "${OUTDIR}/${GROUP}.Chrom3D.tsv"

done


SUMMARY="${OUTDIR}/Chrom3D_4groups_summary.tsv"

echo -e "Group\tCluster_Number\tMean_A_Distance\tMean_B_Distance" \
> "$SUMMARY"


for GROUP in "${GROUP_ORDER[@]}"; do

    FILE="${OUTDIR}/${GROUP}.Chrom3D.tsv"

    awk -v group="$GROUP" '
    BEGIN {
        OFS="\t"
    }

    NR > 1 {
        n++
        sumA += $7
        sumB += $9
    }

    END {

        if (n > 0) {
            printf "%s\t%d\t%.10f\t%.10f\n",
                   group,
                   n,
                   sumA/n,
                   sumB/n
        } else {
            print group,0,0,0
        }
    }
    ' "$FILE" >> "$SUMMARY"

done


WIDE="${OUTDIR}/Chrom3D_mean_distance_4groups.tsv"

echo -e "Feature\tM_H_C_H\tM_H_C_L\tM_L_C_H\tM_L_C_L" \
> "$WIDE"


# A distance
AHH=$(awk '$1=="M_H_C_H" {print $3}' "$SUMMARY")
AHL=$(awk '$1=="M_H_C_L" {print $3}' "$SUMMARY")
ALH=$(awk '$1=="M_L_C_H" {print $3}' "$SUMMARY")
ALL=$(awk '$1=="M_L_C_L" {print $3}' "$SUMMARY")

echo -e "A_Distance\t${AHH}\t${AHL}\t${ALH}\t${ALL}" \
>> "$WIDE"


# B distance
BHH=$(awk '$1=="M_H_C_H" {print $4}' "$SUMMARY")
BHL=$(awk '$1=="M_H_C_L" {print $4}' "$SUMMARY")
BLH=$(awk '$1=="M_L_C_H" {print $4}' "$SUMMARY")
BLL=$(awk '$1=="M_L_C_L" {print $4}' "$SUMMARY")

echo -e "B_Distance\t${BHH}\t${BHL}\t${BLH}\t${BLL}" \
>> "$WIDE"

