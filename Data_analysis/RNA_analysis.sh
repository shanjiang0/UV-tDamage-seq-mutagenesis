#!/bin/bash
#SBATCH --job-name=RNA_RPKM
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --time=6:00:00
#SBATCH --mem=64G
#SBATCH --output=%x_%j.out
#SBATCH --mail-type=END,FAIL

set -euo pipefail


module purge
module load BEDTools/2.31.1-GCC-13.3.0

GROUP_DIR="/scratch/sj31280/06_projects/NHF1_UV_Target_Damage-seq_workflow2_sh/Final/02.Cluster/02.Target_4_group"
RNA="/scratch/sj31280/06_projects/NHF1_UV_Target_Damage-seq_workflow2_sh/Final/03.Multi_Analysis/04.RNA/01.Input_data/HFFc6_rep1_hg19.bed"
OUTDIR="/scratch/sj31280/06_projects/NHF1_UV_Target_Damage-seq_workflow2_sh/Final/03.Multi_Analysis/04.RNA/02.Result"

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


if [ ! -f "$RNA" ]; then
    echo "ERROR: RNA file not found:"
    echo "$RNA"
    exit 1
fi

for GROUP in "${GROUP_ORDER[@]}"; do

    if [ ! -f "${GROUP_FILES[$GROUP]}" ]; then
        echo "ERROR: Group BED file not found:"
        echo "${GROUP_FILES[$GROUP]}"
        exit 1
    fi

done


for GROUP in "${GROUP_ORDER[@]}"; do

    CLUSTER="${GROUP_FILES[$GROUP]}"
    INTERSECT_OUT="${OUTDIR}/${GROUP}.RNA_overlap.bed"
    RESULT_OUT="${OUTDIR}/${GROUP}.RNA_RPKM.tsv"


    bedtools intersect \
        -a "$CLUSTER" \
        -b "$RNA" \
        -wao \
        > "$INTERSECT_OUT"


    awk '
    BEGIN {
        OFS="\t"

        print "Chr","Start","End","Gene",
              "ClusterLength",
              "Weighted_RNA_Sum",
              "Covered_bp",
              "RNA_RPM",
              "RNA_RPKM"
    }

    {
        chr=$1
        start=$2
        end=$3
        gene=$4

        key=chr OFS start OFS end OFS gene

        cluster_len=end-start

        # RNA signal is B column 4 = overall field 10
        signal=$10

        overlap=$NF

        clen[key]=cluster_len


        if (overlap > 0 && signal != ".") {
            weighted[key] += signal * overlap
            covered[key] += overlap
        }

        seen[key]=1
    }

    END {

        for (key in seen) {
            if (clen[key] > 0) {
                rpm = weighted[key] / clen[key]
                rpkm = rpm * 1000

            } else {

                rpm=0
                rpkm=0
            }


            printf "%s\t%d\t%.10f\t%d\t%.10f\t%.10f\n",
                   key,
                   clen[key],
                   weighted[key],
                   covered[key],
                   rpm,
                   rpkm
        }
    }
    ' "$INTERSECT_OUT" \
    | {
        read header
        echo "$header"
        sort -k1,1V -k2,2n -k3,3n
    } \
    > "$RESULT_OUT"


    echo "Created:"
    echo "$RESULT_OUT"

done


SUMMARY="${OUTDIR}/RNA_RPKM_4groups_summary.tsv"

echo -e "Group\tCluster_Number\tMean_RPKM\tMean_RPM\tNonzero_Clusters" \
> "$SUMMARY"


for GROUP in "${GROUP_ORDER[@]}"; do

    FILE="${OUTDIR}/${GROUP}.RNA_RPKM.tsv"

    awk -v group="$GROUP" '
    BEGIN {
        OFS="\t"
    }

    NR > 1 {

        n++

        sum_rpm += $8
        sum_rpkm += $9

        if ($9 > 0) {
            nonzero++
        }
    }

    END {

        if (n > 0) {

            printf "%s\t%d\t%.10f\t%.10f\t%d\n",
                   group,
                   n,
                   sum_rpkm/n,
                   sum_rpm/n,
                   nonzero

        } else {

            print group,0,0,0,0
        }
    }
    ' "$FILE" >> "$SUMMARY"

done


WIDE="${OUTDIR}/RNA_mean_RPKM_4groups.tsv"

echo -e "Feature\tM_H_C_H\tM_H_C_L\tM_L_C_H\tM_L_C_L" \
> "$WIDE"


MHCH=$(awk '$1=="M_H_C_H" {print $3}' "$SUMMARY")
MHCL=$(awk '$1=="M_H_C_L" {print $3}' "$SUMMARY")
MLCH=$(awk '$1=="M_L_C_H" {print $3}' "$SUMMARY")
MLCL=$(awk '$1=="M_L_C_L" {print $3}' "$SUMMARY")


echo -e "RNA_RPKM\t${MHCH}\t${MHCL}\t${MLCH}\t${MLCL}" \
>> "$WIDE"

