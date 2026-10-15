#!/bin/bash
#SBATCH --job-name=XR_RPKM
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --output=%x_%j.out
#SBATCH --mail-type=END,FAIL
set -euo pipefail


# ============================================================
# MODULES
# ============================================================

module purge
module load BEDTools/2.31.1-GCC-13.3.0


# ============================================================
# DIRECTORIES
# ============================================================

INPUT_DIR="/scratch/sj31280/06_projects/NHF1_UV_Target_Damage-seq_workflow2_sh/Final/03.Multi_Analysis/03.XR/01.Input_data"
GROUP_DIR="/scratch/sj31280/06_projects/NHF1_UV_Target_Damage-seq_workflow2_sh/Final/02.Cluster/02.Target_4_group"
OUTDIR="/scratch/sj31280/06_projects/NHF1_UV_Target_Damage-seq_workflow2_sh/Final/03.Multi_Analysis/03.XR/02.Result"

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


declare -A XR_FILES

XR_FILES[CPD12m]="${INPUT_DIR}/CPD12m_merged_replicates.bed"
XR_FILES[CPD1h]="${INPUT_DIR}/CPD1h_merged_replicates.bed"
XR_FILES[CPD4h]="${INPUT_DIR}/CPD4h_merged_replicates.bed"
XR_FILES[CPD8h]="${INPUT_DIR}/CPD8h_merged_replicates.bed"
XR_FILES[CPD16h]="${INPUT_DIR}/CPD16h_merged_replicates.bed"
XR_FILES[CPD24h]="${INPUT_DIR}/CPD24h_merged_replicates.bed"
XR_FILES[CPD48h]="${INPUT_DIR}/CPD48h_merged_replicates.bed"

TIMEPOINT_ORDER=(
    CPD12m
    CPD1h
    CPD4h
    CPD8h
    CPD16h
    CPD24h
    CPD48h
)


# ============================================================
# CHECK INPUT FILES
# ============================================================

echo "Checking input files..."

for GROUP in "${GROUP_ORDER[@]}"; do

    if [ ! -f "${GROUP_FILES[$GROUP]}" ]; then
        echo "ERROR: Missing group file:"
        echo "${GROUP_FILES[$GROUP]}"
        exit 1
    fi

done


for TIMEPOINT in "${TIMEPOINT_ORDER[@]}"; do

    if [ ! -f "${XR_FILES[$TIMEPOINT]}" ]; then
        echo "ERROR: Missing XR file:"
        echo "${XR_FILES[$TIMEPOINT]}"
        exit 1
    fi

done



LIBRARY_SUMMARY="${OUTDIR}/XR_total_mapped_reads.tsv"
echo -e "Timepoint\tTotalReads" > "$LIBRARY_SUMMARY"



for TIMEPOINT in "${TIMEPOINT_ORDER[@]}"; do

    XR_FILE="${XR_FILES[$TIMEPOINT]}"



    TOTAL_READS=$(wc -l < "$XR_FILE")

    echo "Total mapped reads = $TOTAL_READS"

    echo -e "${TIMEPOINT}\t${TOTAL_READS}" \
        >> "$LIBRARY_SUMMARY"



    for GROUP in "${GROUP_ORDER[@]}"; do

        CLUSTER="${GROUP_FILES[$GROUP]}"

        PREFIX="${OUTDIR}/${TIMEPOINT}_${GROUP}"


        bedtools intersect \
            -a "$CLUSTER" \
            -b "$XR_FILE" \
            -c \
            > "${PREFIX}.XR_counts.bed"


        awk -v total="$TOTAL_READS" '
        BEGIN {
            OFS="\t"

            print "Chr","Start","End","Gene",
                  "RawRead","ClusterLength",
                  "TotalReads","RPKM"
        }

        {
            cluster_len=$3-$2
            raw=$NF

            if (total > 0 && cluster_len > 0) {

                rpkm=(raw * 1000000000) / \
                     (total * cluster_len)

            } else {

                rpkm=0

            }

            printf "%s\t%d\t%d\t%s\t%d\t%d\t%d\t%.8f\n",
                   $1,$2,$3,$4,
                   raw,
                   cluster_len,
                   total,
                   rpkm
        }
        ' "${PREFIX}.XR_counts.bed" \
        > "${PREFIX}.RPKM.tsv"



        awk '
        BEGIN {
            OFS="\t"
        }

        NR > 1 {
            n++
            x[n]=$8
            sum += $8

            if ($5 > 0) {
                nonzero++
            }
        }

        END {

            if (n == 0) {
                print "0","0","0","0"
                exit
            }

            mean=sum/n

            # simple insertion sort for median
            for (i=2; i<=n; i++) {
                v=x[i]
                j=i-1

                while (j>=1 && x[j]>v) {
                    x[j+1]=x[j]
                    j--
                }

                x[j+1]=v
            }

            if (n % 2 == 1) {
                median=x[(n+1)/2]
            } else {
                median=(x[n/2] + x[n/2+1])/2
            }

            printf "%d\t%.8f\t%.8f\t%d\n",
                   n,mean,median,nonzero
        }
        ' "${PREFIX}.RPKM.tsv" \
        > "${PREFIX}.summary.tmp"

    done

done



SUMMARY="${OUTDIR}/XR_RPKM_group_summary.tsv"

echo -e "Timepoint\tGroup\tCluster_Number\tMean_RPKM\tMedian_RPKM\tNonzero_Clusters" \
> "$SUMMARY"


for TIMEPOINT in "${TIMEPOINT_ORDER[@]}"; do

    for GROUP in "${GROUP_ORDER[@]}"; do

        TMP="${OUTDIR}/${TIMEPOINT}_${GROUP}.summary.tmp"

        read N MEAN MEDIAN NONZERO < "$TMP"

        echo -e "${TIMEPOINT}\t${GROUP}\t${N}\t${MEAN}\t${MEDIAN}\t${NONZERO}" \
        >> "$SUMMARY"

        rm "$TMP"

    done

done


WIDE="${OUTDIR}/XR_mean_RPKM_4groups.tsv"

echo -e "Timepoint\tM_H_C_H\tM_H_C_L\tM_L_C_H\tM_L_C_L" \
> "$WIDE"


for TIMEPOINT in "${TIMEPOINT_ORDER[@]}"; do

    A=$(awk -v t="$TIMEPOINT" '$1==t && $2=="M_H_C_H" {print $4}' "$SUMMARY")
    B=$(awk -v t="$TIMEPOINT" '$1==t && $2=="M_H_C_L" {print $4}' "$SUMMARY")
    C=$(awk -v t="$TIMEPOINT" '$1==t && $2=="M_L_C_H" {print $4}' "$SUMMARY")
    D=$(awk -v t="$TIMEPOINT" '$1==t && $2=="M_L_C_L" {print $4}' "$SUMMARY")

    echo -e "${TIMEPOINT}\t${A}\t${B}\t${C}\t${D}" \
        >> "$WIDE"

done
