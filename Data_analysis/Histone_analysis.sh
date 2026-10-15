#!/bin/bash
#SBATCH --job-name=Histone_Cluster
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --time=100:00:00
#SBATCH --mem=256G
#SBATCH --output=%x_%j.out
#SBATCH --mail-type=END,FAIL

set -euo pipefail

module purge
module load BEDTools/2.31.1-GCC-13.3.0


INPUT_DIR="/scratch/sj31280/06_projects/NHF1_UV_Target_Damage-seq_workflow2_sh/Final/03.Multi_Analysis/05.Histone/01.Input_data"
GROUP_DIR="/scratch/sj31280/06_projects/NHF1_UV_Target_Damage-seq_workflow2_sh/Final/02.Cluster/02.Target_4_group"
OUTDIR="/scratch/sj31280/06_projects/NHF1_UV_Target_Damage-seq_workflow2_sh/Final/03.Multi_Analysis/05.Histone/02.Result"

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


FEATURES=(
    CTCF
    DNase
    H3K27ac
    H3K4me1
    H3K4me3
    H3K27me3
    H3K9me3
)


for FEATURE in "${FEATURES[@]}"; do

    TRACK="${INPUT_DIR}/${FEATURE}_hg19.bedGraph"
    if [ ! -s "$TRACK" ]; then
        echo "ERROR: missing or empty signal file:"
        echo "$TRACK"
        exit 1
    fi

done

for GROUP in "${GROUP_ORDER[@]}"; do

    if [ ! -f "${GROUP_FILES[$GROUP]}" ]; then
        echo "ERROR: missing group file:"
        echo "${GROUP_FILES[$GROUP]}"
        exit 1
    fi

done



for FEATURE in "${FEATURES[@]}"; do

    TRACK="${INPUT_DIR}/${FEATURE}_hg19.bedGraph"


    for GROUP in "${GROUP_ORDER[@]}"; do

        CLUSTER="${GROUP_FILES[$GROUP]}"
        TMP="${OUTDIR}/${GROUP}.${FEATURE}.overlap.tmp"
        RESULT="${OUTDIR}/${GROUP}.${FEATURE}.signal.tsv"
        echo "Group: $GROUP"


        bedtools intersect \
            -a "$CLUSTER" \
            -b "$TRACK" \
            -wao \
            > "$TMP"


        awk -v feature="$FEATURE" '
        BEGIN {
            OFS="\t"

            print "Chr",
                  "Start",
                  "End",
                  "Gene",
                  "ClusterLength",
                  "Covered_bp",
                  "WeightedSum",
                  feature
        }

        {
            chr=$1
            start=$2
            end=$3
            gene=$4

            key=chr OFS start OFS end OFS gene

            cluster_len=end-start

            signal=$10
            overlap=$11

            clen[key]=cluster_len
            seen[key]=1

            if (overlap > 0 && signal != "." && signal != "") {

                weighted[key] += overlap * signal
                covered[key] += overlap

            }
        }

        END {

            for (key in seen) {

                if (clen[key] > 0) {
                    avg_signal=weighted[key]/clen[key]
                } else {
                    avg_signal=0
                }

                printf "%s\t%d\t%d\t%.10f\t%.10f\n",
                       key,
                       clen[key],
                       covered[key],
                       weighted[key],
                       avg_signal
            }
        }
        ' "$TMP" \
        | {
            read header
            echo "$header"
            sort -k1,1V -k2,2n -k3,3n
        } \
        > "$RESULT"

        rm -f "$TMP"

    done

done



for GROUP in "${GROUP_ORDER[@]}"; do


    FINAL="${OUTDIR}/${GROUP}.All_Histone_signal.tsv"

    paste \
        "${OUTDIR}/${GROUP}.CTCF.signal.tsv" \
        "${OUTDIR}/${GROUP}.DNase.signal.tsv" \
        "${OUTDIR}/${GROUP}.H3K27ac.signal.tsv" \
        "${OUTDIR}/${GROUP}.H3K4me1.signal.tsv" \
        "${OUTDIR}/${GROUP}.H3K4me3.signal.tsv" \
        "${OUTDIR}/${GROUP}.H3K27me3.signal.tsv" \
        "${OUTDIR}/${GROUP}.H3K9me3.signal.tsv" \
    | awk '
    BEGIN {
        OFS="\t"
    }

    NR==1 {
        print "Chr",
              "Start",
              "End",
              "Gene",
              "ClusterLength",
              "CTCF",
              "DNase",
              "H3K27ac",
              "H3K4me1",
              "H3K4me3",
              "H3K27me3",
              "H3K9me3"
        next
    }

    {
        print $1,
              $2,
              $3,
              $4,
              $5,
              $8,
              $16,
              $24,
              $32,
              $40,
              $48,
              $56
    }
    ' \
    > "$FINAL"

    echo "Created:"
    echo "$FINAL"

done

echo
echo "Final combined files:"

for GROUP in "${GROUP_ORDER[@]}"; do
    echo "${OUTDIR}/${GROUP}.All_Histone_signal.tsv"
done
