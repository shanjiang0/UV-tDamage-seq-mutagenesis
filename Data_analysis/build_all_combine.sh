#!/bin/bash
#SBATCH --job-name=Combine_Model
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --time=10:00:00
#SBATCH --mem=32G
#SBATCH --output=%x_%j.out
#SBATCH --mail-type=END,FAIL

set -euo pipefail


module purge
module load Python/3.14.2-GCCcore-15.2.0


WORKDIR="/scratch/sj31280/06_projects/NHF1_UV_Target_Damage-seq_workflow2_sh/Final/03.Multi_Analysis/07.Model"

mkdir -p "$WORKDIR"
cd "$WORKDIR"


TARGET_DIR="/scratch/sj31280/06_projects/NHF1_UV_Target_Damage-seq_workflow2_sh/Final/02.Cluster/02.Target_4_group"
XR_DIR="/scratch/sj31280/06_projects/NHF1_UV_Target_Damage-seq_workflow2_sh/Final/03.Multi_Analysis/03.XR/02.Result"
RNA_DIR="/scratch/sj31280/06_projects/NHF1_UV_Target_Damage-seq_workflow2_sh/Final/03.Multi_Analysis/04.RNA/02.Result"
HISTONE_DIR="/scratch/sj31280/06_projects/NHF1_UV_Target_Damage-seq_workflow2_sh/Final/03.Multi_Analysis/05.Histone/02.Result"
CHROM3D_DIR="/scratch/sj31280/06_projects/NHF1_UV_Target_Damage-seq_workflow2_sh/Final/03.Multi_Analysis/06.Chrom3D/02.Result"


python - <<'PYTHON'
import os
import csv
import sys


WORKDIR = "/scratch/sj31280/06_projects/NHF1_UV_Target_Damage-seq_workflow2_sh/Final/03.Multi_Analysis/07.Model"
TARGET_DIR = "/scratch/sj31280/06_projects/NHF1_UV_Target_Damage-seq_workflow2_sh/Final/02.Cluster/02.Target_4_group"
XR_DIR = "/scratch/sj31280/06_projects/NHF1_UV_Target_Damage-seq_workflow2_sh/Final/03.Multi_Analysis/03.XR/02.Result"
RNA_DIR = "/scratch/sj31280/06_projects/NHF1_UV_Target_Damage-seq_workflow2_sh/Final/03.Multi_Analysis/04.RNA/02.Result"
HISTONE_DIR = "/scratch/sj31280/06_projects/NHF1_UV_Target_Damage-seq_workflow2_sh/Final/03.Multi_Analysis/05.Histone/02.Result"
CHROM3D_DIR = "/scratch/sj31280/06_projects/NHF1_UV_Target_Damage-seq_workflow2_sh/Final/03.Multi_Analysis/06.Chrom3D/02.Result"

# ============================================================
# GROUP MAPPING
# ============================================================

groups = {
    "C_H_M_H": {
        "target": "C_H_M_H.bed",
        "old": "M_H_C_H"
    },

    "C_H_M_L": {
        "target": "C_H_M_L.bed",
        "old": "M_L_C_H"
    },

    "C_L_M_H": {
        "target": "C_L_M_H.bed",
        "old": "M_H_C_L"
    },

    "C_L_M_L": {
        "target": "C_L_M_L.bed",
        "old": "M_L_C_L"
    }
}


# ============================================================
# XR TIMEPOINTS
# ============================================================

xr_timepoints = [
    ("CPD12m", "XR_12m"),
    ("CPD1h",  "XR_1h"),
    ("CPD4h",  "XR_4h"),
    ("CPD8h",  "XR_8h"),
    ("CPD16h", "XR_16h"),
    ("CPD24h", "XR_24h"),
    ("CPD48h", "XR_48h")
]


# ============================================================
# FINAL COLUMN ORDER
# ============================================================

final_header = [
    "chr",
    "start",
    "end",
    "name",
    "cluster_length",
    "CPD_Enrichment",
    "MELA_CPKM",
    "z_CPD",
    "z_MELA",

    "XR_12m",
    "XR_1h",
    "XR_4h",
    "XR_8h",
    "XR_16h",
    "XR_24h",
    "XR_48h",

    "DNase",
    "H3K27ac",
    "H3K4me1",
    "H3K4me3",
    "H3K27me3",
    "H3K9me3",
    "CTCF",

    "RNA_RPKM",

    "Dist_A",
    "Dist_B",

    "Group"
]


def key_from_values(chrom, start, end, gene):
    return (
        str(chrom),
        str(start),
        str(end),
        str(gene)
    )


def require_file(path):

    if not os.path.isfile(path):
        print(f"ERROR: file not found:\n{path}", file=sys.stderr)
        sys.exit(1)

    if os.path.getsize(path) == 0:
        print(f"ERROR: file is empty:\n{path}", file=sys.stderr)
        sys.exit(1)



def read_target(path):

    require_file(path)
    rows = []
    with open(path) as f:
    
        reader = csv.reader(f, delimiter="\t")
        first = next(reader)


        if first[0].lower() in ("chr", "chrom", "chromosome"):
            data_iter = reader
        else:
            data_iter = iter([first] + list(reader))


        for r in data_iter:

            if len(r) < 10:
                continue

            rows.append({
                "chr": r[0],
                "start": r[1],
                "end": r[2],
                "name": r[3],
                "cluster_length": r[4],
                "CPD_Enrichment": r[5],
                "MELA_CPKM": r[6],
                "z_CPD": r[7],
                "z_MELA": r[8],
                "Group": r[9]
            })

    return rows


def read_tsv_dict(path):
    require_file(path)
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        return list(reader)

# ============================================================
# LOAD XR RPKM
# ============================================================

def load_xr(path):
    records = read_tsv_dict(path)
    values = {}
    for r in records:

        key = key_from_values(
            r["Chr"],
            r["Start"],
            r["End"],
            r["Gene"]
        )
        values[key] = r["RPKM"]
    return values


# ============================================================
# LOAD RNA RPKM
# ============================================================

def load_rna(path):
    records = read_tsv_dict(path)
    values = {}
    for r in records:
        key = key_from_values(
            r["Chr"],
            r["Start"],
            r["End"],
            r["Gene"]
        )
        values[key] = r["RNA_RPKM"]

    return values


# ============================================================
# LOAD HISTONE
# ============================================================

def load_histone(path):

    records = read_tsv_dict(path)

    values = {}

    wanted = [
        "DNase",
        "H3K27ac",
        "H3K4me1",
        "H3K4me3",
        "H3K27me3",
        "H3K9me3",
        "CTCF",
    ]

    for r in records:

        key = key_from_values(
            r["Chr"],
            r["Start"],
            r["End"],
            r["Gene"]
        )

        values[key] = {
            feature: r[feature]
            for feature in wanted
        }

    return values


# ============================================================
# LOAD CHROM3D A OR B DISTANCE
# ============================================================

def load_distance(path, distance_column):

    records = read_tsv_dict(path)

    values = {}

    for r in records:

        key = key_from_values(
            r["Chr"],
            r["Start"],
            r["End"],
            r["Gene"]
        )

        values[key] = r[distance_column]

    return values


# ============================================================
# PROCESS EACH GROUP
# ============================================================

all_rows = []

for new_group, info in groups.items():

    old_group = info["old"]

    target_file = os.path.join(
        TARGET_DIR,
        info["target"]
    )

    rows = read_target(target_file)

    print(f"Starting clusters: {len(rows)}")


    # --------------------------------------------------------
    # XR
    # --------------------------------------------------------

    xr_maps = {}

    for prefix, output_name in xr_timepoints:

        path = os.path.join(
            XR_DIR,
            f"{prefix}_{old_group}.RPKM.tsv"
        )

        xr_maps[output_name] = load_xr(path)


    # --------------------------------------------------------
    # RNA
    # --------------------------------------------------------

    rna_file = os.path.join(
        RNA_DIR,
        f"{old_group}.RNA_RPKM.tsv"
    )

    rna_map = load_rna(rna_file)


    # --------------------------------------------------------
    # Histone / DNase
    # --------------------------------------------------------

    histone_file = os.path.join(
        HISTONE_DIR,
        f"{old_group}.All_Histone_signal.tsv"
    )

    histone_map = load_histone(histone_file)


    # --------------------------------------------------------
    # Chrom3D A
    # --------------------------------------------------------

    dist_a_file = os.path.join(
        CHROM3D_DIR,
        f"{old_group}.A_distance.tsv"
    )

    dist_a_map = load_distance(
        dist_a_file,
        "A_Distance"
    )


    # --------------------------------------------------------
    # Chrom3D B
    # --------------------------------------------------------

    dist_b_file = os.path.join(
        CHROM3D_DIR,
        f"{old_group}.B_distance.tsv"
    )

    dist_b_map = load_distance(
        dist_b_file,
        "B_Distance"
    )


    # --------------------------------------------------------
    # Append everything by exact coordinate + gene
    # --------------------------------------------------------

    missing_counts = {
        "XR_12m": 0,
        "XR_1h": 0,
        "XR_4h": 0,
        "XR_8h": 0,
        "XR_16h": 0,
        "XR_24h": 0,
        "XR_48h": 0,
        "Histone": 0,
        "RNA": 0,
        "Dist_A": 0,
        "Dist_B": 0
    }


    combined_rows = []


    for r in rows:

        key = key_from_values(
            r["chr"],
            r["start"],
            r["end"],
            r["name"]
        )


        out = dict(r)


        # --------------------------------------------
        # XR
        # --------------------------------------------

        for output_name in [
            "XR_12m",
            "XR_1h",
            "XR_4h",
            "XR_8h",
            "XR_16h",
            "XR_24h",
            "XR_48h"
        ]:

            value = xr_maps[output_name].get(key)

            if value is None:
                value = "NA"
                missing_counts[output_name] += 1

            out[output_name] = value


        # --------------------------------------------
        # Histone / DNase
        # --------------------------------------------

        hist = histone_map.get(key)

        if hist is None:

            missing_counts["Histone"] += 1

            out["DNase"] = "NA"
            out["H3K27ac"] = "NA"
            out["H3K4me1"] = "NA"
            out["H3K4me3"] = "NA"
            out["H3K27me3"] = "NA"
            out["H3K9me3"] = "NA"
            out["CTCF"] = "NA"

        else:

            out.update(hist)


        # --------------------------------------------
        # RNA
        # --------------------------------------------

        value = rna_map.get(key)

        if value is None:
            value = "NA"
            missing_counts["RNA"] += 1

        out["RNA_RPKM"] = value


        # --------------------------------------------
        # Chrom3D A
        # --------------------------------------------

        value = dist_a_map.get(key)

        if value is None:
            value = "NA"
            missing_counts["Dist_A"] += 1

        out["Dist_A"] = value


        # --------------------------------------------
        # Chrom3D B
        # --------------------------------------------

        value = dist_b_map.get(key)

        if value is None:
            value = "NA"
            missing_counts["Dist_B"] += 1

        out["Dist_B"] = value


        # --------------------------------------------
        # Force final group name to new naming system
        # --------------------------------------------

        out["Group"] = new_group

        combined_rows.append(out)


    # --------------------------------------------------------
    # Write group-specific combined table
    # --------------------------------------------------------

    group_output = os.path.join(
        WORKDIR,
        f"{new_group}.Combined.bed"
    )


    with open(group_output, "w", newline="") as f:

        writer = csv.DictWriter(
            f,
            delimiter="\t",
            fieldnames=final_header,
            lineterminator="\n"
        )

        writer.writeheader()

        for r in combined_rows:
            writer.writerow(r)


    print(f"Created: {group_output}")

    print("Missing-value check:")

    for feature, count in missing_counts.items():
        print(f"  {feature}: {count}")


    all_rows.extend(combined_rows)


# ============================================================
# SORT ALL FOUR GROUPS TOGETHER
# ============================================================

def chromosome_sort_key(chrom):

    c = chrom

    if c.startswith("chr"):
        c = c[3:]

    if c.isdigit():
        return (0, int(c))

    order = {
        "X": 23,
        "Y": 24,
        "M": 25,
        "MT": 25
    }

    if c in order:
        return (0, order[c])

    return (1, c)


all_rows.sort(
    key=lambda r: (
        chromosome_sort_key(r["chr"]),
        int(r["start"]),
        int(r["end"]),
        r["name"]
    )
)


# ============================================================
# WRITE FINAL ALL-COMBINED FILE
# ============================================================

final_output = os.path.join(
    WORKDIR,
    "ALL_Combined.bed"
)


with open(final_output, "w", newline="") as f:

    writer = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=final_header,
        lineterminator="\n"
    )

    writer.writeheader()

    for r in all_rows:
        writer.writerow(r)


print()
print("============================================================")
print("DONE")
print("============================================================")

print(f"Total clusters: {len(all_rows)}")

print()
print("Final file:")
print(final_output)

PYTHON


# ============================================================
# CHECK RESULT
# ============================================================

echo
echo "============================================================"
echo "FINAL OUTPUT"
echo "============================================================"

head -n 5 "${WORKDIR}/ALL_Combined.bed"

echo
echo "Number of rows:"
wc -l "${WORKDIR}/ALL_Combined.bed"

echo
echo "Group counts:"
awk -F '\t' '
NR > 1 {
    count[$NF]++
}
END {
    for (g in count) {
        print g, count[g]
    }
}' "${WORKDIR}/ALL_Combined.bed"
