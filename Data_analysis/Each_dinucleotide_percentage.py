from collections import Counter
import os
import re

# ============================================================
# DIRECTORIES
# ============================================================

WORKDIR = "/scratch/sj31280/06_projects/NHF1_UV_Target_Damage-seq_workflow2_sh/Final/03.Multi_Analysis/02.Dinu_seq"
FASTA_DIR = os.path.join(WORKDIR, "02.Cluster_fasta")
OUTPUT_DIR = os.path.join(WORKDIR, "03.Dinucleotide_results")
os.makedirs(OUTPUT_DIR, exist_ok=True)


groups = [
    "M_H_C_H",
    "M_H_C_L",
    "M_L_C_H",
    "M_L_C_L"
]

dinucleotides = [
    "TT", "TC", "TG", "TA",
    "CT", "CC", "CG", "CA",
    "GT", "GC", "GG", "GA",
    "AT", "AC", "AG", "AA"
]

def read_fasta(fasta_file):

    header = None
    sequence = []
    with open(fasta_file, "r") as f:

        for line in f:
            line = line.strip()
            if line.startswith(">"):

                if header is not None:
                    yield header, "".join(sequence)
                header = line[1:]
                sequence = []
            else:
                sequence.append(line)

        if header is not None:
            yield header, "".join(sequence)

group_results = {}


for group in groups:

    fasta_file = os.path.join(
        FASTA_DIR,
        f"{group}.fasta"
    )

    per_cluster_output = os.path.join(
        OUTPUT_DIR,
        f"{group}.dinucleotide_percentage_per_cluster.tsv"
    )


    group_percentage_sum = {
        d: 0.0 for d in dinucleotides
    }

    cluster_count = 0


    with open(per_cluster_output, "w") as out:

        out.write(
            "Chr\tStart\tEnd\tGene\t" +
            "\t".join(dinucleotides) +
            "\n"
        )

        for header, sequence in read_fasta(fasta_file):
            sequence = sequence.upper()


            match = re.search(
                r"([^:]+)::(chr[0-9XYM]+):(\d+)-(\d+)",
                header
            )

            if not match:
                print(f"WARNING: cannot parse header: {header}")
                continue


            gene, chrom, start, end = match.groups()


            counts = Counter()

            for i in range(len(sequence) - 1):

                dinuc = sequence[i:i+2]

                if len(dinuc) == 2 and all(
                    base in "ACGT" for base in dinuc
                ):
                    counts[dinuc] += 1


            total_dinucleotides = sum(
                counts.get(d, 0)
                for d in dinucleotides
            )


            if total_dinucleotides == 0:
                continue


            percentages = {}

            for d in dinucleotides:

                percentages[d] = (
                    counts.get(d, 0)
                    / total_dinucleotides
                    * 100
                )


            for d in dinucleotides:

                group_percentage_sum[d] += percentages[d]

            cluster_count += 1



            out.write(
                f"{chrom}\t{start}\t{end}\t{gene}\t" +
                "\t".join(
                    f"{percentages[d]:.6f}"
                    for d in dinucleotides
                ) +
                "\n"
            )

    group_average = {}

    for d in dinucleotides:

        if cluster_count > 0:

            group_average[d] = (
                group_percentage_sum[d]
                / cluster_count
            )

        else:

            group_average[d] = 0.0


    group_results[group] = group_average


final_output = os.path.join(
    OUTPUT_DIR,
    "Dinucleotide_percentage_4groups.tsv"
)


with open(final_output, "w") as out:

    out.write(
        "Dinucleotide\t" +
        "\t".join(groups) +
        "\n"
    )

    for d in dinucleotides:

        out.write(
            d + "\t" +
            "\t".join(
                f"{group_results[g][d]:.6f}"
                for g in groups
            ) +
            "\n"
        )
