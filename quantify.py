# Libraries
import pandas as pd
import numpy as np
from scipy.sparse import coo_matrix, csr_matrix
from pathlib import Path
import subprocess


def bdg_to_df(
    bedgraph: str
) -> pd.DataFrame:

    # Data for COO matrix
    row, col, data = [], [], []
    barcodes, regions = {}, {}

    # Begin parsing
    for interval in bedgraph.split("\n"):

        if interval == "":
            continue

        chrom, start, stop, barcode, score, *_ = interval.split("\t")
        region = f"{chrom}:{start}-{stop}"

        # Populate hashmap
        if barcode not in barcodes:
            barcodes[barcode] = len(barcodes)

        if region not in regions:
            regions[region] = len(regions)

        row.append(barcodes[barcode])
        col.append(regions[region])
        data.append(score)

    sparse_mat = coo_matrix((data, (row, col)), dtype=np.float32)

    df = pd.DataFrame.sparse.from_spmatrix(sparse_mat)
    df.index = barcodes.keys()
    df.columns = regions.keys()

    return df


# Given a CSV file, where the index (SampleID) contains a path to BED files
def quantify(
    samples: str,
    binarize: bool = False
) -> csr_matrix:

    df = pd.read_csv(
        samples,
        index_col = 0,
        header = 0
    )

    unsorted_peaks = ""

    for sample, data in df.iterrows():

        _peaks, _ = subprocess.Popen(
            (
                "awk -F'\t' 'BEGIN { OFS = FS } { $4 = "
                f"\"{sample}\""
                "; print }' "
                f"{data['Peaks']}"
            ),
            shell = True, text = True,
            stdout = subprocess.PIPE
        ).communicate()

        unsorted_peaks += _peaks

    subprocess.Popen(
        (
            "bedtools sort | "
            "bedtools merge > "
            "merged.bed"
        ),
        shell = True, text = True,
        stdin = subprocess.PIPE,
        stdout = subprocess.PIPE
    ).communicate(input = unsorted_peaks)

    # Order of new bedgraph
    order = "$1\"\t\"$2\"\t\"$3\"\t\"$7\"\t\"$15"

    peaks, _ = subprocess.Popen(
        (
            "bedtools sort | "
            "bedtools intersect -wao -a merged.bed -b - |"

            # Convert the counts to a fraction of the coverage
            "awk 'BEGIN { FS=OFS = \"\t\" } "
            f"          {{  $15 = $14 / ($6 - $5)"
            f"          {'' if binarize else '* $8'};"
            f"           print {order} }}'"
        ),
        shell = True, text = True,
        stdin = subprocess.PIPE,
        stdout = subprocess.PIPE
    ).communicate(input = unsorted_peaks)

    df = bdg_to_df(peaks)
    return df


def save(
    df: pd.DataFrame,
    output: str
) -> None:

    df = df.T
    df.to_csv(Path(output) / "matrix.csv")


# Parse arguments
if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="comma-delimited file denoting path to peak BED")
    parser.add_argument("output", help="path to output folder to save sparse matrix")
    parser.add_argument("-b", "--binarize", help="binarize matrix counts, stores peaks availability instead of peak counts", action="store_true")
    args = parser.parse_args()

    df = quantify(args.input, binarize = args.binarize)
    save(df, args.output)