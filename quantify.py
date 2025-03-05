# Libraries
import pandas as pd
import numpy as np
from scipy.sparse import coo_matrix, csr_matrix
from pathlib import Path
import subprocess



def bdg_to_df(
    bedgraph: str,
    binarize: bool = False
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
        data.append(1 if binarize else score)

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
    order = "$1\"\t\"$2\"\t\"$3\"\t\"$7\"\t\"$8"

    peaks, _ = subprocess.Popen(
        (
            "bedtools sort | "
            "bedtools intersect -wao -a merged.bed -b - |"

            # Convert the counts to a fraction of the coverage
            "awk 'BEGIN { FS=OFS = \"\t\" } "
            f"         {{ print {order} }}'"
        ),
        shell = True, text = True,
        stdin = subprocess.PIPE,
        stdout = subprocess.PIPE
    ).communicate(input = unsorted_peaks)

    df = bdg_to_df(peaks, binarize = binarize)
    return df


def save(
    df: pd.DataFrame,
    output: str
) -> None:

    df = df.T
    df.to_csv(Path(output) / "matrix.csv")


################################################################################
# Run your functions here!

# df = quantify("/path/to/sample_sheet.csv", binarize = False)
# save(df, "/path/to/output_folder")