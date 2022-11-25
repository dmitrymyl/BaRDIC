from typing import List

import pandas as pd

from ..api.binops import (calculate_bins_coverage, make_interval_centers,
                          make_linear_bins)
from ..api.formats import DnaDataset


def make_background_track(dna_dataset: DnaDataset, rnas: List, binsize: int = 1000) -> pd.DataFrame:
    bg_contacts = list()
    for rna_name in rnas:
        if rna_name in dna_dataset.annotation:
            dna_df = dna_dataset.read_dna_parts_single(rna_name)
            rna_chrom = dna_dataset.annotation[rna_name].chrom
            dna_df = dna_df.query("chrom != @rna_chrom").reset_index(drop=True)
            bg_contacts.append(dna_df)
    bg_df = pd.concat(bg_contacts)
    bg_centers = make_interval_centers(bg_df)
    bg_bins = make_linear_bins(binsize, dna_dataset.chromsizes)
    bg_track = calculate_bins_coverage(bg_bins, bg_centers)
    return bg_track
