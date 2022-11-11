from typing import List

import pandas as pd
from binops import (calculate_bins_coverage, make_interval_centers,
                    make_linear_bins)
from datahandlers import DnaParts


def make_background_track(dna_parts: DnaParts, rnas: List, binsize: int = 1000) -> pd.DataFrame:
    bg_contacts = list()
    for rna_name in rnas:
        if rna_name in dna_parts.annotation:
            dna_df = dna_parts.read_dna_parts_single(rna_name)
            rna_chrom = dna_parts.annotation[rna_name]['chrom']
            dna_df = dna_df.query("chrom != @rna_chrom").reset_index(drop=True)
            bg_contacts.append(dna_df)
    bg_df = pd.concat(bg_contacts)
    bg_centers = make_interval_centers(bg_df)
    bg_bins = make_linear_bins(binsize, dna_parts.chromsizes)
    bg_track = calculate_bins_coverage(bg_bins, bg_centers)
    return bg_track
