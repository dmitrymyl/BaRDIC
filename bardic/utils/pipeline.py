from typing import Dict, List

import pandas as pd

from ..api.convert import annotation_to_dict
from ..api.io import load_bedgraph
from . import (bed2h5, calculate_bin_sizes, calculate_scaling_splines,
               dnadataset_to_rdc, estimate_significance, fetch_peaks,
               format_peaks, make_background_track)


def run_pipeline(dna_parts_fname: str,
                 dna_dataset_fname: str,
                 rdc_fname: str,
                 chromsizes: Dict[str, int],
                 annotation: pd.DataFrame,
                 selection_results_fname: str,
                 bg_fname: str,
                 rna_list: List,
                 bg_binsize: int,
                 peaks_threshold: float,
                 peaks_output: str,
                 binsize_params: Dict = {},
                 rdc_params: Dict = {},
                 scaling_params: Dict = {},
                 peaks_format_params: Dict = {},
                 makebg: bool = True,
                 n_cores: int = 1):
    chromdict = chromsizes
    annotation_dict = annotation_to_dict(annotation)

    dna_dataset = bed2h5(dna_parts_fname, dna_dataset_fname, chromdict, annotation_dict)

    selection_df = calculate_bin_sizes(dna_dataset, n_cores=n_cores, **binsize_params)
    selection_df.to_csv(selection_results_fname, sep='\t', header=True, index=False)

    if makebg:
        bg_track = make_background_track(dna_dataset, rna_list, bg_binsize)
        bg_track.to_csv(bg_fname, header=False, index=False, sep='\t')
    else:
        bg_track = load_bedgraph(bg_fname)

    rdc_data = dnadataset_to_rdc(dna_dataset, bg_track, rdc_fname, n_cores=n_cores, **rdc_params)

    calculate_scaling_splines(rdc_data, n_cores=n_cores, **scaling_params)

    estimate_significance(rdc_data, n_cores)
    peaks = fetch_peaks(rdc_data, peaks_threshold, n_cores)
    formatted_peaks = format_peaks(peaks, **peaks_format_params)
    formatted_peaks.to_csv(peaks_output, sep='\t', header=False, index=False)
