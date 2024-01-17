from typing import Dict, List

import pandas as pd

from ..api.convert import annotation_to_dict
from ..api.io import read_bedgraph
from .background import make_background_track
from .binsizes import optimize_bin_sizes
from .dnadataset import bed2h5
from .peaks import estimate_significance, fetch_peaks, format_peaks
from .rdc import dnadataset_to_rdc
from .scaling import calculate_scaling_splines


def run_pipeline(dna_parts_fname: str,
                 dna_dataset_fname: str,
                 rdc_fname: str,
                 chromsizes: Dict[str, int],
                 annotation: pd.DataFrame,
                 selection_results_fname: str,
                 bg_fname: str,
                 rna_list: List,
                 bg_binsize: int,
                 qval_threshold: float,
                 qval_type: str,
                 peaks_output: str,
                 binsize_params: Dict = {},
                 rdc_params: Dict = {},
                 scaling_params: Dict = {},
                 peaks_format_params: Dict = {},
                 make_bg: bool = True,
                 uniform_bg: bool = False,
                 n_cores: int = 1):
    chromdict = chromsizes
    annotation_dict = annotation_to_dict(annotation)

    dna_dataset = bed2h5(dna_parts_fname, dna_dataset_fname, chromdict, annotation_dict)

    selection_df = optimize_bin_sizes(dna_dataset, n_cores=n_cores, **binsize_params)
    selection_df.to_csv(selection_results_fname, sep='\t', header=True, index=False)

    if make_bg:
        bg_track = make_background_track(dna_dataset, rna_list, bg_binsize, uniform_bg)
        bg_track.to_csv(bg_fname, header=False, index=False, sep='\t')
    else:
        bg_track = read_bedgraph(bg_fname)

    rdc_data = dnadataset_to_rdc(dna_dataset, bg_track, rdc_fname, n_cores=n_cores, **rdc_params)

    calculate_scaling_splines(rdc_data, n_cores=n_cores, **scaling_params)

    estimate_significance(rdc_data, n_cores)
    peaks = fetch_peaks(rdc_data, qval_threshold, qval_type, n_cores)
    formatted_peaks = format_peaks(peaks, **peaks_format_params)
    formatted_peaks.to_csv(peaks_output, sep='\t', header=False, index=False)
