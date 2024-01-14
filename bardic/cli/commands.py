from typing import Union
import os

from ..api.convert import annotation_to_dict
from ..api.formats import DnaDataset, Rdc
from ..api.io import get_chromsizes, read_annotation, read_bedgraph
from ..utils import (bed2h5, optimize_bin_sizes, calculate_scaling_splines,
                     dnadataset_to_rdc, estimate_significance, fetch_peaks,
                     format_peaks, make_background_track, run_pipeline)


def bed2h5_cli(annotation: str, chromsizes: str, dnaparts: str, output: str):
    annotation_dict = annotation_to_dict(read_annotation(annotation))
    chromsizes_dict = get_chromsizes(chromsizes)
    _ = bed2h5(bed_fname=dnaparts,
               h5_fname=output,
               chromsizes=chromsizes_dict,
               annotation=annotation_dict)


def binsizes_cli(dna_dataset: str,
                 output: str,
                 n_contacts: int = 1000,
                 trans_min: int = 10000,
                 trans_max: int = 1000000,
                 trans_step: int = 1000,
                 cis_min: float = 1.1,
                 cis_max: float = 2.,
                 cis_step: float = 0.01,
                 cis_start: int = 5000,
                 tolerance: float = 0.01,
                 w: int = 1,
                 n_cores: int = 1):
    dna_dataset_obj = DnaDataset(dna_dataset)
    selection_results = optimize_bin_sizes(dna_dataset=dna_dataset_obj,
                                           n_contacts=n_contacts,
                                           trans_min=trans_min,
                                           trans_max=trans_max,
                                           trans_step=trans_step,
                                           cis_min=cis_min,
                                           cis_max=cis_max,
                                           cis_step=cis_step,
                                           cis_start=cis_start,
                                           tolerance=tolerance,
                                           w=w,
                                           n_cores=n_cores)
    selection_results.to_csv(output, header=True, index=False, sep='\t')


def background_cli(dna_dataset: str,
                   rnas: str,
                   binsize: int,
                   output: str):
    dna_dataset_obj = DnaDataset(dna_dataset)
    with open(rnas, 'r') as infile:
        rnas_list = [line.strip() for line in infile]
    bg_track = make_background_track(dna_dataset=dna_dataset_obj,
                                     rnas=rnas_list,
                                     binsize=binsize)
    bg_track.to_csv(output, header=False, index=False, sep='\t')


def makerdc_cli(dna_dataset: str,
                bg: str,
                ifactor: float,
                output: str,
                n_cores: int = 1):
    dna_dataset_obj = DnaDataset(dna_dataset)
    bg_track = read_bedgraph(bg)
    _ = dnadataset_to_rdc(dna_dataset=dna_dataset_obj,
                          bg_track=bg_track,
                          fname=output,
                          ifactor=ifactor,
                          n_cores=n_cores)


def scaling_cli(rdc: str,
                degree: int,
                max_threshold: float,
                no_refine: bool = False,
                fill_value: Union[int, float] = 1,
                n_cores: int = 1):
    rdc_data = Rdc(rdc)
    calculate_scaling_splines(rdc_data=rdc_data,
                              degree=degree,
                              no_refine=no_refine,
                              max_threshold=max_threshold,
                              fill_value=fill_value,
                              n_cores=n_cores)


def peaks_cli(rdc: str,
              qval_threshold: float,
              qval_type: 'str',
              output: str,
              format: str,
              score: Union[str, int],
              n_cores: int):
    rdc_data = Rdc(rdc)
    if not rdc_data.are_peaks_estimated:
        estimate_significance(rdc_data, n_cores)
    peaks = fetch_peaks(rdc_data, qval_threshold, qval_type, n_cores)
    peaks = format_peaks(peaks, format, score)
    peaks.to_csv(output, header=False, index=False, sep='\t')


def run_pipeline_cli(dnaparts: str,
                     annotation: str,
                     chromsizes: str,
                     bgdata: str,
                     outdir: str,
                     n_contacts: int = 1000,
                     trans_min: int = 10000,
                     trans_max: int = 1000000,
                     trans_step: int = 1000,
                     cis_min: float = 1.1,
                     cis_max: float = 2.,
                     cis_step: float = 0.01,
                     cis_start: int = 5000,
                     tolerance: float = 0.01,
                     w: int = 1,
                     n_cores: int = 1,
                     binsize: int = 5000,
                     ifactor: float = 0.01,
                     degree: int = 3,
                     max_threshold: float = 0.05,
                     no_refine: bool = False,
                     fill_value: Union[int, float] = 1,
                     qval_threshold: float = 0.05,
                     qval_type: str = 'global',
                     format: str = "narrowPeak",
                     score: Union[str, int] = 0,
                     bg_type: str = "rnas"):

    if not os.path.exists(outdir):
        os.mkdir(outdir)
    dna_dataset_fname = os.path.join(outdir, "DnaDataset.dnah5")
    rdc_fname = os.path.join(outdir, "contacts.rdc")
    selection_results_fname = os.path.join(outdir, "selection.tsv")

    if bg_type == "rnas":
        bg_fname = os.path.join(outdir, 'background.bedGraph')
        with open(bgdata, 'r') as infile:
            rna_list = [line.strip() for line in infile]
        makebg = True
    elif bg_type == "custom":
        bg_fname = bgdata
        rna_list = None
        makebg = False
    else:
        raise ValueError

    peaks_output = os.path.join(outdir, "peaks." + format)

    annotation_df = read_annotation(annotation)
    chromsizes_dict = get_chromsizes(chromsizes)

    run_pipeline(dna_parts_fname=dnaparts,
                 dna_dataset_fname=dna_dataset_fname,
                 rdc_fname=rdc_fname,
                 chromsizes=chromsizes_dict,
                 annotation=annotation_df,
                 binsize_params=dict(n_contacts=n_contacts,
                                     trans_min=trans_min,
                                     trans_max=trans_max,
                                     trans_step=trans_step,
                                     cis_min=cis_min,
                                     cis_max=cis_max,
                                     cis_step=cis_step,
                                     cis_start=cis_start,
                                     tolerance=tolerance,
                                     w=w),
                 selection_results_fname=selection_results_fname,
                 bg_fname=bg_fname,
                 rna_list=rna_list,
                 bg_binsize=binsize,
                 rdc_params=dict(ifactor=ifactor),
                 scaling_params=dict(degree=degree,
                                     max_threshold=max_threshold,
                                     no_refine=no_refine,
                                     fill_value=fill_value),
                 qval_threshold=qval_threshold,
                 qval_type=qval_type,
                 peaks_format_params=dict(format=format, score=score),
                 peaks_output=peaks_output,
                 n_cores=n_cores,
                 makebg=makebg)
