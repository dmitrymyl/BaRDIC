from typing import Union
import os

from ..api.convert import annotation_to_dict
from ..api.formats import DnaDataset, Rdc
from ..api.io import get_chromsizes, read_annotation, read_bedgraph
from ..utils import (bed2h5, optimize_bin_sizes, calculate_scaling_splines,
                     dnadataset_to_rdc, estimate_significance, fetch_peaks,
                     format_peaks, make_background_track, run_pipeline, simulate)


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
                no_scaling: bool = False,
                n_cores: int = 1):
    rdc_data = Rdc(rdc)
    calculate_scaling_splines(rdc_data=rdc_data,
                              degree=degree,
                              no_refine=no_refine,
                              max_threshold=max_threshold,
                              fill_value=fill_value,
                              no_scaling=no_scaling,
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
                     no_scaling: bool = False,
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
        make_bg = True
        uniform_bg = False
    elif bg_type == "custom":
        bg_fname = bgdata
        rna_list = None
        make_bg = False
        uniform_bg = True
    elif bg_type == "uniform":
        bg_fname = os.path.join(outdir, 'background.bedGraph')
        rna_list = None
        make_bg = True
        uniform_bg = True
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
                                     fill_value=fill_value,
                                     no_scaling=no_scaling),
                 qval_threshold=qval_threshold,
                 qval_type=qval_type,
                 peaks_format_params=dict(format=format, score=score),
                 peaks_output=peaks_output,
                 n_cores=n_cores,
                 make_bg=make_bg,
                 uniform_bg=uniform_bg)


Numeric = Union[int, float]


def simulate_cli(outdir: str,
                 rna_name: str = 'simulRNA',
                 L_cis: int = 100_000_000,
                 L_genome: int = 2_000_000_000,
                 bin_size: int = 1_000,
                 bg_exp: Numeric = 8,
                 imputation_factor: float = 0.01,
                 A: Numeric = 3,
                 D_high: Numeric = 3,
                 D_low: Numeric = 7,
                 F_high: float = 0.9,
                 F_low: float = 0.1,
                 B_left: float = 0.05,
                 B_right: float = 0.95,
                 N_total: int = 20_000,
                 frac_cis: float = 0.75,
                 N_peaks: int = 5,
                 frac_cis_peaks: float = 0.2,
                 peak_exp: int = 100,
                 sd: int = 5_000,
                 ):

    simulate(outdir=outdir,
             rna_name=rna_name,
             L_cis=L_cis,
             L_genome=L_genome,
             bin_size=bin_size,
             bg_exp=bg_exp,
             imputation_factor=imputation_factor,
             A=A,
             D_high=D_high,
             D_low=D_low,
             F_high=F_high,
             F_low=F_low,
             B_left=B_left,
             B_right=B_right,
             N_total=N_total,
             frac_cis=frac_cis,
             N_peaks=N_peaks,
             frac_cis_peaks=frac_cis_peaks,
             peak_exp=peak_exp,
             sd=sd,
             )
