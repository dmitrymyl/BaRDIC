import argparse
from pathlib import Path
from urllib.error import HTTPError, URLError

from ..api.convert import annotation_to_dict, chromsizes_to_dict
from ..api.io import fetch_chromsizes, load_annotation, load_chromsizes
from ..utils import (bed2h5, calculate_bin_sizes, calculate_scaling_splines,
                     dnadataset_to_rdc, estimate_significance, fetch_peaks,
                     format_peaks, make_background_track)


def bed2h5_cli(args: argparse.Namespace):
    annotation = annotation_to_dict(load_annotation(args.annotation))
    if Path(args.chromsizes).exists():
        chromsizes = chromsizes_to_dict(load_chromsizes(args.chromsizes))
    else:
        try:
            chromsizes = chromsizes_to_dict(fetch_chromsizes(args.chromsizes))
        except HTTPError:
            raise ValueError(f'{args.chromsizes} file does not exist and is not a valid UCSC genome name.')
        except URLError:
            raise Exception(f"Couldn't fetch chromsizes for {args.chromsizes}, check internet connection.")
    _ = bed2h5(bed_fname=args.dnaparts,
               h5_fname=args.output,
               chromsizes=chromsizes,
               annotation=annotation)


def binsizes_cli(args: argparse.Namespace):
    output_filename = args.output
    del args.output
    selection_results = calculate_bin_sizes(**dict(args._get_kwargs()))
    selection_results.to_csv(output_filename, header=True, index=False, sep='\t')


def background_cli(args: argparse.Namespace):
    output_filename = args.output
    del args.output
    bg_track = make_background_track(...)
    bg_track.to_csv(output_filename, header=True, index=False, sep='\t')


def makerdc_cli(args: argparse.Namespace):
    _ = dnadataset_to_rdc(...)


def scaling_cli(args: argparse.Namespace):
    calculate_scaling_splines(...)


def peaks_cli(args: argparse.Namespace):
    estimate_significance(...)
    peaks = fetch_peaks(...)
    peaks = format_peaks(peaks, ...)
    peaks.to_csv(..., header=False, index=False, sep='\t')
