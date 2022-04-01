from typing import Dict, Tuple
from .binops import make_cis_bins, make_linear_bins, make_trans_bins, make_track
import pandas as pd
import bioframe as bf


def encode_chromsizes(chromsizes_df: pd.DataFrame) -> Dict[str, int]:
    return chromsizes_df.set_index('chrom').to_dict()['length']


def decode_chromsizes(chrom_dict: Dict[str, int]) -> pd.DataFrame:
    return pd.DataFrame.from_dict(chrom_dict, orient='index').reset_index().rename(columns={'index': 'chrom', 0: 'length'})


def encode_bg_track(bg_track: pd.DataFrame) -> Dict[str, Tuple[int, ...]]:
    bg_dict = {chrom: tuple(vector.values)
               for chrom, vector in bg_track.groupby('chrom')['score']}
    return bg_dict


def decode_bg_track(bg_dict: Dict[str, Tuple[int, ...]],
                    binsize: int,
                    chromsizes: Dict[str, int]) -> pd.DataFrame:
    bg_bins = make_linear_bins(binsize, decode_chromsizes(chromsizes))
    bg_subtracks = list()
    for chrom, bg_subbins in bg_bins.groupby('chrom'):
        bg_subbins['score'] = bg_dict[chrom]
        bg_subtracks.append(bg_subbins)
    bg_track = pd.concat(bg_subtracks, ignore_index=True)
    return bg_track


def encode_track(track: pd.DataFrame) -> Dict:
    data_cols = set(track.columns) - set(('chrom', 'start', 'end'))
    return {chrom: {key: tuple(value)
                    for key, value in data.loc[:, data_cols].to_dict('list').items()}
            for chrom, data in track.groupby('chrom')}


def encode_rna_track(dna_parts: pd.DataFrame,
                     rna_name: str,
                     chrom_dict: Dict[str, int],
                     selection_dict: Dict,
                     annot_dict: Dict,
                     bg_track: pd.DataFrame,
                     impute: bool = False,
                     ifactor: float = 0):
    params = {"rna_cis_factor": selection_dict[rna_name]['cis_factor'],
              "rna_trans_size": selection_dict[rna_name]['trans_bin_size'],
              "cis_start_size": selection_dict[rna_name]['cis_start'],
              "rna_chrom" : annot_dict[rna_name]['chrom'],
              "rna_chrom_length": chrom_dict[annot_dict[rna_name]['chrom']],
              "rna_start": annot_dict[rna_name]['start'],
              "rna_end": annot_dict[rna_name]['end']}
    track = make_track(dna_parts,
                       rna_name,
                       chrom_dict,
                       selection_dict,
                       annot_dict,
                       bg_track,
                       impute,
                       ifactor)
    track_dict = encode_track(track)
    track_dict['params'] = params
    return track_dict


def decode_rna_track(track_dict: Dict,
                     trans_bin_size: int,
                     cis_factor: float,
                     start_size: int,
                     chrom_dict: Dict[str, int],
                     chrom_name: str,
                     start: int, 
                     end: int,
                     fillgene: bool = False) -> pd.DataFrame:
    trans_bins = make_trans_bins(trans_bin_size, decode_chromsizes(chrom_dict), chrom_name)
    cis_bins = make_cis_bins(cis_factor, start_size, chrom_name, chrom_dict[chrom_name], start, end, trans_bin_size, fillgene=fillgene)
    total_bins = pd.concat((cis_bins, trans_bins), ignore_index=True)
    total_bins = bf.sort_bedframe(total_bins)
    sub_bins_list = list()
    for chrom, chrom_bins in total_bins.groupby('chrom'):
        sub_bins = pd.concat((chrom_bins.reset_index(drop=True), pd.DataFrame.from_dict(track_dict[chrom])), axis=1)
        sub_bins_list.append(sub_bins)
    return pd.concat(sub_bins_list, ignore_index=True)
