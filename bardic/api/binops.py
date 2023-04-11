from typing import Dict, Optional, Union, cast

import bioframe as bf
import numpy as np
import pandas as pd

from .schemas import GeneCoord, RnaAttrs


def _make_geom_bins(length: int, start_size: int, factor: Union[int, float]) -> np.ndarray:
    """Makes geometrically increasing bins.

    Args:
        length (int): chromosome lenth.
        start_size (int): start bin size, base of the geometric progression.
        factor (float, int): a common ratio between two successive bin sizes.
            Must be not less than 1.

    Returns:
        np.array: bins edges.

    Raises:
        ValueError: in case factor value is less than 1.
    """
    if factor < 1:
        raise ValueError(f"factor value {factor} is less than 1.")
    elif factor == 1:
        bin_edges = np.arange(0, length, start_size, dtype='int')
    else:
        final_index = round(np.log(1 - length * (1 - factor) / start_size) / np.log(factor) - 1)  # from sum of geometric series
        bin_edges = np.cumsum(start_size * factor ** np.arange(0, final_index + 1)).astype('int')
        bin_edges = np.insert(bin_edges, 0, 0)
    bin_edges[-1] = length
    return bin_edges


def _prune_geom_bins(geom_bins: np.ndarray, max_linear_size: int) -> np.ndarray:
    """Prunes geometric bins to maximal linear size.

    Given output from `make_geom_bins` removes all bins
    that are longer than `max_linear_size` and fills the
    remaining chromosome length with bins of `max_linear_size`
    in length.

    Args:
        geom_bins (np.array): bins edges from `make_geom_bins`.
        max_linear_size (int): maximal linear size.
    """
    chrom_length = geom_bins[-1]
    scaled_idx = cast(int, np.sum(np.diff(geom_bins) <= max_linear_size)) + 1
    scaled_bins = geom_bins[: scaled_idx]
    if len(scaled_bins) == len(geom_bins):
        return geom_bins
    else:
        linear_bins = np.array(range(scaled_bins[-1], chrom_length, max_linear_size))
        linear_bins[-1] = chrom_length
        return np.concatenate((scaled_bins, linear_bins[1:]))


def make_cis_bins(factor: Union[int, float],
                  start_size: int,
                  chrom_name: str,
                  chrom_length: int,
                  gene_start: int,
                  gene_end: int,
                  max_linear_size: Optional[int] = None,
                  fillgene: bool = False) -> pd.DataFrame:
    """Makes cis bins of increasing size from gene borders in geometric progression.

    Bin i is of size start_size * factor ** i.

    Args:
        factor (int, float): geometric progression factor.
        start_size (int): start bin size.
        chrom_name (str): chromosome name to put in the dataframe.
        chrom_length (int): length of the chromosome.
        gene_start (int): start coordinate of the gene.
        gene_end (int): end coordinate of the gene.
        min_bin_size (int): minimal linear bin size to
            preserve the bin.
        max_linear_size (int): maximal linear bin size to stop
            geometric growth of bin size and use given bin size
            hereafter (default: None).
        fillgene (bool): if True, will add a single bin of gene
            coordinates, otherwise will exclude gene from binning
            (default: False).

    Returns:
        pd.DataFrame: cis bins as a dataframe with 3 columns:
            chrom, start, end.

    Raises:
        ValueError: in case factor < 1.
    """
    if factor < 1:
        raise ValueError(f"factor value {factor} is less than 1.")
    upstream_size = gene_start
    downstream_size = chrom_length - gene_end

    upstream_edges = _make_geom_bins(upstream_size, start_size, factor)
    if max_linear_size is not None:
        upstream_edges = _prune_geom_bins(upstream_edges, max_linear_size)
    upstream_edges = gene_start - upstream_edges[::-1]
    upstream_df = pd.DataFrame({'chrom': chrom_name, 'start': upstream_edges[:-1], 'end': upstream_edges[1:]})

    downstream_edges = _make_geom_bins(downstream_size, start_size, factor)
    if max_linear_size is not None:
        downstream_edges = _prune_geom_bins(downstream_edges, max_linear_size)
    downstream_edges = gene_end + downstream_edges
    downstream_df = pd.DataFrame({'chrom': chrom_name, 'start': downstream_edges[:-1], 'end': downstream_edges[1:]})

    if fillgene:
        gene_df = pd.DataFrame({'chrom': chrom_name, 'start': gene_start, 'end': gene_end}, index=[0])
        bins_df = pd.concat((upstream_df, gene_df, downstream_df), ignore_index=True)
    else:
        bins_df = pd.concat((upstream_df, downstream_df), ignore_index=True)
    return bins_df


def make_linear_bins(bin_size: int, chromsizes: Dict[str, int]) -> pd.DataFrame:
    """Makes equal sized bins in linear scale for given chromosomes.

    Args:
        bin_size (int): linear bin size.
        chromsizes_df (pd.DataFrame): chromosome sizes dataframe
            with 2 columns: "chrom" (chromosome name) and "length"
            (chromosome length).

    Returns:
        pd.DataFrame: a dataframe with bins of 3 columns:
            "chrom", "start", "end".
    """
    chrom_dfs = list()
    bin_size = round(bin_size)
    for chrom, length in chromsizes.items():
        break_points = list(range(0, length, bin_size))
        starts = break_points[:-1]
        ends = break_points[1:]
        if len(break_points) <= 1:
            starts = [0]
            ends = [length]
        else:
            if ends[-1] < length:
                starts.append(ends[-1])
                ends.append(length)
        chrom_df = pd.DataFrame(((chrom, start, end)
                                 for start, end in zip(starts, ends)),
                                columns=['chrom', 'start', 'end'])
        chrom_dfs.append(chrom_df)
    total_bins = pd.concat(chrom_dfs, ignore_index=True)
    total_bins['chrom'] = total_bins['chrom'].astype('category')
    return total_bins


def make_trans_bins(bin_size: int, chromsizes: Dict[str, int], gene_chrom: str) -> pd.DataFrame:
    """Makes trans bins of equal size.

    Uses `make_linear_bins` for binning.

    Args:
        bin_size (int): linear bin size.
        chromsizes_df (pd.DataFrame): chromosome sizes dataframe
            with 2 columns: "chrom" (chromosome name) and "length"
            (chromosome length).
        gene_chrom (str): a gene chromosome to exclude it
            from binning.

    Returns:
        pd.DataFrame: a dataframe with bins of 3 columns:
            "chrom", "start", "end".
    """
    _chromsizes = chromsizes.copy()
    del _chromsizes[gene_chrom]
    return make_linear_bins(bin_size, _chromsizes)


def calculate_bins_coverage(bins_df: pd.DataFrame, contacts_df: pd.DataFrame) -> pd.DataFrame:
    """Calculates how many contacts overlap every bin.

    Args:
        contacts_df (pd.DataFrame): a contacts dataframe
            of 3 columns: "chrom", "start", "end".
        bins_df (pd.DataFrame): a bins dataframe of
            3 columns: "chrom", "start", "end".

    Returns:
        pd.DataFrame: a 4-column dataframe
            ("chrom", "start", "end", "count") with
            bins from bins_df and amount of contacts
            overlapping with each of them.
    """
    return bf.count_overlaps(bins_df, contacts_df)


def make_interval_centers(intervals_df: pd.DataFrame) -> pd.DataFrame:
    """Make annotation of centers of provided intervals.

    Args:
        intervals_df (pd.DataFrame): a 3-column dataframe with
            genomic intervals: "chrom", "start", "end".

    Returns:
        pd.DataFrame: a 3-column dataframe ("chrom", "start", "end")
            of interval centers.
    """
    intervals_centers = (intervals_df['end'] + intervals_df['start']) // 2
    intervals_centers_df = pd.DataFrame({'chrom': intervals_df['chrom'],
                                         'start': intervals_centers,
                                         'end': intervals_centers + 1})
    return intervals_centers_df


def calculate_dists_to_points(points: np.ndarray,
                              start: int,
                              end: int) -> np.ndarray:
    """Finds relative distance of points to gene with [start; end] coordinates.

    Args:
        points (np.array): genomic points.
        start (int): gene start.
        end (int): gene end.

    Returns:
        np.array: relative distances.
    """
    return np.where(points < start, start - points, np.where(points > end, points - end, 0))


def calculate_dists_to_centers(cis_bins_df: pd.DataFrame,
                               gene_start: int,
                               gene_end: int) -> np.ndarray:
    """Finds relative distance of bins centers to gene with [start; end] coordinates.

    Args:
        cis_bins_df (pd.DataFrame): genomic bins on the same chromosome as gene.
        start (int): gene start.
        end (int): gene end.

    Returns:
        np.array: relative distances.
    """
    bins_centers = (cis_bins_df['start'] + cis_bins_df['end']) // 2
    rel_dists = calculate_dists_to_points(bins_centers, gene_start, gene_end)
    return rel_dists


def make_track(bins_df: pd.DataFrame,
               centers_df: pd.DataFrame,
               bg_track: pd.DataFrame,
               ivalue: Optional[Union[int, float]] = None,
               bg_count_name: str = 'count') -> pd.DataFrame:
    """Computes signal and bg track in given bins.

    Args:
        bins_df (pd.DataFrame): a bed3 df with bins.
        centers_df (pd.DataFrame): a contacts centers df.
        bg_track (pd.DataFrame): a bed4 df with binned bg counts
            ('chrom', 'start', 'end', 'count').
        impute (bool): whether to impute missing bg values or not
            (default: False).
        ivalue (float): an imputation value as a bg count per nt
            (default: None).
        bg_count_name (str): 4th column name in bg_track (default: 'count').

    Returns:
        pd.DataFrame: a df with 7 columns: 'chrom', 'start', 'end',
            'count', 'signal_prob', 'bg_count', 'bg_prob'.

    Raises:
        ValueError: in case `impute` is True and `ivalue` is None.
    """
    bins_coverage = calculate_bins_coverage(bins_df, centers_df)
    bins_coverage.rename(columns={'count': 'signal_count'}, inplace=True)
    bins_coverage['signal_prob'] = bins_coverage['signal_count'] / bins_coverage['signal_count'].sum()
    overlap_with_bg = bf.overlap(bins_coverage,
                                 bg_track,
                                 how="left",
                                 return_input=True,
                                 keep_order=True,
                                 suffixes=('', '_bg'),
                                 return_index=True)\
                        .groupby('index')\
                        .agg({'start': 'min',
                              'end': 'max',
                              'start_bg': 'min',
                              'end_bg': 'max',
                              f'{bg_count_name}_bg': 'sum'})
    bin_sizes = bins_coverage['end'] - bins_coverage['start']
    bg_bin_sizes = overlap_with_bg['end_bg'].astype('int64') - overlap_with_bg['start_bg'].astype('int64')
    rescaled_bg_counts = overlap_with_bg[f'{bg_count_name}_bg'].astype('int64') / bg_bin_sizes * bin_sizes
    if ivalue is not None:
        if not (isinstance(ivalue, int) or isinstance(ivalue, float)):
            raise ValueError("ivalue must be int or float.")
        dropout = (bins_coverage['signal_count'] > 0) & (rescaled_bg_counts == 0)
        imputed_rescaled_bg_counts = np.where(dropout,
                                              bin_sizes * ivalue,
                                              rescaled_bg_counts)
        bins_coverage['impute'] = dropout
        bins_coverage['bg_count'] = imputed_rescaled_bg_counts
    else:
        bins_coverage['bg_count'] = rescaled_bg_counts
        bins_coverage['impute'] = False
    bins_coverage['bg_prob'] = bins_coverage['bg_count'] / bins_coverage['bg_count'].sum()
    return bins_coverage


def make_genomic_track(dna_contacts: pd.DataFrame,
                       bg_track: pd.DataFrame,
                       chromdict: Dict[str, int],
                       gene_coord: GeneCoord,
                       rna_attrs: RnaAttrs,
                       ivalue: Optional[Union[int, float]] = None) -> pd.DataFrame:
    gene_chrom, gene_start, gene_end = gene_coord.chrom, gene_coord.start, gene_coord.end
    trans_bin_size = rna_attrs.trans_bin_size
    cis_factor = rna_attrs.cis_factor
    cis_start = rna_attrs.cis_start
    dna_parts_centers = make_interval_centers(dna_contacts)

    # trans bins
    trans_bins = make_trans_bins(trans_bin_size, chromdict, gene_chrom)
    trans_coverage = make_track(trans_bins,
                                dna_parts_centers,
                                bg_track,
                                ivalue=ivalue)
    trans_coverage['raw_bg_prob'] = trans_coverage['bg_prob']
    trans_coverage['scaling_factor'] = 1

    # Cis bins
    cis_bins = make_cis_bins(cis_factor,
                             cis_start,
                             gene_chrom,
                             chromdict[gene_chrom],
                             gene_start,
                             gene_end,
                             max_linear_size=trans_bin_size)
    cis_coverage = make_track(cis_bins,
                              dna_parts_centers,
                              bg_track,
                              ivalue=ivalue)
    cis_coverage['raw_bg_prob'] = cis_coverage['bg_prob']
    cis_coverage['scaling_factor'] = 1

    # combine
    total_track = bf.sort_bedframe(pd.concat((trans_coverage, cis_coverage), ignore_index=True))
    return total_track
