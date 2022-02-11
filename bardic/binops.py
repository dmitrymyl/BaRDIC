import bioframe as bf
import numpy as np
import pandas as pd


def make_log_bins(length, log_bin_size):
    """Makes log-scaled bins.
    
    Args:
        length (int, float): region length in log10 scale.
        log_bin_size (int, float): bin size in log10 scale.
    
    Returns:
        np.array: bin edges in linear scale.
    """
    bin_edges = (10 ** np.arange(0, np.log10(length), log_bin_size)).round().astype('int')
    bin_edges[0] = 0
    bin_edges[-1] = length
    return bin_edges


def make_cis_bins(log_bin_size, chrom_name, chrom_length, gene_start, gene_end, min_bin_size=20):
    """Makes cis bins of increasing size from gene borders.
    
    Bins are of equal size in a logarithmic scale, but expanding
    in a linear scale. Bins with linear size less than min_bin_size
    are discarded.

    Args:
        log_bin_size (int, float): bin size in log10 scale.
        chrom_name (str): chromosome name to put in the dataframe.
        chrom_length (int): length of the chromosome.
        gene_start (int): start coordinate of the gene.
        gene_end (int): end coordinate of the gene.
        min_bin_size (int): minimal linear bin size to
            preserve the bin.
    
    Returns:
        pd.DataFrame: cis bins as a dataframe with 3 columns:
            chrom, start, end.
    """
    upstream_size = gene_start
    downstream_size = chrom_length - gene_end

    upstream_edges = make_log_bins(upstream_size, log_bin_size)
    upstream_edges = gene_start - upstream_edges[::-1]
    upstream_df = pd.DataFrame({'chrom': chrom_name, 'start': upstream_edges[:-1], 'end': upstream_edges[1:]})

    downstream_edges = make_log_bins(downstream_size, log_bin_size)
    downstream_edges = gene_end + downstream_edges
    downstream_df = pd.DataFrame({'chrom': chrom_name, 'start': downstream_edges[:-1], 'end': downstream_edges[1:]})

    bins_df = pd.concat((upstream_df, downstream_df), ignore_index=True)
    bins_df = bins_df[bins_df['end'] - bins_df['start'] >= min_bin_size].copy().reset_index(drop=True)
    return bins_df


def make_geom_bins(length, start_size, factor):
    final_index = round(np.log(1 - length * (1 - factor) / start_size) / np.log(factor) - 1)
    bin_edges = np.cumsum(start_size * factor ** np.arange(0, final_index + 1)).astype('int')
    bin_edges = np.insert(bin_edges, 0, 0)
    bin_edges[-1] = length
    return bin_edges


def prune_geom_bins(geom_bins, max_linear_size):
    chrom_length = geom_bins[-1]
    scaled_bins = geom_bins[:np.sum(np.diff(geom_bins) <= max_linear_size) + 1]
    if len(scaled_bins) == len(geom_bins):
        return geom_bins
    else:
        linear_bins = np.array(range(scaled_bins[-1], chrom_length, max_linear_size))
        linear_bins[-1] = chrom_length
        return np.concatenate((scaled_bins, linear_bins[1:]))


def make_cis_bins3(start_size, factor, chrom_name, chrom_length, gene_start, gene_end, max_linear_size=None, fillgene=False):
    """Makes cis bins of increasing size from gene borders in geometric progression.
    
    Bin i is of size start_size * factor ** i.

    Args:
        start_size (int): start bin size.
        factor (int, float): geometric progression factor.
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

    upstream_edges = make_geom_bins(upstream_size, start_size, factor)
    if max_linear_size is not None:
        upstream_edges = prune_geom_bins(upstream_edges, max_linear_size)
    upstream_edges = gene_start - upstream_edges[::-1]
    upstream_df = pd.DataFrame({'chrom': chrom_name, 'start': upstream_edges[:-1], 'end': upstream_edges[1:]})

    downstream_edges = make_geom_bins(downstream_size, start_size, factor)
    if max_linear_size is not None:
        downstream_edges = prune_geom_bins(downstream_edges, max_linear_size)
    downstream_edges = gene_end + downstream_edges
    downstream_df = pd.DataFrame({'chrom': chrom_name, 'start': downstream_edges[:-1], 'end': downstream_edges[1:]})
    
    if fillgene:
        gene_df = pd.DataFrame({'chrom': chrom_name, 'start': gene_start, 'end': gene_end}, index=[0])
        bins_df = pd.concat((upstream_df, gene_df, downstream_df), ignore_index=True)
    else:
        bins_df = pd.concat((upstream_df, downstream_df), ignore_index=True)
    return bins_df


def make_linear_bins(bin_size, chromsizes_df):
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
    for i, (chrom, length) in chromsizes_df.iterrows():
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


def make_trans_bins(bin_size, chromsizes_df, gene_chrom):
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
    return make_linear_bins(bin_size, chromsizes_df.query('chrom != @gene_chrom'))


def calculate_bins_coverage(contacts_df, bins_df):
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


def make_interval_centers(intervals_df):
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


def make_rel_dist(point, start, end):
    if point < start:
        return start - point
    elif point > end:
        return point - end
    else:
        return 0


def make_rel_dist_vector(points, start, end):
    return np.where(points < start, start - points, np.where(points > end, points - end, 0))


