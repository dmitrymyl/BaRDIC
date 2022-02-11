import numpy as np
import warnings

from .binops import calculate_bins_coverage, make_interval_centers


def calculate_cost_function_from_bins(bins_coverage, bs):
    """Caculates cost function from bin coverage dataframe.
    
    Args:
        bins_coverage (pd.DataFrame): an output from
            `calculate_bins_coverage` function.
        bs (int, float): a bin size used to create
            bins dataframe.
    
    Returns:
        float: a cost function value.
    """
    mean_cov = bins_coverage['count'].mean()
    var_cov = bins_coverage['count'].var(ddof=0)
    cost_func_value = (2 * mean_cov - var_cov) / (bs) ** 2
    return cost_func_value


def calculate_cost_function(contacts_df, bin_size, binner, **kwargs):
    """Calculates cost function for given contacts and bin size.
    
    Args:
        contacts_df (pd.DataFrame): a contacts dataframe
            of 3 columns: "chrom", "start", "end".
        bin_size (int): bin size.
        binnner (callable): a binning function (one of `make_cis_bins`, `make_trans_bins`).
        kwargs: keyword arguments used for specific binner.
    
    Returns:
        float: a cost function value.
    """
    bins_df = binner(bin_size, **kwargs)
    bins_coverage = calculate_bins_coverage(contacts_df, bins_df)
    return calculate_cost_function_from_bins(bins_coverage, bin_size)


def optimize_cost_function(contacts_df, binner, start=1000, end=10000000, step=1000, tolerance=0.01, **kwargs):
    """Optimizes cost function for given contacts.

    Tests multiple bin sizes from the smallest to the biggest
    with given step and returns the one which
    minimizes cost function. There can be 3 outcomes:
    1. cost function declines too slow (abs value of [c(i+1) - c(i)] / c(i) is
        less than tolerance), which means optimization has converged;
    2. cost function values start rising, which means the local minimum was
        found;
    3. The first two outcomes are not met up to the biggest bin size,
        which means optimization has not converged.

    Cis and trans bins are optimized separately with specific binners.
    
    Args:
        contacts_df (pd.DataFrame): a contacts dataframe
            of 3 columns: "chrom", "start", "end".
        binnner (callable): a binning function (one of `make_cis_bins`, `make_trans_bins`).
        start (int, float): the smallest bin size to test (default: 1000 for trans bins).
        end (int, float): the biggest bin size to test (default: 10000000 for trans bins).
        tolerance (float): ...
        kwargs: keyword arguments used for specific binner.
    
    Returns:
        int, float: the best bin size.
        str: optimization status, one of "converged", "rising", "not converged".
    """
    warnings.filterwarnings('ignore')
    prev_bin_size = 0
    prev_cost_value = float('inf')
    contacts_centers_df = make_interval_centers(contacts_df)
    bin_sizes = np.around(np.arange(start, end + step, step), -int(np.log10(min(start, step))))
    for cur_bin_size in bin_sizes:
        cur_cost_value = calculate_cost_function(contacts_centers_df,
                                                 cur_bin_size,
                                                 binner,
                                                 **kwargs)
        diff = cur_cost_value - prev_cost_value
        abs_relative_diff = abs(diff / prev_cost_value)

        if prev_cost_value == float('inf') or diff == float('-inf'):
            prev_cost_value = cur_cost_value
            prev_bin_size = cur_bin_size
            continue
        elif diff > 0:
            return prev_bin_size, 'rising'
        elif abs_relative_diff <= tolerance:
            return prev_bin_size, 'converged'
        else:
            prev_bin_size = cur_bin_size
            cur_bin_size = prev_bin_size + step
            prev_cost_value = cur_cost_value
    return end, 'not converged'
