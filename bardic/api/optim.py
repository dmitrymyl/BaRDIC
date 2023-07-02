from collections import deque
from typing import Callable, Deque, Tuple, Union

import numpy as np
import pandas as pd

from .binops import calculate_bins_coverage, make_interval_centers


def _calculate_cost_function_from_bins(bins_coverage: pd.DataFrame,
                                       bs: Union[int, float]) -> float:
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


def _calculate_cost_function(contacts_df: pd.DataFrame,
                             bin_size: Union[int, float],
                             binner: Callable,
                             **kwargs) -> float:
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
    bins_coverage = calculate_bins_coverage(bins_df, contacts_df)
    return _calculate_cost_function_from_bins(bins_coverage, bin_size)


def optimize_cost_function(contacts_df: pd.DataFrame,
                           binner: Callable,
                           start: Union[int, float] = 1000,
                           end: Union[int, float] = 100_000,
                           step: Union[int, float] = 1000,
                           tolerance: float = 0.01,
                           w: int = 5,
                           **kwargs) -> Tuple[Union[int, float], str]:
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
        end (int, float): the biggest bin size to test (default: 100000 for trans bins).
        step (int, float): bin size increasing step (default: 1000 for trans bins).
        tolerance (float): a tolerance parameter (default: 0.01).
        w (int): cost function smoothing window size (default: 5).
        kwargs: keyword arguments used for specific binner.

    Returns:
        int, float: the best bin size.
        str: optimization status, one of "converged", "rising", "not converged".
    """
    contacts_centers_df = make_interval_centers(contacts_df)
    bin_sizes = np.around(np.arange(start, end + step, step), -int(np.log10(min(start, step))))
    deq: Deque[float] = deque()
    for i, cur_bin_size in enumerate(bin_sizes):
        cur_cost_value = _calculate_cost_function(contacts_centers_df,
                                                  cur_bin_size,
                                                  binner,
                                                  **kwargs)
        if i < w:
            deq.append(cur_cost_value)
            prev_bin_size = cur_bin_size
            continue
        else:
            if i == w:
                prev_mean_cost = sum(deq) / len(deq)
            left_cost_value = deq.popleft()
            deq.append(cur_cost_value)
            cur_mean_cost = prev_mean_cost + 1 / w * (cur_cost_value - left_cost_value)
            diff = cur_mean_cost - prev_mean_cost
            if prev_mean_cost == 0:
                prev_bin_size = cur_bin_size
                prev_mean_cost = cur_mean_cost
                continue
            abs_relative_diff = abs(diff / prev_mean_cost)
        if diff > 0:
            return prev_bin_size, 'rising'
        elif abs_relative_diff <= tolerance:
            return prev_bin_size, 'converged'
        else:
            prev_bin_size = cur_bin_size
            prev_mean_cost = cur_mean_cost
    return end, 'not converged'
