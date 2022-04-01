import os
import sys
from typing import Any, Callable, Dict

import pandas as pd
import ray

from .binops import calculate_bins_coverage, make_cis_bins, make_trans_bins
from .mp import process_by_chunks, add_pbar
from .optim import optimize_cost_function
from .schemas import Coords


os.environ['RAY_OBJECT_STORE_ALLOW_SLOW_STORAGE'] = "1"


def make_selection_results(gene_name: str,
                           contacts_data: pd.DataFrame,
                           chromsizes_dict: Dict[str, int],
                           annot_dict: Dict[str, Coords],
                           make_trans_bins: Callable = make_trans_bins,
                           trans_min: int = 10000,
                           trans_max: int = 1000000,
                           trans_step: int = 1000,
                           make_cis_bins: Callable = make_cis_bins,
                           cis_min: float = 1.1,
                           cis_max: float = 2,
                           cis_step: float = 0.1,
                           cis_start: int = 5000,
                           tolerance: float = 0.01,
                           w: int = 1) -> Dict:
    gene_contacts = contacts_data.query('name == @gene_name').reset_index(drop=True)
    gene_annot = annot_dict.get(gene_name)
    if gene_annot is None:
        print(gene_name, 'no annotation found')
        sys.stdout.flush()
        return dict()
    gene_chrom = gene_annot.chrom
    gene_start = gene_annot.start
    gene_end = gene_annot.end
    chrom_length = chromsizes_dict.get(gene_chrom)
    if chrom_length is None:
        print(gene_name, 'no chromosome found')
        sys.stdout.flush()
        return dict()
    cis_contacts = gene_contacts.query('chrom == @gene_chrom').shape[0]
    trans_contacts = gene_contacts.query('chrom != @gene_chrom').shape[0]
    total_contacts = gene_contacts.shape[0]

    trans_bin_size, trans_status = optimize_cost_function(gene_contacts,
                                                          make_trans_bins,
                                                          start=trans_min,
                                                          end=trans_max,
                                                          step=trans_step,
                                                          tolerance=tolerance,
                                                          w=w,
                                                          chromsizes=chromsizes_dict,
                                                          gene_chrom=gene_chrom)

    opt_trans_bins = make_trans_bins(trans_bin_size,
                                     chromsizes_dict,
                                     gene_chrom)
    opt_trans_cov = calculate_bins_coverage(opt_trans_bins, gene_contacts)
    trans_mean = opt_trans_cov['count'].mean()
    trans_var = opt_trans_cov['count'].var(ddof=0)
    trans_zeroes_share = (opt_trans_cov['count'] == 0).sum() / opt_trans_cov.shape[0]
    
    cis_factor, cis_status = optimize_cost_function(gene_contacts,
                                                    make_cis_bins,
                                                    start=cis_min,
                                                    end=cis_max,
                                                    step=cis_step,
                                                    tolerance=tolerance,
                                                    w=w,
                                                    chrom_name=gene_chrom,
                                                    chrom_length=chrom_length,
                                                    gene_start=gene_start,
                                                    gene_end=gene_end,
                                                    start_size=cis_start,
                                                    max_linear_size=trans_bin_size)

    opt_cis_bins = make_cis_bins(cis_factor,
                                 start_size=cis_start,
                                 chrom_name=gene_chrom,
                                 chrom_length=chrom_length,
                                 gene_start=gene_start,
                                 gene_end=gene_end,
                                 max_linear_size=trans_bin_size)
    opt_cis_cov = calculate_bins_coverage(opt_cis_bins, gene_contacts)
    cis_mean = opt_cis_cov['count'].mean()
    cis_var = opt_cis_cov['count'].var(ddof=0)
    cis_zeroes_share = (opt_cis_cov['count'] == 0).sum() / opt_cis_cov.shape[0]
    

    results_dict = {'gene_name': gene_name,
                    'total_contacts': total_contacts,
                    'cis_contacts': cis_contacts,
                    'trans_contacts': trans_contacts,
                    'trans_bin_size': trans_bin_size,
                    'trans_status': trans_status,
                    'trans_mean': trans_mean,
                    'trans_var': trans_var,
                    'trans_zeroes_share': trans_zeroes_share,
                    'cis_factor': cis_factor,
                    'cis_status': cis_status,
                    'cis_mean': cis_mean,
                    'cis_var': cis_var,
                    'cis_zeroes_share': cis_zeroes_share,
                    'cis_start': cis_start,
                    }
    return results_dict


def select_bin_sizes(contacts_data: pd.DataFrame,
                     chromsizes_dict: Dict[str, int],
                     annot_dict: Dict[str, Coords],
                     n_contacts: int = 1000,
                     make_trans_bins: Callable = make_trans_bins,
                     trans_min: int = 10000,
                     trans_max: int = 1000000,
                     trans_step: int = 1000, 
                     make_cis_bins: Callable = make_cis_bins,
                     cis_min: float = 1.1,
                     cis_max: float = 2.,
                     cis_step: float = 0.01,
                     cis_start: int = 5000,
                     tolerance: float = 0.01,
                     w: int = 1,
                     n_cores: int = 1,
                     n_tasks: int = 1000) -> pd.DataFrame:
    ray.init(num_cpus=n_cores)
    remote_func = ray.remote(add_pbar(make_selection_results))
    rna_contact_amount = contacts_data['name'].value_counts()
    selected_rnas = (rna_contact_amount[rna_contact_amount >= n_contacts]).index.to_list()

    shared_chromsizes_dict = ray.put(chromsizes_dict)
    shared_annot_dict = ray.put(annot_dict)
    shared_contacts_data = ray.put(contacts_data)
    
    config = {"make_trans_bins": make_trans_bins,
              "trans_min": trans_min,
              "trans_max": trans_max,
              "trans_step": trans_step,
              "tolerance": tolerance,
              "w": w,
              "make_cis_bins": make_cis_bins,
              "cis_min": cis_min,
              "cis_max": cis_max,
              "cis_step": cis_step,
              "cis_start": cis_start}

    bin_selection_results = process_by_chunks(selected_rnas,
                                              remote_func,
                                              n_tasks,
                                              shared_contacts_data,
                                              shared_chromsizes_dict,
                                              shared_annot_dict,
                                              **config)

    selection_df = pd.DataFrame.from_records(bin_selection_results).dropna(how='all')
    return selection_df
