from functools import partial
from typing import Callable, Dict

import bioframe as bf
import pandas as pd
from tqdm.contrib.concurrent import process_map

from ..api.binops import (calculate_bins_coverage, make_cis_bins,
                          make_trans_bins)
from ..api.formats import DnaDataset
from ..api.mp import adjust_chunksize
from ..api.optim import optimize_cost_function


def _optimize_bin_size_single(gene_name: str,
                              dna_dataset: DnaDataset,
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
    dna_df = dna_dataset.read_dna_parts(gene_name)
    chromsizes = dna_dataset.chromsizes
    annotation = dna_dataset.annotation
    gene_annot = annotation.get(gene_name)
    if gene_annot is None:
        print(gene_name, 'no annotation found')
        return dict()
    gene_chrom = gene_annot.chrom
    gene_start = gene_annot.start
    gene_end = gene_annot.end
    chrom_length = chromsizes.get(gene_chrom)
    if chrom_length is None:
        print(gene_name, 'no chromosome found')
        return dict()
    gene_chrom_contacts = len(dna_df.query('chrom == @gene_chrom'))
    genic_contacts = len(bf.select(dna_df, (gene_chrom, gene_start, gene_end)))
    cis_contacts = gene_chrom_contacts - genic_contacts
    trans_contacts = len(dna_df.query('chrom != @gene_chrom'))
    total_contacts = len(dna_df)

    trans_bin_size, trans_status = optimize_cost_function(dna_df,
                                                          make_trans_bins,
                                                          start=trans_min,
                                                          end=trans_max,
                                                          step=trans_step,
                                                          tolerance=tolerance,
                                                          w=w,
                                                          chromsizes=chromsizes,
                                                          gene_chrom=gene_chrom)

    opt_trans_bins = make_trans_bins(trans_bin_size,
                                     chromsizes,
                                     gene_chrom)
    opt_trans_cov = calculate_bins_coverage(opt_trans_bins, dna_df)
    trans_mean = opt_trans_cov['count'].mean()
    trans_var = opt_trans_cov['count'].var(ddof=0)
    trans_zeroes = (opt_trans_cov['count'] == 0).sum()
    n_trans_bins = opt_trans_cov.shape[0]
    trans_zeroes_share = trans_zeroes / n_trans_bins

    cis_factor, cis_status = optimize_cost_function(dna_df,
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
    opt_cis_cov = calculate_bins_coverage(opt_cis_bins, dna_df)
    cis_mean = opt_cis_cov['count'].mean()
    cis_var = opt_cis_cov['count'].var(ddof=0)
    cis_zeroes = (opt_cis_cov['count'] == 0).sum()
    n_cis_bins = opt_cis_cov.shape[0]
    cis_zeroes_share = cis_zeroes / n_cis_bins

    results_dict = {'gene_name': gene_name,
                    'total_contacts': total_contacts,
                    'genic_contacts': genic_contacts,
                    'cis_contacts': cis_contacts,
                    'trans_contacts': trans_contacts,
                    'trans_bin_size': trans_bin_size,
                    'trans_status': trans_status,
                    'n_trans_bins': n_trans_bins,
                    'trans_mean': trans_mean,
                    'trans_var': trans_var,
                    'n_trans_zeroes': trans_zeroes,
                    'trans_zeroes_share': trans_zeroes_share,
                    'cis_factor': cis_factor,
                    'cis_status': cis_status,
                    'n_cis_bins': n_cis_bins,
                    'cis_mean': cis_mean,
                    'cis_var': cis_var,
                    'n_cis_zeroes': cis_zeroes,
                    'cis_zeroes_share': cis_zeroes_share,
                    'cis_start': cis_start,
                    }
    return results_dict


def optimize_bin_sizes(dna_dataset: DnaDataset,
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
                       chunksize: int = 100) -> pd.DataFrame:
    rna_contact_amount = dna_dataset.get_num_contacts()
    rna_eligibility = {rna_name: rna_num_contacts >= n_contacts
                       for rna_name, rna_num_contacts in rna_contact_amount.items()}
    dna_dataset.write_rna_attribute_batch('eligible', rna_eligibility)
    selected_rnas = [rna_name
                     for rna_name, eligible in rna_eligibility.items()
                     if eligible]

    config = {"dna_dataset": dna_dataset,
              "make_trans_bins": make_trans_bins,
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

    func = partial(_optimize_bin_size_single, **config)
    chunksize = adjust_chunksize(len(selected_rnas), n_cores, chunksize)
    bin_selection_results = process_map(func,
                                        selected_rnas,
                                        max_workers=n_cores,
                                        chunksize=chunksize,
                                        desc='Selecting bin sizes',
                                        unit='RNA')

    attr_names = ['genic_contacts',
                  'cis_contacts',
                  'trans_contacts',
                  'trans_bin_size',
                  'cis_factor',
                  'cis_start']
    for attr_name in attr_names:
        attr_vals = {item['gene_name']: item[attr_name] for item in bin_selection_results}
        dna_dataset.write_rna_attribute_batch(attr_name, attr_vals)
    dna_dataset.are_binsizes_selected = True

    selection_df = pd.DataFrame.from_records(bin_selection_results).dropna(how='all')
    return selection_df
