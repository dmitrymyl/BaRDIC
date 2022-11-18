import argparse
import json
import math
import os
from collections import namedtuple

import numpy as np
import pandas as pd
import ray
import scipy.interpolate as si
import scipy.stats as ss
import statsmodels.api as sm

from bardic.api.binops import (calculate_rel_dist_from_centers, compute_track,
                               make_cis_bins3, make_interval_centers,
                               make_rel_dist_vector, make_trans_bins)
from bardic.api.io import load_bedgraph, load_dna_parts_of_rdc2
from bardic.legacy.mp import process_by_chunks
from bardic.api.schemas import rdc_gene_dtypes, rdc_gene_schema

os.environ['RAY_OBJECT_STORE_ALLOW_SLOW_STORAGE'] = "1"

parser = argparse.ArgumentParser(description='Calculate peaks')
input_group = parser.add_argument_group('Input')
input_group.add_argument('contacts',
                         nargs='?',
                         type=str,
                         help='RDC contacts file')
input_group.add_argument('annotation',
                         nargs='?',
                         type=str,
                         help='Annotation in rdc-compatible format')
input_group.add_argument('chromsizes',
                         nargs='?',
                         type=str,
                         help='chromsizes file')
input_group.add_argument('selection',
                         nargs='?',
                         type=str,
                         help='selection file')
input_group.add_argument('scaling',
                         nargs='?',
                         type=str,
                         help='scaling json file')
input_group.add_argument('background',
                         nargs='?',
                         type=str,
                         help='Binarized background coverage in bed3 format')
output_group = parser.add_argument_group('Output')
output_group.add_argument('output',
                          nargs='?',
                          type=str,
                          help='output tsv filename')
processing_group = parser.add_argument_group('Processing')
processing_group.add_argument('-ncores',
                              nargs='?',
                              type=int,
                              default=1,
                              help='Number of cores for parallelization')
processing_group.add_argument('-ntasks',
                              nargs='?',
                              type=int,
                              default=1000,
                              help='Number of tasks in a single batch to process simultaniously')
processing_group.add_argument('-cisStartSize',
                              nargs='?',
                              type=int,
                              default=1000,
                              help='Cis bins start size')
processing_group.add_argument('-fillValue',
                              nargs='?',
                              type=float,
                              default=1,
                              help='0/0 fold-change fill value')
processing_group.add_argument('-imputationFactor',
                              nargs='?',
                              type=float,
                              default=0.01,
                              help='imputation value multiplication factor')

args = parser.parse_args()
contacts_path = args.contacts
annotation_path = args.annotation
chromsizes_path = args.chromsizes
selection_path = args.selection
scaling_path = args.scaling
background_path = args.background
output_path = args.output
n_cores = args.ncores
n_tasks = args.ntasks
cis_start_size = args.cisStartSize
fill_value = args.fillValue
imp_factor = args.imputationFactor


Coords = namedtuple('Coords', ['chrom', 'start', 'end'])
contacts_df = load_dna_parts_of_rdc2(contacts_path)
chromsizes_df = pd.read_csv(chromsizes_path,
                            header=None,
                            index_col=None,
                            sep='\t',
                            names=['chrom', 'length'])
chromsizes_df = chromsizes_df.query('chrom != "chrM"').reset_index(drop=True)
chromsizes_dict = {item[1]: item[2]
                   for item in chromsizes_df.to_records()}
annotation = pd.read_csv(annotation_path,
                         header=None,
                         index_col=None,
                         sep='\t',
                         names=rdc_gene_schema,
                         dtype=rdc_gene_dtypes)
gene_coords = {item[1]: Coords(item[2], item[3], item[4])
               for item in annotation.loc[:, ['gene_name', 'gene_chr', 'gene_start', 'gene_end']].to_records()}
selection_results = pd.read_csv(selection_path,
                                header=0,
                                index_col=None,
                                sep='\t')
bin_sizes = {item[1]: (item[2], item[3])
             for item in selection_results.loc[:, ['gene_name', 'trans_bin_size', 'cis_factor']].to_records()}
with open(scaling_path, 'r') as infile:
    scaling_dict = json.load(infile)
background_track = load_bedgraph(background_path)
bg_bin_sizes = background_track['end'] - background_track['start']
mean_bg_count_per_nt = (background_track['score'] / bg_bin_sizes).mean()
imputation_bg = mean_bg_count_per_nt * imp_factor

selected_rnas = [rna_name
                 for rna_name in bin_sizes
                 if rna_name in contacts_df['name'].values and gene_coords[rna_name].chrom in scaling_dict]

ray.init(num_cpus=n_cores)

shared_contacts_df = ray.put(contacts_df)
shared_bin_sizes = ray.put(bin_sizes)
shared_annotation_dict = ray.put(gene_coords)
shared_chromsizes_dict = ray.put(chromsizes_dict)
shared_chromsizes_df = ray.put(chromsizes_df)
shared_scaling_dict = ray.put(scaling_dict)
shared_background_df = ray.put(background_track)


@ray.remote
def calculate_peaks(rna_name,
                    bin_sizes_r,
                    contacts_df_r,
                    annotation_dict_r,
                    chromsizes_df_r,
                    chromsizes_dict_r,
                    background_df_r,
                    scaling_dict_r,
                    imputation_bg_r,
                    cis_start_size_r,
                    pba):
    try:
        # preparation
        gene_chrom, gene_start, gene_end = annotation_dict_r.get(rna_name)
        trans_bin_size, cis_bin_factor = bin_sizes_r.get(rna_name)
        rna_contacts = contacts_df_r.query("name == @rna_name").reset_index(drop=True)
        rna_contacts_centers = make_interval_centers(rna_contacts)
        # # trans bins
        trans_bins = make_trans_bins(trans_bin_size, chromsizes_df_r, gene_chrom)
        trans_coverage = compute_track(trans_bins,
                                       rna_contacts_centers,
                                       background_df_r,
                                       impute=True,
                                       imputation_bg=imputation_bg_r,
                                       bg_count_name='score')
        trans_contacts_num = trans_coverage['count'].sum()
        trans_coverage['raw_bg_prob'] = trans_coverage['bg_prob']
        trans_coverage['scaling_factor'] = 1

        # Cis bins
        cis_bins = make_cis_bins3(cis_bin_factor,
                                  cis_start_size_r,
                                  gene_chrom,
                                  chromsizes_dict_r[gene_chrom],
                                  gene_start,
                                  gene_end,
                                  max_linear_size=trans_bin_size)
        cis_coverage = compute_track(cis_bins,
                                     rna_contacts_centers,
                                     background_df_r,
                                     impute=True,
                                     imputation_bg=imputation_bg_r,
                                     bg_count_name='score')
        cis_contacts_num = cis_coverage['count'].sum()
        cis_coverage['raw_bg_prob'] = cis_coverage['bg_prob']

        chrom_spls = scaling_dict_r.get(gene_chrom)
        if chrom_spls is None:
            pba.update.remote(1)
            return None
        spl = chrom_spls.get(rna_name)
        if spl is None or any(math.isnan(item) for record in spl[:2] for item in record):
            spl = chrom_spls.get(gene_chrom)
            if spl is None:
                pba.update.remote(1)
                return None

        rel_cis_bin_centers = calculate_rel_dist_from_centers(cis_coverage, gene_start, gene_end)
        rel_cis_bin_starts = make_rel_dist_vector(cis_coverage['start'].astype('int'), gene_start, gene_end)
        rel_cis_bin_ends = make_rel_dist_vector(cis_coverage['end'].astype('int'), gene_start, gene_end)
        scaling_factors = 10 ** ((si.splev(np.log10(rel_cis_bin_starts + 1), spl, ext=3)
                                  + si.splev(np.log10(rel_cis_bin_ends + 1), spl, ext=3)
                                  + si.splev(np.log10(rel_cis_bin_centers + 1), spl, ext=3)) / 3)
        # scaling_factors = 10 ** si.splev(np.log10(rel_cis_bin_centers), spl, ext=3)
        cis_coverage['scaling_factor'] = scaling_factors
        cis_scaled_bg_probs = cis_coverage['raw_bg_prob'] * cis_coverage['scaling_factor']
        cis_coverage['bg_prob'] = cis_scaled_bg_probs / cis_scaled_bg_probs.sum()

        # Probs renormalization
        total_contacts_num = trans_contacts_num + cis_contacts_num
        cis_share = cis_contacts_num / total_contacts_num
        trans_share = trans_contacts_num / total_contacts_num

        trans_coverage['raw_bg_prob'] *= trans_share
        trans_coverage['bg_prob'] *= trans_share
        trans_coverage['signal_prob'] *= trans_share

        cis_coverage['raw_bg_prob'] *= cis_share
        cis_coverage['bg_prob'] *= cis_share
        cis_coverage['signal_prob'] *= cis_share

        total_coverage = pd.concat((trans_coverage, cis_coverage), ignore_index=True)
        cleaned_coverage = total_coverage.query('count > 0').reset_index(drop=True)
        cleaned_coverage['fc'] = (cleaned_coverage['signal_prob'] / cleaned_coverage['bg_prob']).fillna(fill_value)

        # statistics
        pvals = cleaned_coverage.apply(lambda row: ss.binomtest(k=row['count'], n=total_contacts_num, p=row['bg_prob'], alternative='greater').pvalue, axis=1)
        cleaned_coverage['pvalue'] = pvals

        cleaned_coverage['rna_name'] = rna_name
    except Exception:
        print('\n', 'Problem with', rna_name, '\n')
        cleaned_coverage = None
    pba.update.remote(1)
    return cleaned_coverage


peaks_results = process_by_chunks(selected_rnas,
                                  calculate_peaks,
                                  n_tasks,
                                  shared_bin_sizes,
                                  shared_contacts_df,
                                  shared_annotation_dict,
                                  shared_chromsizes_df,
                                  shared_chromsizes_dict,
                                  shared_background_df,
                                  shared_scaling_dict,
                                  imputation_bg,
                                  cis_start_size)

concat_peaks = pd.concat(peaks_results, ignore_index=True)
dummy1, qvals, *dummy2 = sm.stats.multipletests(concat_peaks['pvalue'].values, method='fdr_bh')
concat_peaks['qvalue'] = qvals
concat_peaks.to_csv(output_path, sep='\t', header=True, index=False)
