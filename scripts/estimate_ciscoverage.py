from collections import namedtuple

import numpy as np
import pandas as pd
import ray
import argparse
import os


from bardic.binops import make_cis_bins3, make_interval_centers, make_rel_dist_vector, compute_track
from bardic.io import load_dna_parts_of_rdc2, load_bedgraph
from bardic.schemas import rdc_gene_dtypes, rdc_gene_schema
from bardic.mp import process_by_chunks


os.environ["RAY_OBJECT_STORE_ALLOW_SLOW_STORAGE"] = "1"

parser = argparse.ArgumentParser('Estimate cis coverage v5')
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
input_group.add_argument('background',
                         nargs='?',
                         type=str,
                         help='Binarized background coverage in bed3 format')
parameters_group = parser.add_argument_group('Parameters')
parameters_group.add_argument('-cisStartSize',
                              nargs="?",
                              type=int,
                              default=1000,
                              help='Start cis bin size')
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
args = parser.parse_args()
contacts_path = args.contacts
annotation_path = args.annotation
chromsizes_path = args.chromsizes
selection_path = args.selection
background_path = args.background
cis_start_size = args.cisStartSize
output_path = args.output
n_cores = args.ncores
n_tasks = args.ntasks

Coords = namedtuple('Coords', ['chrom', 'start', 'end'])
contacts_df = load_dna_parts_of_rdc2(contacts_path)
chromsizes_df = pd.read_csv(chromsizes_path,
                            header=None,
                            index_col=None,
                            sep='\t',
                            names=['chrom', 'length'])
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
                                header=[0],
                                index_col=None,
                                sep='\t')

bin_sizes = {item[1]: (int(item[2]), float(item[3]))
             for item in selection_results.loc[:, ['gene_name', 'trans_bin_size', 'cis_factor']].to_records()}

background_track = pd.read_csv(background_path,
                               header=None,
                               index_col=None,
                               sep='\t',
                               names=['chrom', 'start', 'end', 'count'],
                               dtype={'chrom': 'category',
                                      'start': 'int',
                                      'end': 'int',
                                      'count': 'int'})
background_track = load_bedgraph(background_path)
background_track['prob'] = background_track['score'] / background_track['score'].sum()


ray.init(num_cpus=n_cores)

shared_contacts_df = ray.put(contacts_df)
shared_bin_sizes = ray.put(bin_sizes)
shared_annotation_dict = ray.put(gene_coords)
shared_chromsizes_dict = ray.put(chromsizes_dict)
shared_background_df = ray.put(background_track)


@ray.remote
def get_processed_cov_for_scaling(rna_name,
                                  bin_sizes_r,
                                  contacts_df_r,
                                  annotation_dict_r,
                                  chromsizes_dict_r,
                                  background_df_r,
                                  cis_start_size,
                                  pba):
    gene_annotation = annotation_dict_r.get(rna_name)
    if gene_annotation is None:
        pba.update.remote(1)
        return None
    gene_start = gene_annotation.start
    gene_end = gene_annotation.end
    gene_chrom = gene_annotation.chrom
    chrom_size = chromsizes_dict_r.get(gene_chrom)
    if chrom_size is None:
        pba.update.remote(1)
        return None
    trans_bin_size, cis_factor = bin_sizes_r.get(rna_name, (None, None))
    if cis_factor is None:
        pba.update.remote(1)
        return None
    cis_bins = make_cis_bins3(cis_factor,
                              start_size=cis_start_size,
                              chrom_name=gene_chrom,
                              chrom_length=chrom_size,
                              gene_start=gene_start,
                              gene_end=gene_end,
                              max_linear_size=trans_bin_size)
    rna_contacts = contacts_df_r.query("name == @rna_name").reset_index(drop=True)
    rna_contacts_centers = make_interval_centers(rna_contacts)
    cis_coverage = compute_track(cis_bins, rna_contacts_centers, background_df_r, bg_count_name='score')
    undefined = (cis_coverage['count'] == 0) | (cis_coverage['bg_count'] == 0)
    cis_coverage['fc'] = cis_coverage['signal_prob'] / cis_coverage['bg_prob']
    cis_coverage = cis_coverage.loc[~undefined, :].reset_index(drop=True)

    bin_centers = (cis_coverage['start'] + cis_coverage['end']) // 2
    gene_dist = make_rel_dist_vector(bin_centers, gene_start, gene_end)
    cis_coverage['log_gene_dist'] = np.log10(gene_dist)
    cis_coverage['log_fc'] = np.log10(cis_coverage['fc']).astype('float64')

    cis_coverage = cis_coverage.query('count > 0 and bg_count > 0').loc[:, ['log_gene_dist', 'log_fc', 'count', 'bg_prob']].reset_index(drop=True)
    cis_coverage['chrom'] = gene_chrom
    cis_coverage['rna'] = rna_name
    pba.update.remote(1)
    return cis_coverage


selected_rnas = [rna_name for rna_name in bin_sizes if rna_name in contacts_df['name'].values]

scaling_results = process_by_chunks(selected_rnas,
                                    get_processed_cov_for_scaling,
                                    n_tasks,
                                    shared_bin_sizes,
                                    shared_contacts_df,
                                    shared_annotation_dict,
                                    shared_chromsizes_dict,
                                    shared_background_df,
                                    cis_start_size)

concat_scaling = pd.concat(scaling_results, ignore_index=True)
concat_scaling.to_csv(output_path, sep='\t', header=True, index=False)
