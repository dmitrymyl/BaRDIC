import argparse
import os
from collections import namedtuple

import pandas as pd
import ray
from bardic.api.io import load_dna_parts_of_rdc2
from bardic.api.schemas import rdc_gene_dtypes, rdc_gene_schema
from bardic.utils.binsizes import select_bin_sizes
from bardic.api.convert import encode_chromsizes


parser = argparse.ArgumentParser(description='select bin size for specific RNAs')

input_group = parser.add_argument_group('Input')
input_group.add_argument('contacts',
                         nargs='?',
                         type=str,
                         help='rdc filename')
input_group.add_argument('annotation',
                         nargs='?',
                         type=str,
                         help='RNA annotation, bed')
input_group.add_argument('chromsizes',
                         nargs='?',
                         type=str,
                         help='chromsizes filename')

output_group = parser.add_argument_group('Output')
output_group.add_argument('output',
                          nargs='?',
                          type=str,
                          help='output filename')

filtering_group = parser.add_argument_group('Filtering')
filtering_group.add_argument('-ncontacts',
                    nargs='?',
                    type=int,
                    default=1000,
                    help='Minimal number of RNA\'s contacts to consider it for binning')

common_bin_group = parser.add_argument_group('Common selection parameters')
common_bin_group.add_argument('-tolerance',
                              nargs='?',
                              type=float,
                              default=0.01,
                              help='Tolerance')
common_bin_group.add_argument('-w',
                              nargs='?',
                              type=int,
                              default=1,
                              help='Cost function smoothing window size')

trans_bin_group = parser.add_argument_group('Trans bin selection parameters')
trans_bin_group.add_argument('-minTransBin',
                             nargs='?',
                             type=int,
                             default=10000,
                             help='Minimal trans bin size')
trans_bin_group.add_argument('-maxTransBin',
                             nargs='?',
                             type=int,
                             default=1000000,
                             help='Maximal trans bin size')
trans_bin_group.add_argument('-stepTransBin',
                             nargs='?',
                             type=int,
                             default=1000,
                             help='Trans bin size step')

cis_bin_group = parser.add_argument_group('Cis bin selection parameters')
cis_bin_group.add_argument('-minCisFactor',
                             nargs='?',
                             type=float,
                             default=1.1,
                             help='Minimal cis logarithmic bin size')
cis_bin_group.add_argument('-maxCisFactor',
                             nargs='?',
                             type=float,
                             default=2,
                             help='Maximal cis logarithmic bin size')
cis_bin_group.add_argument('-stepCisFactor',
                             nargs='?',
                             type=float,
                             default=0.01,
                             help='Cis logarithmic bin size step')
cis_bin_group.add_argument('-cisStartSize',
                           nargs='?',
                           type=int,
                           default=5000,
                           help='Cis bins start size')

processing_group = parser.add_argument_group('Processing')
processing_group.add_argument('-ncores',
                              nargs='?',
                              type=int,
                              default=1,
                              help='Number of cores to use')
processing_group.add_argument('-ntasks',
                              nargs='?',
                              type=int,
                              default=1000,
                              help='Number of simultanious tasks to be scheduled')

args = parser.parse_args()

contacts_filename = args.contacts
annot_filename = args.annotation
chromsizes_filename = args.chromsizes
output_filename = args.output

n_cores = args.ncores
n_tasks = args.ntasks
n_contacts = args.ncontacts

tolerance = args.tolerance
w = args.w

trans_min = args.minTransBin
trans_step = args.stepTransBin
trans_max = args.maxTransBin

cis_min = args.minCisFactor
cis_max = args.maxCisFactor
cis_step = args.stepCisFactor
cis_start = args.cisStartSize

contacts_data = load_dna_parts_of_rdc2(contacts_filename)
chromsizes_df = pd.read_csv(chromsizes_filename,
                            header=None,
                            index_col=None,
                            sep='\t',
                            names=['chrom', 'length']).query('chrom != "chrM"').copy().reset_index(drop=True)
rna_annot_df = pd.read_csv(annot_filename,
                           index_col=None,
                           sep='\t',
                           header=None,
                           names=rdc_gene_schema,
                           dtype=rdc_gene_dtypes)

chromsizes_dict = encode_chromsizes(chromsizes_df)
Coords = namedtuple('Coords', ['chrom', 'start', 'end'])
annot_dict = {item[1]: Coords(item[2], item[3], item[4])
              for item in rna_annot_df.loc[:, ['gene_name', 'gene_chr', 'gene_start', 'gene_end']].to_records()}


selection_df = select_bin_sizes(contacts_data,
                                chromsizes_dict,
                                annot_dict,
                                n_contacts=n_contacts,
                                trans_min=trans_min,
                                trans_max=trans_max,
                                trans_step=trans_step,
                                cis_min=cis_min,
                                cis_max=cis_max,
                                cis_step=cis_step,
                                cis_start=cis_start,
                                tolerance=tolerance,
                                w=w,
                                n_cores=n_cores,
                                n_tasks=n_tasks)

selection_df.to_csv(output_filename, header=True, index=False, sep='\t')
