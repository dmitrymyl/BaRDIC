import pandas as pd
from bardic.io import load_rdc2
from bardic.schemas import rdc_to_dna
from bardic.binops import make_linear_bins, make_interval_centers, calculate_bins_coverage
import argparse


parser = argparse.ArgumentParser(description='Makes bg track from contacts and given RNAs')
input_group = parser.add_argument_group('Input')
input_group.add_argument('contacts',
                         nargs='?',
                         type=str,
                         help='rdc file')
input_group.add_argument('rnas',
                         nargs='?',
                         type=str,
                         help='file with bg RNAs')
input_group.add_argument('chromsizes',
                         nargs='?',
                         type=str,
                         help='chromsizes file')
params_group = parser.add_argument_group('Parameters')
params_group.add_argument('-binsize',
                          nargs='?',
                          type=int,
                          default=1000,
                          help='bin size')
output_group = parser.add_argument_group('Output')
output_group.add_argument('output',
                          nargs='?',
                          type=str,
                          help='output bedgraph file')

args = parser.parse_args()
contacts_file = args.contacts
rnas_file = args.rnas
chromsizes_file = args.chromsizes
binsize = args.binsize
output_file = args.output

contacts = load_rdc2(contacts_file)
rnas = pd.read_csv(rnas_file, sep='\t', header=None, index_col=None)[0]
chromsizes_df = pd.read_csv(chromsizes_file, header=None, index_col=None, sep='\t', names=['chrom', 'length']).query('chrom != "chrM"').reset_index(drop=True)

bins_df = make_linear_bins(binsize, chromsizes_df)
trans_dna_parts = contacts[contacts['rna_chr'] != contacts['dna_chr']].loc[:, rdc_to_dna.keys()].rename(columns=rdc_to_dna).reset_index(drop=True)
contacts_centers = make_interval_centers(trans_dna_parts.query('name in @rnas'))
bins_coverage = calculate_bins_coverage(bins_df, contacts_centers)
bins_coverage.to_csv(output_file, sep='\t', header=False, index=False)
