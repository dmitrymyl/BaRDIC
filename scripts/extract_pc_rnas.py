import pandas as pd
import argparse


parser = argparse.ArgumentParser(description='extract names of protein-coding RNAs')
input_group = parser.add_argument_group('Input')
input_group.add_argument('input',
                         nargs='?',
                         type=str,
                         help='input')
output_group = parser.add_argument_group('Output')
output_group.add_argument('output',
                          nargs='?',
                          type=str,
                          help='output')
args = parser.parse_args()
input_file = args.input
output_file = args.output

contacts_stats = pd.read_csv(input_file, sep='\t', header=0, index_col=None)
contacts_stats.query('coding == "coding"')['gene_name'].to_csv(output_file, sep='\t', header=False, index=False)
