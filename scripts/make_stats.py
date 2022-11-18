import argparse
import pandas as pd
from bardic.api.io import load_rdc2
from bardic.api.schemas import rdc_gene_dtypes, rdc_gene_schema


parser = argparse.ArgumentParser(description='make contacts stats')
input_group = parser.add_argument_group('Input')
input_group.add_argument('contacts',
                         nargs='?',
                         type=str,
                         help='rdc contacts file')
input_group.add_argument('annotation',
                         nargs='?',
                         type=str,
                         help='Annotation file')
output_group = parser.add_argument_group('Output')
output_group.add_argument('output',
                          nargs='?',
                          type=str,
                          help='output file')

args = parser.parse_args()
contacts_file = args.contacts
annotation_file = args.annotation
output_file = args.output

contacts_data = load_rdc2(contacts_file)
annotation = pd.read_csv(annotation_file, header=None, index_col=None, sep='\t', names=rdc_gene_schema, dtype=rdc_gene_dtypes)

trans_contacts_stats = contacts_data['gene_name'][contacts_data['rna_chr'] != contacts_data['dna_chr']].value_counts()
cis_contacts_stats = contacts_data['gene_name'][contacts_data['rna_chr'] == contacts_data['dna_chr']].value_counts()
contacts_stats = pd.DataFrame({'trans': trans_contacts_stats, 'cis': cis_contacts_stats})
contacts_stats['total'] = contacts_stats['cis'] + contacts_stats['trans']
contacts_stats = contacts_stats.reset_index().rename(columns={'index': 'gene_name'})

annotation['coding'] = annotation['gene_type'].map(lambda rna: 'coding' if rna == 'protein_coding' else 'noncoding')
stats_w_biotype = pd.merge(contacts_stats, annotation.loc[:, ['gene_name', 'gene_type', 'coding']], how='left', on='gene_name')

stats_w_biotype.to_csv(output_file, sep='\t', header=True, index=False)
