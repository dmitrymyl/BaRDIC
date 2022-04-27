import argparse
import json
from math import ceil, isnan

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.interpolate as si
import scipy.stats as ss
import seaborn as sns
from natsort import natsorted
from tqdm import tqdm


parser = argparse.ArgumentParser(description='Calculate scaling')
input_group = parser.add_argument_group('Input')
input_group.add_argument('ciscoverage',
                         nargs='?',
                         type=str,
                         help='cis coverage tsv file')
param_group = parser.add_argument_group('Parameters')
param_group.add_argument('-degree',
                         nargs='?',
                         type=int,
                         default=3,
                         help='spline degree')
param_group.add_argument('--noRefine',
                         action='store_true',
                         help='include to not refine splines')
param_group.add_argument('-maxThreshold',
                         nargs='?',
                         type=float,
                         default=0.05,
                         help='maximal p-value threshold for spline refinement')
output_group = parser.add_argument_group('Output')
output_group.add_argument('density',
                          nargs='?',
                          type=str,
                          help='png output of raw scaling')
output_group.add_argument('fitted',
                          nargs='?',
                          type=str,
                          help='png output of fitted scaling')
output_group.add_argument('linreg',
                          nargs='?',
                          type=str,
                          help='linreg tsv output')
output_group.add_argument('splines',
                          nargs='?',
                          type=str,
                          help='splines json output')
output_group.add_argument('indsplines',
                          nargs='?',
                          type=str,
                          help='individual splines png output')

args = parser.parse_args()
ciscoverage_filename = args.ciscoverage
density_filename = args.density
density_fitted_filename = args.fitted
linreg_filename = args.linreg
splines_filename = args.splines
indsplines_filename = args.indsplines
degree = args.degree
no_refine = args.noRefine
max_threshold = args.maxThreshold

cis_coverage_data = pd.read_csv(ciscoverage_filename, header=0, index_col=None, sep='\t')

with sns.axes_style(rc={'figure.facecolor': 'white'}):
    g = sns.FacetGrid(data=cis_coverage_data, col='chrom', col_wrap=4, height=4, col_order=natsorted(cis_coverage_data['chrom'].unique()))
    g.map(sns.histplot, 'log_gene_dist', 'log_fc', cbar=True)
    g.set_titles("{col_name}")
    g.set_xlabels('Log10 distance from gene')
    g.set_ylabels('Log10 fc')
    g.savefig(density_filename, dpi=300)

linreg_data = dict()
spline_data = dict()
chrom_counts = cis_coverage_data['chrom'].value_counts()
fittable_chroms = natsorted(chrom_counts[chrom_counts > degree].index)
fittable_chroms_amount = len(fittable_chroms)

for chrom in tqdm(fittable_chroms, desc='Fitting'):
    chrom_data = cis_coverage_data.query('chrom == @chrom').sort_values('log_gene_dist', ignore_index=True)
    x = chrom_data['log_gene_dist'].values
    y = chrom_data['log_fc'].values
    linreg_results = ss.linregress(x, y)
    spline_results = si.splrep(x, y,
                               w=np.ones(x.shape),
                               k=degree)
    linreg_data[chrom] = linreg_results
    spline_data[chrom] = dict()
    spline_data[chrom][chrom] = (list(spline_results[0]), list(spline_results[1]), spline_results[2])
    for rna_name, rna_data in chrom_data.groupby('rna'):
        x = rna_data['log_gene_dist'].values
        y = rna_data['log_fc'].values
        if len(x) > degree:
            spl = si.splrep(x, y, k=degree, w=np.ones(x.shape))
            if not any(isnan(item) for record in spl[:2] for item in record) and not no_refine:
                scaling_factors = 10 ** si.splev(x, spl)
                cis_scaled_prob = rna_data['bg_prob'] * scaling_factors
                rna_data['scaled_bg_prob'] = cis_scaled_prob / cis_scaled_prob.sum()
                n_total = rna_data['count'].sum()
                pvals = rna_data.apply(lambda row: ss.binomtest(k=row['count'],
                                                                n=n_total,
                                                                p=row['scaled_bg_prob'],
                                                                alternative='greater').pvalue,
                                    axis=1)
                threshold = min(1 / len(x), max_threshold)
                rna_data['outlier'] = (pvals < (1 / len(x)))
                cleaned_rna_data = rna_data.query('~outlier')

                if cleaned_rna_data.shape[0] > degree:
                    x = cleaned_rna_data['log_gene_dist'].values
                    y = cleaned_rna_data['log_fc'].values
                    spl = si.splrep(x, y, k=degree, w=np.ones(x.shape))

            spline_data[chrom][rna_name] = (list(spl[0]), list(spl[1]), spl[2])
        else:
            spline_data[chrom][rna_name] = None

cols = 4
width = 6
rows = ceil(fittable_chroms_amount / cols)
plt.figure(figsize=(cols * width, rows * width), facecolor='white')
for i, chrom in tqdm(enumerate(fittable_chroms), desc='Second plot', total=fittable_chroms_amount):
    chrom_data = cis_coverage_data.query('chrom == @chrom').sort_values('log_gene_dist', ignore_index=True)
    x = chrom_data['log_gene_dist'].values
    y = chrom_data['log_fc'].values

    plt.subplot(rows, cols, i + 1)
    sns.histplot(data=chrom_data, x='log_gene_dist', y='log_fc', cbar=True)
    plt.plot(x, x * linreg_data[chrom].slope + linreg_data[chrom].intercept,
             color='red',
             alpha=.7,
             label='LinReg')
    plt.plot(x, si.splev(x, spline_data[chrom][chrom]), 
             color='yellow',
             alpha=.7,
             label='Spline')
    plt.title(chrom)
    plt.xlabel('Log10 distance from gene')
    plt.ylabel('Log10 fc')
    plt.legend(title='Fitting')
plt.tight_layout()
plt.savefig(density_fitted_filename, dpi=300)

with open(linreg_filename, 'w') as outfile:
    outfile.write('Chrom\tSlope\tIntercept\n')
    for chrom, result in linreg_data.items():
        outfile.write(f"{chrom}\t{result.slope}\t{result.intercept}\n")

with open(splines_filename, 'w') as outfile:
    json.dump(spline_data, outfile)

plt.figure(figsize=(cols * width, rows * width), facecolor='white')
for i, chrom in tqdm(enumerate(fittable_chroms), desc='Third plot', total=fittable_chroms_amount):
    chrom_data = cis_coverage_data.query('chrom == @chrom').sort_values('log_gene_dist', ignore_index=True)
    x = chrom_data['log_gene_dist'].values
    sampled_x = np.linspace(x.min(), x.max(), 50)
    predicted_for_chrom = list()
    for rna_name, rna_data in chrom_data.groupby('rna'):
        if rna_data.shape[0] > degree:
            predicted = si.splev(sampled_x, spline_data[chrom][rna_name])
            predicted_for_chrom.append(pd.DataFrame({'x': sampled_x, 'pred': predicted, 'rna': rna_name}))
    predicted_df = pd.concat(predicted_for_chrom, ignore_index=True)
    plt.subplot(rows, cols, i + 1)
    sns.lineplot(data=predicted_df,
                 x='x',
                 y='pred',
                 hue='rna',
                 ci=None,
                 palette='gnuplot2',
                 legend=False,
                 alpha=.2)
    plt.plot(x,
             si.splev(x, spline_data[chrom][chrom]), 
             color='black',
             alpha=.7,
             label='general',
             lw=2,
             marker=None)
    plt.title(chrom)
    plt.xlabel('Log10 distance from gene')
    plt.ylabel('Log10 fc')
    plt.yscale('symlog')
    plt.legend(loc='lower right')
plt.tight_layout()
plt.savefig(indsplines_filename, dpi=300)
