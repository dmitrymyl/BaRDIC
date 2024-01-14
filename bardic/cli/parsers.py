import argparse
from typing import Union

from .commands import (background_cli, bed2h5_cli, binsizes_cli, makerdc_cli,
                       peaks_cli, scaling_cli, run_pipeline_cli)


class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.MetavarTypeHelpFormatter):
    def _get_default_metavar_for_positional(self, action):
        return action.dest


def numeric(value: str) -> Union[int, float]:
    transformed_val: float = float(value)
    if transformed_val.is_integer():
        transformed_val = int(transformed_val)
    return transformed_val


bardic_parser = argparse.ArgumentParser(prog='bardic',
                                        description='Binomial RNA-DNA interaction caller.',
                                        formatter_class=CustomFormatter,
                                        fromfile_prefix_chars='@')
bardic_subparsers = bardic_parser.add_subparsers(title='Subcommands',
                                                 metavar='SUBCOMMAND',
                                                 required=True)


bed2h5_parser = bardic_subparsers.add_parser('bed2h5',
                                             help='convert DNA parts to custom HDF5 format',
                                             description='Convert DNA parts to custom HDF5 format.',
                                             formatter_class=CustomFormatter)
bed2h5_parser.set_defaults(func=bed2h5_cli)

bed2h5_input_group = bed2h5_parser.add_argument_group('Input')
bed2h5_input_group.add_argument('dnaparts',
                                type=str,
                                help='BED6 file with coordinates of DNA parts. Names of corresponding RNAs are in the "name" column.')
bed2h5_input_group.add_argument('annotation',
                                type=str,
                                help='RNA annotation in BED format.')
bed2h5_input_group.add_argument('chromsizes',
                                type=str,
                                help='If filename, then it is a UCSC headerless chromsizes file; if genome abbreviation, then will fetch chromsizes from UCSC')

bed2h5_output_group = bed2h5_parser.add_argument_group('Output')
bed2h5_output_group.add_argument('output',
                                 type=str,
                                 help='output file in dnah5 format')


binsizes_parser = bardic_subparsers.add_parser('binsizes',
                                               help='Select bin size for each RNA and save it into dnah5 file',
                                               description='Select bin size for each RNA and save it into dnah5 file',
                                               formatter_class=CustomFormatter)
binsizes_parser.set_defaults(func=binsizes_cli)

binsizes_input_group = binsizes_parser.add_argument_group('Input')
binsizes_input_group.add_argument('dna_dataset',
                                  type=str,
                                  help='DNA parts dnah5 file.')

binsizes_params_group = binsizes_parser.add_argument_group('Parameters')
binsizes_params_group.add_argument('-mcon', '--min_contacts',
                                   type=int,
                                   nargs='?',
                                   dest='n_contacts',
                                   default=1000,
                                   help='Minimal number of contacts to consider an RNA. Any RNA with less contacts will be discarded from further processing.')
binsizes_params_group.add_argument('-tmin', '--trans_min',
                                   type=int,
                                   nargs='?',
                                   default=10_000,
                                   help='Minimal trans bin size.')
binsizes_params_group.add_argument('-tmax', '--trans_max',
                                   type=int,
                                   nargs='?',
                                   default=1_000_000,
                                   help='Maximal trans bin size.')
binsizes_params_group.add_argument('-tstep', '--trans_step',
                                   type=int,
                                   nargs='?',
                                   default=1_000,
                                   help='Step for increasing trans bin size.')
binsizes_params_group.add_argument('-cmin', '--cis_min',
                                   type=float,
                                   nargs='?',
                                   default=1.1,
                                   help='Minimal cis factor.')
binsizes_params_group.add_argument('-cmax', '--cis_max',
                                   type=float,
                                   nargs='?',
                                   default=2.,
                                   help='Maximal cis factor.')
binsizes_params_group.add_argument('-cstep', '--cis_step',
                                   type=float,
                                   nargs='?',
                                   default=0.01,
                                   help='Step for inreasing cis factor.')
binsizes_params_group.add_argument('-cstart', '--cis_start',
                                   type=int,
                                   nargs='?',
                                   default=5000,
                                   help='Starting cis bin size.')
binsizes_params_group.add_argument('-tol', '--tolerance',
                                   type=float,
                                   nargs='?',
                                   default=0.01,
                                   help='Maximal absolute difference between two consecutive cost function values to consider optimization converged.')
binsizes_params_group.add_argument('-w', '--window',
                                   type=float,
                                   nargs='?',
                                   dest='w',
                                   default=1,
                                   help='Window size to average cost function values over.')

binsizes_output_group = binsizes_parser.add_argument_group('Output')
binsizes_output_group.add_argument('output',
                                   type=str,
                                   help='Output tsv table with extended bin size selection results.')

binsizes_processing_group = binsizes_parser.add_argument_group('Processing')
binsizes_processing_group.add_argument('-c', '--cores',
                                       type=int,
                                       nargs='?',
                                       dest='n_cores',
                                       default=1,
                                       help='Maximal number of cores to use.')


background_parser = bardic_subparsers.add_parser('background',
                                                 help='Create a bedGraph background track from DNA parts of selected RNAs',
                                                 description='Create a bedGraph background track from DNA parts of selected RNAs',
                                                 formatter_class=CustomFormatter)
background_parser.set_defaults(func=background_cli)

background_input_group = background_parser.add_argument_group('Input')
background_input_group.add_argument('dna_dataset',
                                    type=str,
                                    help='DNA parts dnah5 file.')
background_input_group.add_argument('rnas',
                                    type=str,
                                    help='A file with a list of RNAs with one RNA name per line.')

background_params_group = background_parser.add_argument_group('Parameters')
background_params_group.add_argument('-bs', '--binsize',
                                     type=int,
                                     nargs='?',
                                     default=1000,
                                     help='Bin size of the background track.')

background_output_group = background_parser.add_argument_group('Output')
background_output_group.add_argument('output',
                                     type=str,
                                     help='Filename of a background track in a bedGraph format.')


makerdc_parser = bardic_subparsers.add_parser('makerdc',
                                              help='Create RDC file from dnah5 DNA parts and bedGraph background track.',
                                              description='Create RDC file from dnah5 DNA parts and bedGraph background track.',
                                              formatter_class=CustomFormatter)
makerdc_parser.set_defaults(func=makerdc_cli)

makerdc_input_group = makerdc_parser.add_argument_group('Input')
makerdc_input_group.add_argument('dna_dataset',
                                 type=str,
                                 help='dnah5 DNA parts file.')
makerdc_input_group.add_argument('bg',
                                 type=str,
                                 help='bedGraph background track.')

makerdc_params_group = makerdc_parser.add_argument_group('Parameters')
makerdc_params_group.add_argument('-i', '--ifactor',
                                  type=float,
                                  nargs='?',
                                  default=0.01,
                                  help='Imputation factor: if background coverage of a bin is 0, this value is a multiplier of an average background coverage to impute zero background coverage.')

makerdc_output_group = makerdc_parser.add_argument_group('Output')
makerdc_output_group.add_argument('output',
                                  type=str,
                                  help='Output .rdc filename.')

makerdc_processing_group = makerdc_parser.add_argument_group('Processing')
makerdc_processing_group.add_argument('-c', '--cores',
                                      type=int,
                                      nargs='?',
                                      default=1,
                                      dest='n_cores',
                                      help='Maximal number of cores to use.')


scaling_parser = bardic_subparsers.add_parser('scaling',
                                              help='Estimate scaling by fitting splines and adjust background probabilities in RDC file.',
                                              description='Estimate scaling by fitting splines and adjust background probabilities in RDC file.',
                                              formatter_class=CustomFormatter)
scaling_parser.set_defaults(func=scaling_cli)

scaling_input_group = scaling_parser.add_argument_group('Input')
scaling_input_group.add_argument('rdc',
                                 type=str,
                                 help='Input RDC file name.')

scaling_params_group = scaling_parser.add_argument_group('Parameters')
scaling_params_group.add_argument('-d', '--degree',
                                  type=int,
                                  nargs='?',
                                  default=3,
                                  help='Spline degree.')
scaling_params_group.add_argument('-mt', '--max_threshold',
                                  type=float,
                                  nargs='?',
                                  default=0.05,
                                  help='Maximal binomial test p-value to consider a point as an outlier in a spline refinement procedure.')
scaling_params_group.add_argument('-ns', '--no_scaling',
                                  action='store_true',
                                  help='If included, do not estimate scaling.')
scaling_params_group.add_argument('-nr', '--no_refine',
                                  action='store_true',
                                  help='If included, do not apply a spline refinement procedure.')
scaling_params_group.add_argument('-fv', '--fill_value',
                                  type=numeric,
                                  nargs='?',
                                  default=1,
                                  help='Fold-change fill ratio in case of 0/0.')

scaling_processing_group = scaling_parser.add_argument_group('Processing')
scaling_processing_group.add_argument('-c', '--cores',
                                      type=int,
                                      nargs='?',
                                      default=1,
                                      dest='n_cores',
                                      help='Maximal number of cores to use.')


valid_score_fields = ("bg_count",
                      "raw_bg_prob",
                      "scaling_factor",
                      "bg_prob",
                      "signal_count",
                      "signal_prob",
                      "impute",
                      "fc",
                      "pvalue",
                      "qvalue")


def score_field(value):
    try:
        transformed_val = int(value)
    except ValueError:
        transformed_val = str(value)
        if transformed_val not in valid_score_fields:
            raise ValueError(f'Provided score field {value} is not one of valid RDC fields ({", ".join(valid_score_fields)})')
    return transformed_val


peaks_parser = bardic_subparsers.add_parser('peaks',
                                            help='Estimate significance and fetch peaks at specified FDR level.',
                                            description='Estimate significance and fetch peaks at specified FDR level. Significance is estimated only once and all p- and q-values are store in the RDC.',
                                            formatter_class=CustomFormatter)
peaks_parser.set_defaults(func=peaks_cli)

peaks_input_group = peaks_parser.add_argument_group('Input')
peaks_input_group.add_argument('rdc',
                               type=str,
                               help='RDC filename.')

peaks_params_group = peaks_parser.add_argument_group('Parameters')
peaks_params_group.add_argument('-q', '--qval_threshold',
                                type=float,
                                nargs='?',
                                default=0.05,
                                help='BH q-value threshold to consider bin a peak.')

peaks_params_group.add_argument('-qt', '--qval_type',
                                type=str,
                                nargs='?',
                                choices=['global', 'rna'],
                                default='global',
                                help='BH q-value type to use for peak calling. '
                                     'If "global" (default), will use q-values calculated for all RNAs; '
                                     'if "rna", will use q-values calculated for each RNA separately.')

peaks_output_group = peaks_parser.add_argument_group('Output')
peaks_output_group.add_argument('output',
                                type=str,
                                help='Output peaks filename.')
peaks_output_group.add_argument('-f', '--format',
                                type=str,
                                nargs='?',
                                choices=['narrowPeak', 'bed'],
                                default='narrowPeak',
                                help='Output peaks file format.')
peaks_output_group.add_argument('-s', '--score',
                                type=score_field,
                                nargs='?',
                                default=0,
                                help='If --format=bed, which value to fill the score field with. '
                                     'If int, will fill every peak score with it; '
                                     'if str, will take corresponding values from the column in RDC '
                                     f'(choices: {", ".join(valid_score_fields)})')

peaks_processing_group = peaks_parser.add_argument_group('Processing')
peaks_processing_group.add_argument('-c', '--cores',
                                    type=int,
                                    nargs='?',
                                    default=1,
                                    dest='n_cores',
                                    help='Maximal number of cores to use.')


run_pipeline_parser = bardic_subparsers.add_parser('run',
                                                   help='Run pipeline with a single command.',
                                                   description='Run pipeline with a single command.',
                                                   formatter_class=CustomFormatter)
run_pipeline_parser.set_defaults(func=run_pipeline_cli)

run_input_group = run_pipeline_parser.add_argument_group('Input')
run_input_group.add_argument('dnaparts',
                             type=str,
                             help='BED6 file with coordinates of DNA parts. Names of corresponding RNAs are in the "name" column.')
run_input_group.add_argument('annotation',
                             type=str,
                             help='RNA annotation in BED format.')
run_input_group.add_argument('chromsizes',
                             type=str,
                             help='If filename, then it is a UCSC headerless chromsizes file; if genome abbreviation, then will fetch chromsizes from UCSC')
run_input_group.add_argument('bgdata',
                             type=str,
                             help='A file with data on background. If --bgtype="rnas", '
                                  'this is a file with a list of RNAs with one RNA name per line. '
                                  'If --bgtype="custom", this is a bedGraph file with '
                                  'background signal in equally-sized bins.')

run_output_group = run_pipeline_parser.add_argument_group('Output')
run_output_group.add_argument('outdir',
                              type=str,
                              help='Output directory name.')
run_output_group.add_argument('-f', '--format',
                              type=str,
                              nargs='?',
                              choices=['narrowPeak', 'bed'],
                              default='narrowPeak',
                              help='Output peaks file format.')
run_output_group.add_argument('-s', '--score',
                              type=score_field,
                              nargs='?',
                              default=0,
                              help='If --format=bed, which value to fill the score field with. '
                                   'If int, will fill every peak score with it; '
                                   'if str, will take corresponding values from the column in RDC '
                                   f'(choices: {", ".join(valid_score_fields)})')

run_binsizes_group = run_pipeline_parser.add_argument_group('Binsize selection parameters')
run_binsizes_group.add_argument('-mcon', '--min_contacts',
                                type=int,
                                nargs='?',
                                dest='n_contacts',
                                default=1000,
                                help='Minimal number of contacts to consider an RNA. Any RNA with less contacts will be discarded from further processing.')
run_binsizes_group.add_argument('-tmin', '--trans_min',
                                type=int,
                                nargs='?',
                                default=10_000,
                                help='Minimal trans bin size.')
run_binsizes_group.add_argument('-tmax', '--trans_max',
                                type=int,
                                nargs='?',
                                default=1_000_000,
                                help='Maximal trans bin size.')
run_binsizes_group.add_argument('-tstep', '--trans_step',
                                type=int,
                                nargs='?',
                                default=1_000,
                                help='Step for increasing trans bin size.')
run_binsizes_group.add_argument('-cmin', '--cis_min',
                                type=float,
                                nargs='?',
                                default=1.1,
                                help='Minimal cis factor.')
run_binsizes_group.add_argument('-cmax', '--cis_max',
                                type=float,
                                nargs='?',
                                default=2.,
                                help='Maximal cis factor.')
run_binsizes_group.add_argument('-cstep', '--cis_step',
                                type=float,
                                nargs='?',
                                default=0.01,
                                help='Step for inreasing cis factor.')
run_binsizes_group.add_argument('-cstart', '--cis_start',
                                type=int,
                                nargs='?',
                                default=5000,
                                help='Starting cis bin size.')
run_binsizes_group.add_argument('-tol', '--tolerance',
                                type=float,
                                nargs='?',
                                default=0.01,
                                help='Maximal absolute difference between two consecutive cost function values to consider optimization converged.')
run_binsizes_group.add_argument('-w', '--window',
                                type=float,
                                nargs='?',
                                dest='w',
                                default=1,
                                help='Window size to average cost function values over.')

run_background_group = run_pipeline_parser.add_argument_group('Background parameters')
run_background_group.add_argument('-bs', '--binsize',
                                  type=int,
                                  nargs='?',
                                  default=1000,
                                  help='Bin size of the background track.')
run_background_group.add_argument('-bt', '--bgtype',
                                  type=str,
                                  nargs='?',
                                  dest='bg_type',
                                  default='rnas',
                                  choices=('rnas', 'custom'),
                                  help='Type of backround. If "rnas", then will calculate background '
                                       'from trans-contacts of RNAs supplied as "bgdata". If "custom", '
                                       'will use bedgraph track provided as "bgdata".')

run_makerdc_group = run_pipeline_parser.add_argument_group('RDC creation parameters')
run_makerdc_group.add_argument('-i', '--ifactor',
                               type=float,
                               nargs='?',
                               default=0.01,
                               help='Imputation factor: if background coverage of a bin is 0, this value is a multiplier of an average background coverage to impute zero background coverage.')

run_scaling_group = run_pipeline_parser.add_argument_group('Scaling parameters')
run_scaling_group.add_argument('-d', '--degree',
                               type=int,
                               nargs='?',
                               default=3,
                               help='Spline degree.')
run_scaling_group.add_argument('-mt', '--max_threshold',
                               type=float,
                               nargs='?',
                               default=0.05,
                               help='Maximal binomial test p-value to consider a point as an outlier in a spline refinement procedure.')
run_scaling_group.add_argument('-ns', '--no_scaling',
                               action='store_true',
                               help='If included, do not estimate scaling.')
run_scaling_group.add_argument('-nr', '--no_refine',
                               action='store_true',
                               help='If included, do not apply a spline refinement procedure.')
run_scaling_group.add_argument('-fv', '--fill_value',
                               type=numeric,
                               nargs='?',
                               default=1,
                               help='Fold-change fill ratio in case of 0/0.')

run_peaks_group = run_pipeline_parser.add_argument_group('Peaks parameters')
run_peaks_group.add_argument('-q', '--qval_threshold',
                             type=float,
                             nargs='?',
                             default=0.05,
                             help='BH q-value threshold to consider bin a peak.')
run_peaks_group.add_argument('-qt', '--qval_type',
                             type=str,
                             nargs='?',
                             choices=['global', 'rna'],
                             default='global',
                             help='BH q-value type to use for peak calling. '
                                  'If "global" (default), will use q-values calculated for all RNAs; '
                                  'if "rna", will use q-values calculated for each RNA separately.')

run_processing_group = run_pipeline_parser.add_argument_group('Processing')
run_processing_group.add_argument('-c', '--cores',
                                  type=int,
                                  nargs='?',
                                  default=1,
                                  dest='n_cores',
                                  help='Maximal number of cores to use.')
