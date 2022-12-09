import argparse


bardic_parser = argparse.ArgumentParser(prog='bardic',
                                        description='Binomial RNA-DNA interaction caller.')
bardic_subparsers = bardic_parser.add_subparsers(title='Subcommands',
                                                 metavar='SUBCOMMAND',
                                                 required=True)


bed2h5_parser = bardic_subparsers.add_parser('bed2h5',
                                             help='convert DNA parts to custom HDF5 format',
                                             description='Convert DNA parts to custom HDF5 format.')

bed2h5_input_group = bed2h5_parser.add_argument_group('Input')
bed2h5_input_group.add_argument('dnaparts',
                                type=str,
                                nargs='?',
                                help='BED6 file with coordinates of DNA parts. Names of corresponding RNAs are in the "name" column.')
bed2h5_input_group.add_argument('annotation',
                                type=str,
                                nargs='?',
                                help='RNA annotation in BED format.')
bed2h5_input_group.add_argument('chromsizes',
                                type=str,
                                nargs='?',
                                help='If filename, then it is a UCSC headerless chromsizes file; if genome abbreviation, then will fetch chromsizes from UCSC')

bed2h5_output_group = bed2h5_parser.add_argument_group('Output')
bed2h5_output_group.add_argument('output',
                                 type=str,
                                 nargs='?',
                                 help='output file in dnah5 format')


binsizes_parser = bardic_subparsers.add_parser('binsizes',
                                               help='Select bin size for each RNA and save it into dnah5 file',
                                               description='Select bin size for each RNA and save it into dnah5 file',
                                               formatter_class=argparse.ArgumentDefaultsHelpFormatter)

binsizes_input_group = binsizes_parser.add_argument_group('Input')
binsizes_input_group.add_argument('dnadataset',
                                  type=str,
                                  nargs='?',
                                  help='DNA parts dnah5 file.')

binsizes_params_group = binsizes_parser.add_argument_group('Parameters')
binsizes_params_group.add_argument('-min_contacts',
                                   type=int,
                                   nargs='?',
                                   default=1000,
                                   help='Minimal number of contacts to consider an RNA. Any RNA with less contacts will be discarded from further processing.')
binsizes_params_group.add_argument('-trans_min',
                                   type=int,
                                   nargs='?',
                                   default=10_000,
                                   help='Minimal trans bin size.')
binsizes_params_group.add_argument('-trans_max',
                                   type=int,
                                   nargs='?',
                                   default=1_000_000,
                                   help='Maximal trans bin size.')
binsizes_params_group.add_argument('-trans_step',
                                   type=int,
                                   nargs='?',
                                   default=1_000,
                                   help='Step for increasing trans bin size.')
binsizes_params_group.add_argument('-cis_min',
                                   type=float,
                                   nargs='?',
                                   default=1.1,
                                   help='Minimal cis factor.')
binsizes_params_group.add_argument('-cis_max',
                                   type=float,
                                   nargs='?',
                                   default=2.,
                                   help='Maximal cis factor.')
binsizes_params_group.add_argument('-cis_step',
                                   type=float,
                                   nargs='?',
                                   default=0.01,
                                   help='Step for inreasing cis factor.')
binsizes_params_group.add_argument('-cis_start',
                                   type=int,
                                   nargs='?',
                                   default=5000,
                                   help='Starting cis bin size.')
binsizes_params_group.add_argument('-tolerance',
                                   type=float,
                                   nargs='?',
                                   default=0.01,
                                   help='Maximal absolute difference between two consecutive cost function values to consider optimization converged.')
binsizes_params_group.add_argument('-w',
                                   type=float,
                                   nargs='?',
                                   default=1,
                                   help='Window size to average cost function values over.')

binsizes_output_group = binsizes_parser.add_argument_group('Output')
binsizes_output_group.add_argument('output',
                                   type=str,
                                   nargs='?',
                                   help='Output tsv table with extended bin size selection results.')

binsizes_processing_group = binsizes_parser.add_argument_group('Processing')
binsizes_processing_group.add_argument('-cores',
                                       type=int,
                                       nargs='?',
                                       default=1,
                                       help='Maximal number of cores to use.')
