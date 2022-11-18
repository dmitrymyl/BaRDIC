import pandas as pd

from .schemas import bedgraph_schema, bedgraph_dtypes


def load_bedgraph(bg_filename):
    """Loads bedgraph file with specific dtypes into pd.DataFrame."""
    return pd.read_csv(bg_filename,
                       header=None,
                       index_col=None,
                       sep='\t',
                       names=bedgraph_schema,
                       dtype=bedgraph_dtypes)


def load_chromsizes(chromsizes_filename):
    return pd.read_csv(chromsizes_filename,
                       header=None,
                       index_col=None,
                       sep='\t',
                       names=('chrom', 'length'))
