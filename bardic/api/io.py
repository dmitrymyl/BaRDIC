import bioframe as bf
import pandas as pd

from .schemas import bedgraph_dtypes, bedgraph_schema


def load_bedgraph(bg_filename: str) -> pd.DataFrame:
    """Loads bedgraph file with specific dtypes into pd.DataFrame."""
    return pd.read_csv(bg_filename,
                       header=None,
                       index_col=None,
                       sep='\t',
                       names=bedgraph_schema,
                       dtype=bedgraph_dtypes)


def load_annotation(filename, schema='bed6'):
    return bf.read_table(filename, schema=schema)


def load_chromsizes(chromsizes_filename: str) -> pd.Series:
    return pd.read_csv(chromsizes_filename,
                       header=None,
                       index_col=None,
                       sep='\t',
                       names=('chrom', 'length')).set_index('chrom').squeeze()


def fetch_chromsizes(genome: str) -> pd.Series:
    return bf.fetch_chromsizes(genome)
