import pandas as pd

from .schemas import rdc_contact_dtypes, rdc_contact_schema2, rdc_contact_dtypes2


def load_dna_parts_of_rdc(rdc_filename, name='gene_name'):
    """Loads DNA parts of contacts from rdc file and returns it as bed-6 df.
    
    Args:
        rdc_filename (str): rdc filename.
        name (str): a column name used as a "name" column in bed-6 df (default: 'gene_name').
    
    Returns:
        pd.DataFrame: DNA parts of contacts in bed-6 manner.
    """
    rdc_data = pd.read_csv(rdc_filename,
                           header=0,
                           index_col=None,
                           sep='\t',
                           dtype=rdc_contact_dtypes,
                           usecols=['dna_chr', 'dna_start', 'dna_end', name, 'score', 'dna_strand'])
    dna_parts = rdc_data.loc[:, ['dna_chr', 'dna_start', 'dna_end', name, 'score', 'dna_strand']]
    dna_parts.columns = ['chrom', 'start', 'end', 'name', 'score', 'strand']
    return dna_parts


def load_rna_parts_of_rdc(rdc_filename, name='gene_name'):
    """Loads RNA parts of contacts from rdc file and returns it as bed-6 df.
    
    Args:
        rdc_filename (str): rdc filename.
        name (str): a column name used as a "name" column in bed-6 df (default: 'gene_name').
    
    Returns:
        pd.DataFrame: RNA parts of contacts in bed-6 manner.
    """
    rdc_data = pd.read_csv(rdc_filename, header=0, index_col=None, sep='\t', dtype=rdc_contact_dtypes)
    rna_parts = rdc_data.loc[:, ['rna_chr', 'rna_start', 'rna_end', name, 'score', 'rna_strand']]
    rna_parts.columns = ['chrom', 'start', 'end', 'name', 'score', 'strand']
    return rna_parts


def load_rdc(rdc_filename):
    """Loads rdc file with specific dtypes into pd.DataFrame."""
    return pd.read_csv(rdc_filename, header=0, index_col=None, sep='\t', dtype=rdc_contact_dtypes)


def load_dna_parts_of_rdc2(rdc_filename, name='gene_name'):
    """Loads DNA parts of contacts from rdc file and returns it as bed-6 df.
    
    Args:
        rdc_filename (str): rdc filename.
        name (str): a column name used as a "name" column in bed-6 df (default: 'gene_name').
    
    Returns:
        pd.DataFrame: DNA parts of contacts in bed-6 manner.
    """
    rdc_data = pd.read_csv(rdc_filename,
                           header=None,
                           index_col=None,
                           sep='\t',
                           names=rdc_contact_schema2,
                           dtype=rdc_contact_dtypes2,
                           usecols=['dna_chr', 'dna_start', 'dna_end', name, 'score', 'dna_strand'])
    dna_parts = rdc_data.loc[:, ['dna_chr', 'dna_start', 'dna_end', name, 'score', 'dna_strand']]
    dna_parts.columns = ['chrom', 'start', 'end', 'name', 'score', 'strand']
    return dna_parts


def load_rna_parts_of_rdc2(rdc_filename, name='gene_name'):
    """Loads RNA parts of contacts from rdc file and returns it as bed-6 df.
    
    Args:
        rdc_filename (str): rdc filename.
        name (str): a column name used as a "name" column in bed-6 df (default: 'gene_name').
    
    Returns:
        pd.DataFrame: RNA parts of contacts in bed-6 manner.
    """
    rdc_data = pd.read_csv(rdc_filename,
                           header=None,
                           index_col=None,
                           sep='\t',
                           names=rdc_contact_schema2,
                           dtype=rdc_contact_dtypes2,
                           usecols=['rna_chr', 'rna_start', 'rna_end', name, 'score', 'rna_strand'])
    rna_parts = rdc_data.loc[:, ['rna_chr', 'rna_start', 'rna_end', name, 'score', 'rna_strand']]
    rna_parts.columns = ['chrom', 'start', 'end', 'name', 'score', 'strand']
    return rna_parts


def load_rdc2(rdc_filename):
    """Loads rdc file with specific dtypes into pd.DataFrame."""
    return pd.read_csv(rdc_filename,
                       header=None,
                       index_col=None,
                       sep='\t',
                       names=rdc_contact_schema2,
                       dtype=rdc_contact_dtypes2)