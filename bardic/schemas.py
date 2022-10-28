from collections import namedtuple
from typing import TypedDict

# Schema for Nastya's 6-vote contacts
annot_contact_schema = ["rna_chr",
                        "rna_start",
                        "rna_end",
                        "rna_strand",
                        "rna_cigar",
                        "id",
                        "dna_chr",
                        "dna_start",
                        "dna_end",
                        "sample",
                        "gene_name_un",
                        "gene_type",
                        "gene_strand",
                        "source"]

# Dtypes for Nastya's 6-vote contacts
annot_contact_dtypes = {"rna_chr": "category",
                        "rna_start": "int",
                        "rna_end": "int",
                        "rna_strand": "category",
                        "rna_cigar": "str",
                        "id": "str",
                        "dna_chr": "category",
                        "dna_start": "int",
                        "dna_end": "int",
                        "sample": "category",
                        "gene_name_un": "category",
                        "gene_type": "category",
                        "gene_strand": "category",
                        "source": "category"}

# Schema for Nastya's 3-tiss contacts
unannot_contact_schema = ["id",
                          "rna_chr",
                          "rna_start",
                          "rna_end",
                          "rna_strand",
                          "rna_cigar",
                          "dna_chr",
                          "dna_start",
                          "dna_end",
                          "dna_strand",
                          "dna_cigar",
                          "sample"]

# Dtypes for Nastya's 3-contacts_tiss contacts
unannot_contact_dtypes = {"id": "str",
                          "rna_chr": "category",
                          "rna_start": "int",
                          "rna_end": "int",
                          "rna_strand": "category",
                          "rna_cigar": "str",
                          "dna_chr": "category",
                          "dna_start": "int",
                          "dna_end": "int",
                          "dna_strand": "category",
                          "dna_cigar": "str",
                          "sample": "category"}

# Nastya's schema for genes
gene_schema = ['gene_chr',
               'gene_start',
               'gene_end',
               'gene_strand',
               'gene_type',
               'gene_name',
               'gene_name_un',
               'source']

# Schema for transition between Nastya's and my (RDC) contacts schemas
rdc_transition_contact_schema = ['rna_chr',
                                 'rna_start',
                                 'rna_end',
                                 'dna_chr',
                                 'dna_start',
                                 'dna_end',
                                 'id',
                                 'score',
                                 'rna_strand',
                                 'dna_strand',
                                 'gene_name_un']

# My (RDC) contacts schema
rdc_contact_schema = ['rna_chr',
                      'rna_start',
                      'rna_end',
                      'dna_chr',
                      'dna_start',
                      'dna_end',
                      'id',
                      'score',
                      'rna_strand',
                      'dna_strand',
                      'gene_name']

rdc_contact_schema2 = ['rna_chr',
                       'rna_start',
                       'rna_end',
                       'dna_chr',
                       'dna_start',
                       'dna_end',
                       'gene_name',
                       'score',
                       'rna_strand',
                       'dna_strand']

# My (RDC) contacts dtypes
rdc_contact_dtypes = {'rna_chr': "category",
                      'rna_start': 'int',
                      'rna_end': 'int',
                      'dna_chr': 'category',
                      'dna_start': 'int',
                      'dna_end': 'int',
                      'id': 'str',
                      'score': 'category',
                      'rna_strand': 'category',
                      'dna_strand': 'category',
                      'gene_name': 'category'}

rdc_contact_dtypes2 = {'rna_chr': "category",
                       'rna_start': 'int',
                       'rna_end': 'int',
                       'dna_chr': 'category',
                       'dna_start': 'int',
                       'dna_end': 'int',
                       'gene_name': 'category',
                       'score': 'category',
                       'rna_strand': 'category',
                       'dna_strand': 'category'}

rdc_to_dna = {'dna_chr': 'chrom',
              'dna_start': 'start',
              'dna_end': 'end',
              'gene_name': 'name',
              'score': 'score',
              'dna_strand': 'strand'}

# change of colnames between transition and my schema
transition_to_final = {'gene_name_un': "gene_name"}

# gene table schema for conversion between Nastya's and mine (RDC)
gene_table_schema = ['gene_name',
                     'gene_name_un']

# change of colnames between Nastya's gene schema and mine (RDC)
transition_to_genes = {'gene_name_un': 'gene_name',
                       'source': 'gene_source'}

# My (RDC) gene annotation schema
rdc_gene_schema = ['gene_chr',
                   'gene_start',
                   'gene_end',
                   'gene_name',
                   'gene_score',
                   'gene_strand',
                   'gene_type',
                   'gene_source']

rdc_gene_dtypes = {'gene_chr': 'category',
                   'gene_start': 'int',
                   'gene_end': 'int',
                   'gene_name': 'str',
                   'gene_score': 'str',
                   'gene_strand': 'category',
                   'gene_type': 'category',
                   'gene_source': 'category'}

bed_schema = ['chrom',
              'start',
              'end',
              'name',
              'score',
              'strand']

bed_dtypes_noscore = {'chrom': 'category',
                      'start': 'int',
                      'end': 'int',
                      'name': 'category',
                      'score': 'category',
                      'strand': 'category'}

bedgraph_schema = ['chrom',
                   'start',
                   'end',
                   'score']

bedgraph_dtypes = {'chrom': 'category',
                   'start': 'int',
                   'end': 'int'}

Coords = namedtuple("Coords", "chrom start end")


class GeneCoord(TypedDict):
    chrom: str
    start: int
    end: int
