from dataclasses import dataclass
from typing import Dict, Optional, Union

import numpy as np
import pandas as pd


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
                   'count']

bedgraph_dtypes = {'chrom': 'category',
                   'start': 'int',
                   'end': 'int'}


@dataclass
class GeneCoord:
    chrom: str
    start: int
    end: int


@dataclass
class RnaAttrs:
    eligible: bool
    cis_factor: float
    cis_start: int
    trans_bin_size: int
    total_contacts: int
    genic_contacts: int
    cis_contacts: int
    trans_contacts: int
    spline_params: Optional[Dict] = None


@dataclass
class RnaPixelRecord:
    pixels: pd.DataFrame
    gene_coord: GeneCoord
    rna_attrs: RnaAttrs


@dataclass
class SplineResult:
    t: np.ndarray
    c: np.ndarray
    k: Union[int, float]
