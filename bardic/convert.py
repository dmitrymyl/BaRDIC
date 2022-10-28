from typing import Dict
from .schemas import GeneCoord

import pandas as pd


def annotation_to_dict(annotation_df: pd.DataFrame) -> Dict[str, GeneCoord]:
    annotation_dummy: Dict[str, GeneCoord] = annotation_df[['name', 'chrom', 'start', 'end']].set_index('name').to_dict('index')
    # because mypy doesn't accept any other solution, we have to use kwargs in GeneCoord
    annotation_dict = {rna_name: GeneCoord(chrom=data['chrom'], start=data['start'], end=data['end'])
                       for rna_name, data in annotation_dummy.items()}
    return annotation_dict


def chromsizes_to_dict(chromsizes_series: pd.Series) -> Dict[str, int]:
    return chromsizes_series.to_dict()
