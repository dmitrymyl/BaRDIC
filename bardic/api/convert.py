from typing import Dict

import pandas as pd

from .schemas import GeneCoord


def annotation_to_dict(annotation: pd.DataFrame) -> Dict[str, GeneCoord]:
    annotation_dummy: Dict[str, Dict] = annotation[['name', 'chrom', 'start', 'end']].set_index('name').to_dict('index')
    annotation_dict = {rna_name: GeneCoord(chrom=data["chrom"], start=data["start"], end=data["end"])
                       for rna_name, data in annotation_dummy.items()}
    return annotation_dict


def chromsizes_to_dict(chromsizes_series: pd.Series) -> Dict[str, int]:
    return chromsizes_series.to_dict()
