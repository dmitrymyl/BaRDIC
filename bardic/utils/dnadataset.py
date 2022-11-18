from typing import Dict

import bioframe as bf

from ..api.formats import DnaDataset
from ..api.schemas import GeneCoord


def bed2h5(bed_fname: str,
           h5_fname: str,
           chromsizes: Dict[str, int],
           annotation: Dict[str, GeneCoord]) -> DnaDataset:
    dataset = DnaDataset(h5_fname, chromsizes, annotation)
    dna_frame = bf.read_table(bed_fname, schema='bed6')
    dataset.write_dna_parts(dna_frame)
    return dataset
