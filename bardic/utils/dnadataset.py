from typing import Dict

import bioframe as bf

from ..api.formats import DnaDataset
from ..api.schemas import GeneCoord


def bed2h5(bed_fname: str,
           h5_fname: str,
           chromsizes: Dict[str, int],
           annotation: Dict[str, GeneCoord]) -> DnaDataset:
    dna_frame = bf.read_table(bed_fname, schema='bed6')
    present_rnas = dna_frame['name'].unique()
    refined_annotation = {rna_name: gene_coord
                          for rna_name, gene_coord in annotation.items()
                          if rna_name in present_rnas}
    dataset = DnaDataset(h5_fname, chromsizes, refined_annotation)
    dataset.write_dna_parts_batch(dna_frame)
    return dataset
