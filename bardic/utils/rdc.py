from functools import partial
from typing import Any, Dict, Optional

import pandas as pd
from tqdm.contrib.concurrent import process_map

from ..api.binops import make_genomic_track
from ..api.formats import DnaDataset, Rdc
from ..api.mp import adjust_chunksize
from ..api.schemas import RnaAttrs, RnaPixelRecord


def _imputation_value(bg_track: pd.DataFrame, ifactor=0.01) -> float:
    bg_binsizes = bg_track['end'] - bg_track['start']
    mean_bg_count_per_nt = (bg_track['count'] / bg_binsizes).mean()
    ivalue = mean_bg_count_per_nt * ifactor
    return ivalue


def _get_rna_attrs(rna_name: str, attrs_vector: Dict[str, Dict[str, Any]]) -> RnaAttrs:
    rna_attr_names = ('eligible',
                      'cis_factor', 'cis_start', 'trans_bin_size',
                      'total_contacts', 'genic_contacts', 'cis_contacts', 'trans_contacts')
    data = {attr: attrs_vector[attr][rna_name] for attr in rna_attr_names}
    return RnaAttrs(**data)


def _extract_dna_contacts(dna_dataset: DnaDataset, rna_name: str) -> pd.DataFrame:
    return dna_dataset.read_dna_parts(rna_name)


def _cook_pixels(rna_name, dna_contacts, rna_annot, rna_attrs, bg_track, chromdict, ivalue=None):
    dna_track = make_genomic_track(dna_contacts, bg_track, chromdict, rna_annot, rna_attrs, ivalue=ivalue)
    return rna_name, RnaPixelRecord(**{'pixels': dna_track, 'gene_coord': rna_annot, 'rna_attrs': rna_attrs})


def dnadataset_to_rdc(dna_dataset: DnaDataset,
                      bg_track: pd.DataFrame,
                      fname: str,
                      ifactor: Optional[float] = 0.01,
                      n_cores: int = 1,
                      chunksize: int = 50) -> Rdc:
    if not dna_dataset.are_binsizes_selected:
        raise Exception
    if ifactor is None:
        ivalue = None
    else:
        ivalue = _imputation_value(bg_track, ifactor)
    chromdict = dna_dataset.chromsizes
    annotation = dna_dataset.annotation
    rna_attr_names = ('eligible',
                      'cis_factor', 'cis_start', 'trans_bin_size',
                      'total_contacts', 'genic_contacts', 'cis_contacts', 'trans_contacts')
    rnas_attrs_vector = {attr: dna_dataset.read_rna_attribute_batch(attr)
                         for attr in rna_attr_names}
    eligible_rnas = [rna_name
                     for rna_name, eligible in rnas_attrs_vector['eligible'].items()
                     if eligible]

    cook_pixels_prep = partial(_cook_pixels, bg_track=bg_track, chromdict=chromdict, ivalue=ivalue)
    dna_contacts_gen = (_extract_dna_contacts(dna_dataset, rna_name)
                        for rna_name in eligible_rnas)
    gene_coords_gen = (annotation[rna_name]
                       for rna_name in eligible_rnas)
    rna_attrs_gen = (_get_rna_attrs(rna_name, rnas_attrs_vector)
                     for rna_name in eligible_rnas)

    chunksize = adjust_chunksize(len(eligible_rnas), n_cores, chunksize)
    results = process_map(cook_pixels_prep,
                          eligible_rnas,
                          dna_contacts_gen,
                          gene_coords_gen,
                          rna_attrs_gen,
                          max_workers=n_cores,
                          chunksize=chunksize,
                          desc='Creating RDC',
                          unit='RNA')

    rdc = Rdc(fname, chromdict)
    rdc.write_bg_track(bg_track)
    rdc.write_pixels_batch(dict(results))
    return rdc
