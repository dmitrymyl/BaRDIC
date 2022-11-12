from .datahandlers import DnaParts, Rdc
import pandas as pd
from .schemas import RnaAttrs, RnaPixelRecord
from .binops import make_track
from typing import Dict, Any
from functools import partial
from tqdm.contrib.concurrent import process_map


def imputation_value(bg_track: pd.DataFrame, ifactor=0.01) -> float:
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


def _extract_dna_contacts(dna_parts: DnaParts, rna_name: str) -> pd.DataFrame:
    return dna_parts.read_dna_parts_single(rna_name)


def _cook_pixels(rna_name, dna_contacts, rna_annot, rna_attrs, bg_track, chromdict, impute=False, ivalue=0):
    dna_track = make_track(dna_contacts, bg_track, chromdict, rna_annot, rna_attrs, impute=impute, ivalue=ivalue)
    return rna_name, RnaPixelRecord(**{'pixels': dna_track, 'gene_coord': rna_annot, 'rna_attrs': rna_attrs})


def dnaparts_to_rdc(dna_parts: DnaParts, bg_track: pd.DataFrame, fname: str, impute: bool = False, ifactor: float = 0., max_workers=2) -> Rdc:
    ivalue = imputation_value(bg_track, ifactor)
    chromdict = dna_parts.chromsizes
    annotation = dna_parts.annotation
    rna_attr_names = ('eligible',
                      'cis_factor', 'cis_start', 'trans_bin_size',
                      'total_contacts', 'genic_contacts', 'cis_contacts', 'trans_contacts')
    rnas_attrs_vector = {attr: dna_parts.read_attribute(attr)
                         for attr in rna_attr_names}
    eligible_rnas = [rna_name
                     for rna_name, eligible in rnas_attrs_vector['eligible'].items()
                     if eligible]

    cook_pixels_prep = partial(_cook_pixels, bg_track=bg_track, chromdict=chromdict, impute=impute, ivalue=ivalue)
    dna_contacts_gen = (_extract_dna_contacts(dna_parts, rna_name)
                        for rna_name in eligible_rnas)
    gene_coords_gen = (annotation[rna_name]
                       for rna_name in eligible_rnas)
    rna_attrs_gen = (_get_rna_attrs(rna_name, rnas_attrs_vector)
                     for rna_name in eligible_rnas)

    results = process_map(cook_pixels_prep, eligible_rnas, dna_contacts_gen, gene_coords_gen, rna_attrs_gen, max_workers=max_workers)

    rdc = Rdc(fname, chromdict)
    rdc.write_pixels(dict(results))
    return rdc
