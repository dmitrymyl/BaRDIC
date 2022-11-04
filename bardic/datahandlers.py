from pathlib import Path
from typing import Dict, Optional

import bioframe as bf
import h5py
import numpy as np
import pandas as pd
from .schemas import GeneCoord


class DnaParts:

    def __init__(self,
                 fname: str,
                 chromsizes: Optional[Dict[str, int]] = None,
                 annotation: Optional[Dict[str, GeneCoord]] = None,
                 chrom_groupname: str = 'chrom_sizes',
                 dna_groupname: str = 'dna_parts') -> None:
        self.fname: Path = Path(fname)
        if not self.fname.exists():
            with h5py.File(self.fname, 'w'):
                pass
        self.chrom_groupname: str = chrom_groupname
        self.dna_groupname: str = dna_groupname

        self._chromsizes: Optional[Dict[str, int]] = chromsizes
        if chromsizes is not None:
            self.write_chromsizes()
        if annotation is not None:
            self.validate_annotation(annotation, self.chromsizes)
        self._annotation: Optional[Dict[str, GeneCoord]] = annotation

    @property
    def chromsizes(self) -> Dict[str, int]:
        if self._chromsizes is None:
            try:
                self._chromsizes = self.read_chromsizes()
            except Exception:
                raise Exception
        return self._chromsizes

    @chromsizes.setter
    def chromsizes(self, chromsizes_dict: Dict[str, int]) -> None:
        self._chromsizes = chromsizes_dict
        self.write_chromsizes()

    @property
    def annotation(self) -> Dict[str, GeneCoord]:
        if self._annotation is None:
            try:
                self._annotation = self._get_annotation()
            except Exception:
                raise Exception
        return self._annotation

    @annotation.setter
    def annotation(self, annotation_dict: Dict[str, GeneCoord]) -> None:
        self.validate_annotation(annotation_dict, self.chromsizes)
        self._annotation = annotation_dict

    @staticmethod
    def validate_annotation(annotation_dict: Dict[str, GeneCoord],
                            chromsizes_dict: Dict[str, int]) -> bool:
        for rna_name, rna_annot in annotation_dict.items():
            if rna_annot['start'] >= rna_annot['end']:
                raise Exception  # rna_name, start, end
            if rna_annot['chrom'] not in chromsizes_dict:
                raise Exception  # rna_name, chrom
            if rna_annot['start'] < 0:
                raise Exception  # rna_name, start
            if rna_annot['end'] > chromsizes_dict[rna_annot['chrom']]:
                raise Exception  # rna_name, end, chromsize
        return True

    @staticmethod
    def validate_dna_frame(dna_frame: pd.DataFrame,
                           annotation_dict: Dict[str, GeneCoord],
                           chromsizes_dict: Dict[str, int]) -> bool:
        if not bf.is_bedframe(dna_frame):
            raise Exception
        for rna_name in dna_frame['name'].unique():
            if rna_name not in annotation_dict:
                raise Exception
        for chrom_name in dna_frame['chrom'].unique():
            if chrom_name not in chromsizes_dict:
                raise Exception
        for chrom_name, chrom_df in dna_frame.groupby('chrom'):
            if (chrom_df['end'] > chromsizes_dict[chrom_name]).sum() > 0:
                raise Exception
            if (chrom_df['start'] < 0).sum() > 0:
                raise Exception
        return True

    def write_chromsizes(self) -> None:
        chromsizes_dict = self.chromsizes
        with h5py.File(self.fname, 'a') as f:
            names = np.array(list(chromsizes_dict.keys()), dtype='O')
            sizes = np.array(list(chromsizes_dict.values()), dtype='int64')
            if self.chrom_groupname in f:
                del f[self.chrom_groupname]
            chromsizes_group = f.create_group(self.chrom_groupname)
            chromsizes_group.create_dataset('chrom', data=names)
            chromsizes_group.create_dataset('size', data=sizes)

    def read_chromsizes(self) -> Dict[str, int]:
        if not self.fname.exists():
            raise Exception

        with h5py.File(self.fname, 'r') as f:
            if self.chrom_groupname not in f:
                raise Exception
            chromsizes_group = f[self.chrom_groupname]
            if 'chrom' not in chromsizes_group:
                raise Exception
            names = chromsizes_group['chrom'].asstr()[()]
            if 'size' not in chromsizes_group:
                raise Exception
            sizes = chromsizes_group['size'][()]

        return dict(zip(names, sizes))

    def _write_dna_parts_single(self,
                                f: h5py.File,
                                rna_name: str,
                                dna_parts: pd.DataFrame) -> None:
        annotation_dict = self.annotation
        rna_annot = annotation_dict[rna_name]

        if self.dna_groupname not in f:
            dna_parts_group = f.create_group(self.dna_groupname)
        else:
            dna_parts_group = f[self.dna_groupname]

        if rna_name in dna_parts_group:
            del dna_parts_group[rna_name]
        rna_group = dna_parts_group.create_group(rna_name)

        for key, value in rna_annot.items():
            rna_group.attrs[key] = value
        rna_group.attrs['ncontacts'] = dna_parts.size
        for chrom_name, chrom_df in dna_parts.groupby('chrom'):
            rna_chrom_group = rna_group.create_group(chrom_name)
            starts = chrom_df['start'].values.astype('int64')
            ends = chrom_df['end'].values.astype('int64')
            rna_chrom_group.create_dataset('start', data=starts)
            rna_chrom_group.create_dataset('end', data=ends)

    def write_dna_parts_single(self,
                               rna_name: str,
                               dna_parts: pd.DataFrame) -> None:
        with h5py.File(self.fname, 'a') as f:
            self._write_dna_parts_single(f, rna_name, dna_parts)

    def write_dna_parts(self, dna_frame: pd.DataFrame) -> None:
        annotation_dict = self.annotation
        self.validate_dna_frame(dna_frame, annotation_dict, self.chromsizes)
        self.validate_annotation(annotation_dict, self.chromsizes)
        with h5py.File(self.fname, 'a') as f:
            for rna_name, dna_parts in dna_frame.groupby('name'):
                self._write_dna_parts_single(f, rna_name, dna_parts)

    def _read_dna_parts_single(self, f: h5py.File, rna_name: str) -> pd.DataFrame:
        if self.dna_groupname not in f:
            raise Exception
        dna_parts_group = f[self.dna_groupname]
        if rna_name not in dna_parts_group:
            raise Exception
        rna_group = dna_parts_group[rna_name]
        dna_parts_list = list()
        for chrom_name, chrom_group in rna_group.items():
            starts = chrom_group['start'][()]
            ends = chrom_group['end'][()]
            dna_parts_list.append(pd.DataFrame({'chrom': chrom_name, 'start': starts, 'end': ends}))
        return pd.concat(dna_parts_list, ignore_index=True)

    def read_dna_parts_single(self, rna_name: str) -> pd.DataFrame:
        with h5py.File(self.fname, 'r') as f:
            dna_parts = self._read_dna_parts_single(f, rna_name)
        return dna_parts

    def _get_coordinates(self, f: h5py.File, rna_name: str) -> GeneCoord:
        if self.dna_groupname not in f:
            raise Exception
        dna_group = f[self.dna_groupname]
        if rna_name not in dna_group:
            raise Exception
        rna_group = dna_group[rna_name]
        return GeneCoord(chrom=rna_group.attrs['chrom'],
                         start=rna_group.attrs['start'],
                         end=rna_group.attrs['end'])

    def get_coordinates(self, rna_name: str) -> GeneCoord:
        with h5py.File(self.fname, 'r') as f:
            annot = self._get_coordinates(f, rna_name)
        return annot

    def _get_annotation(self) -> Dict[str, GeneCoord]:
        with h5py.File(self.fname, 'r') as f:
            if self.dna_groupname not in f:
                raise Exception
            dna_group = f[self.dna_groupname]
            annotation_dict = {rna_name: GeneCoord(chrom=rna_group.attrs['chrom'],
                                                   start=rna_group.attrs['start'],
                                                   end=rna_group.attrs['end'])
                               for rna_name, rna_group in dna_group.items()}
        return annotation_dict

    def get_num_contacts(self) -> Dict[str, int]:
        with h5py.File(self.fname, 'r') as f:
            if self.dna_groupname not in f:
                raise Exception
            dna_group = f[self.dna_groupname]
            sizes = {rna_name: rna_group.attrs['ncontacts'] for rna_name, rna_group in dna_group.items()}
        return sizes


def bed2h5(bed_fname: str,
           h5_fname: str,
           chromsizes: Dict[str, int],
           annotation: Dict[str, GeneCoord]) -> DnaParts:
    dataset = DnaParts(h5_fname, chromsizes, annotation)
    dna_frame = bf.read_table(bed_fname, schema='bed6')
    dataset.write_dna_parts(dna_frame)
    return dataset
