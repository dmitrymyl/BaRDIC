from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
from typing import Any, Dict, Optional

import bioframe as bf
import h5py
import numpy as np
import pandas as pd

from .schemas import GeneCoord, RnaAttrs, RnaPixelRecord


class DnaDataset:
    chrom_groupname: str = "chrom_sizes"
    dna_groupname: str = "dna_parts"

    def __init__(self,
                 fname: str,
                 chromsizes: Optional[Dict[str, int]] = None,
                 annotation: Optional[Dict[str, GeneCoord]] = None) -> None:
        self.fname: Path = Path(fname)
        if not self.fname.exists():
            with h5py.File(self.fname, 'w') as f:
                f.attrs["binsizes_selected"] = False

        self._binsizes_selected = None

        self._chromsizes: Optional[Dict[str, int]] = chromsizes
        if chromsizes is not None:
            self.write_chromsizes()
        if annotation is not None:
            self.validate_annotation(annotation, self.chromsizes)
        self._annotation: Optional[Dict[str, GeneCoord]] = annotation

    @property
    def binsizes_selected(self) -> bool:
        if self._binsizes_selected is None:
            with h5py.File(self.fname, 'r') as f:
                binsizes_selected = f.attrs['binsizes_selected']
                if binsizes_selected in (True, False):
                    self._binsizes_selected = binsizes_selected
                else:
                    raise Exception
        return self._binsizes_selected

    @binsizes_selected.setter
    def binsizes_selected(self, indicator: bool) -> None:
        if indicator not in (True, False):
            raise ValueError
        with h5py.File(self.fname, 'w') as f:
            f.attrs['binsizes_selected'] = indicator
        self._binsizes_selected = indicator

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
        for _rna_name, rna_annot in annotation_dict.items():
            if rna_annot.start >= rna_annot.end:
                raise Exception  # rna_name, start, end
            if rna_annot.chrom not in chromsizes_dict:
                raise Exception  # rna_name, chrom
            if rna_annot.start < 0:
                raise Exception  # rna_name, start
            if rna_annot.end > chromsizes_dict[rna_annot.chrom]:
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
        rna_group.attrs['total_contacts'] = dna_parts.shape[0]
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
            sizes = {rna_name: rna_group.attrs['total_contacts'] for rna_name, rna_group in dna_group.items()}
        return sizes

    def write_attribute(self, name: str, data: Dict[str, Any]) -> None:
        with h5py.File(self.fname, 'a') as f:
            if self.dna_groupname not in f:
                raise Exception
            dna_group = f[self.dna_groupname]
            for rna_name, value in data.items():
                if rna_name in dna_group:
                    dna_group[rna_name].attrs[name] = value

    def read_attribute(self, attrname: str) -> Dict[str, Any]:
        with h5py.File(self.fname, 'r') as f:
            if self.dna_groupname not in f:
                raise Exception
            dna_group = f[self.dna_groupname]
            data = {rna_name: rna_group.attrs.get(attrname) for rna_name, rna_group in dna_group.items()}
        return data


class Rdc:
    pixels_cols = {'start': 'int64',
                   'end': 'int64',
                   'signal_count': 'int64',
                   'signal_prob': 'float',
                   'impute': 'bool',
                   'bg_count': 'float',
                   'bg_prob': 'float',
                   'raw_bg_prob': 'float',
                   'scaling_factor': 'float',
                   'fc': 'float',
                   'pvalue': 'float',
                   'qvalue': 'float'}

    chrom_groupname: str = "chrom_sizes"
    bg_groupname: str = "background"
    pixels_groupname: str = "pixels"

    def __init__(self,
                 fname: str,
                 chromsizes: Optional[Dict[str, int]] = None) -> None:
        self.fname: Path = Path(fname)
        if not self.fname.exists():
            with h5py.File(self.fname, 'w'):
                pass

        self._chromsizes: Optional[Dict[str, int]] = chromsizes
        if chromsizes is not None:
            self.write_chromsizes()

    @property
    def chromsizes(self) -> Dict[str, int]:
        if self._chromsizes is None:
            try:
                self._chromsizes = self.read_chromsizes()
            except Exception:
                raise Exception
        return self._chromsizes

    @chromsizes.setter
    def chromsizes(self, chromdict: Dict[str, int]) -> None:
        self._chromsizes = chromdict
        self.write_chromsizes()

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

    def write_bg_track(self, bg_track: pd.DataFrame) -> None:
        with h5py.File(self.fname, 'a') as f:
            if self.bg_groupname not in f:
                del self.bg_groupname
            bg_group = f.create_group(self.bg_groupname)
            for chrom_name, chrom_df in bg_track.groupby('chrom'):
                chrom_group = bg_group.create_group(chrom_name)
                for key in ('start', 'end', 'count'):
                    chrom_group.create_dataset(key, data=chrom_df[key].values)

    def read_bg_track(self) -> pd.DataFrame:
        with h5py.File(self.fname, 'a') as f:
            if self.bg_groupname not in f:
                raise Exception
            bg_group = f[self.bg_groupname]
            bg_dfs = list()
            for chrom_name, chrom_group in bg_group.items():
                chrom_df = pd.DataFrame({key: chrom_group[key][()]
                                         for key in ('start', 'end', 'count')})
                chrom_df['chrom'] = chrom_name
                chrom_df = chrom_df[['chrom', 'start', 'end', 'count']]
                bg_dfs.append(chrom_df)
            return pd.concat(bg_dfs, ignore_index=True)

    def _write_pixels_single(self,
                             f: h5py.File,
                             rna_name: str,
                             pixels_df: pd.DataFrame,
                             rna_coords: GeneCoord,
                             rna_attrs: RnaAttrs) -> None:
        if self.pixels_groupname not in f:
            pixels_group = f.create_group(self.pixels_groupname)
        else:
            pixels_group = f[self.pixels_groupname]

        if rna_name in pixels_group:
            del pixels_group[rna_name]
        rna_group = pixels_group.create_group(rna_name)

        rna_dict = asdict(rna_attrs)
        rna_dict.update(asdict(rna_coords))
        for key, value in rna_dict.items():
            if value is not None:
                rna_group.attrs[key] = value

        for chrom_name, chrom_df in pixels_df.groupby('chrom'):
            chrom_group = rna_group.create_group(chrom_name)
            for col_name, col_dtype in self.pixels_cols.items():
                if col_name in chrom_df:
                    col_data = chrom_df[col_name].values.astype(col_dtype)
                    chrom_group.create_dataset(col_name, data=col_data)

    def write_pixels_single(self, rna_name, pixels_df, rna_coords, rna_attrs) -> None:
        with h5py.File(self.fname) as f:
            self._write_pixels_single(f, rna_name, pixels_df, rna_coords, rna_attrs)

    def write_pixels(self, data: Dict[str, RnaPixelRecord]) -> None:
        with h5py.File(self.fname, 'a') as f:
            for rna_name, rna_record in data.items():
                pixels_df = rna_record.pixels
                rna_coord = rna_record.gene_coord
                rna_attrs = rna_record.rna_attrs
                self._write_pixels_single(f, rna_name, pixels_df, rna_coord, rna_attrs)

    def write_array(self, name: str, arrays: Dict[str: np.ndarray]) -> None:
        raise NotImplementedError

    def write_attribute(self, name: str, data: Dict[str, Any]) -> None:
        with h5py.File(self.fname, 'a') as f:
            if self.dna_groupname not in f:
                raise Exception
            dna_group = f[self.dna_groupname]
            for rna_name, value in data.items():
                if rna_name in dna_group:
                    dna_group[rna_name].attrs[name] = value

    def read_attribute(self, attrname: str) -> Dict[str, Any]:
        with h5py.File(self.fname, 'r') as f:
            if self.dna_groupname not in f:
                raise Exception
            dna_group = f[self.dna_groupname]
            data = {rna_name: rna_group.attrs.get(attrname) for rna_name, rna_group in dna_group.items()}
        return data
