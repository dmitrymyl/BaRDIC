from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
from typing import Any, Dict, List, Optional

import bioframe as bf
import h5py
import numpy as np
import pandas as pd

from .schemas import GeneCoord, RnaAttrs, RnaPixelRecord, SplineResult


class StatusProperty:

    def __set_name__(self, owner, name: str) -> None:
        self.private_name = '_' + name
        self.public_name = name

    def __get__(self, obj, objtype=None) -> bool:
        value = getattr(obj, self.private_name)
        if value is None:
            with h5py.File(obj.fname, 'r') as f:
                value = f.attrs[self.public_name]
                if value in (True, False):
                    setattr(obj, self.private_name, value)
                else:
                    raise Exception
        return value

    def __set__(self, obj, value: bool) -> None:
        if value not in (True, False):
            raise ValueError
        with h5py.File(obj.fname, 'a') as f:
            f.attrs[self.public_name] = value
        setattr(obj, self.private_name, value)


class VersionProperty:

    def __set_name__(self, owner, name: str) -> None:
        self.private_name = '_' + name
        self.public_name = name

    def __get__(self, obj, objtype=None) -> bool:
        value = getattr(obj, self.private_name)
        if value is None:
            with h5py.File(obj.fname, 'r') as f:
                value = f.attrs[self.public_name]
                if value in ("1", ):
                    setattr(obj, self.private_name, value)
                else:
                    raise Exception
        return value

    def __set__(self, obj, value: bool) -> None:
        if value not in ("1", ):
            raise ValueError
        with h5py.File(obj.fname, 'a') as f:
            f.attrs[self.public_name] = value
        setattr(obj, self.private_name, value)


class DnaDataset:
    chrom_groupname: str = "chrom_sizes"
    dna_groupname: str = "dna_parts"

    are_binsizes_selected = StatusProperty()
    version = VersionProperty()

    def __init__(self,
                 fname: str,
                 chromsizes: Optional[Dict[str, int]] = None,
                 annotation: Optional[Dict[str, GeneCoord]] = None) -> None:
        self.fname: Path = Path(fname)
        if not self.fname.exists():
            with h5py.File(self.fname, 'w'):
                pass
            self.are_binsizes_selected = False
            self.version = "1"

        self._are_binsizes_selected: Optional[bool] = None
        self._version: Optional[str] = None
        if self.version != "1":
            raise ValueError

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

        for key, value in asdict(rna_annot).items():
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

    is_scaling_fitted = StatusProperty()
    are_peaks_estimated = StatusProperty()
    version = VersionProperty()

    def __init__(self,
                 fname: str,
                 chromsizes: Optional[Dict[str, int]] = None) -> None:
        self.fname: Path = Path(fname)
        if not self.fname.exists():
            with h5py.File(self.fname, 'w'):
                pass
            self.is_scaling_fitted = False
            self.are_peaks_estimated = False
            self.version = "1"

        self._chromsizes: Optional[Dict[str, int]] = chromsizes
        if chromsizes is not None:
            self.write_chromsizes()

        self._annotation: Optional[Dict[str, GeneCoord]] = None
        self._is_scaling_fitted: Optional[bool] = None
        self._are_peaks_estimated: Optional[bool] = None
        self._version: Optional[str] = None
        if self.version != "1":
            raise Exception

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
            if self.bg_groupname in f:
                del f[self.bg_groupname]
            bg_group = f.create_group(self.bg_groupname)
            for chrom_name, chrom_df in bg_track.groupby('chrom'):
                chrom_group = bg_group.create_group(chrom_name)
                for key in ('start', 'end', 'count'):
                    chrom_group.create_dataset(key, data=chrom_df[key].values)

    def read_bg_track(self) -> pd.DataFrame:
        with h5py.File(self.fname, 'r') as f:
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
        with h5py.File(self.fname, 'a') as f:
            self._write_pixels_single(f, rna_name, pixels_df, rna_coords, rna_attrs)

    def write_pixels(self, data: Dict[str, RnaPixelRecord]) -> None:
        with h5py.File(self.fname, 'a') as f:
            for rna_name, rna_record in data.items():
                pixels_df = rna_record.pixels
                rna_coord = rna_record.gene_coord
                rna_attrs = rna_record.rna_attrs
                self._write_pixels_single(f, rna_name, pixels_df, rna_coord, rna_attrs)

    def _write_array_single(self, f: h5py.File, rna_name: str, arr_name: str, array: pd.DataFrame) -> None:
        if self.pixels_groupname not in f:
            raise Exception
        pixels_group = f[self.pixels_groupname]
        if rna_name not in pixels_group:
            raise Exception
        rna_group = pixels_group[rna_name]

        for chrom_name, chrom_df in array.groupby('chrom'):
            if chrom_name not in rna_group:
                raise Exception
            chrom_group = rna_group[chrom_name]
            if arr_name not in self.pixels_cols:
                raise Exception
            arr_dtype = self.pixels_cols[arr_name]
            arr_data = chrom_df[arr_name].values.astype(arr_dtype)
            if arr_name in chrom_group:
                chrom_group[arr_name][...] = arr_data
            else:
                chrom_group.create_dataset(arr_name, data=arr_data)

    def write_array_single(self, rna_name: str, arr_name: str, array: pd.DataFrame) -> None:
        with h5py.File(self.fname, 'a') as f:
            self._write_array_single(f, rna_name, arr_name, array)

    def write_array(self, arr_name: str, arrays: Dict[str, pd.DataFrame]) -> None:
        with h5py.File(self.fname, 'a') as f:
            for rna_name, array in arrays.items():
                self._write_array_single(f, rna_name, arr_name, array)

    def write_attribute(self, attrname: str, data: Dict[str, Any]) -> None:
        with h5py.File(self.fname, 'a') as f:
            if self.pixels_groupname not in f:
                raise Exception
            pixels_group = f[self.pixels_groupname]
            for rna_name, value in data.items():
                if rna_name in pixels_group:
                    pixels_group[rna_name].attrs[attrname] = value

    def read_attribute(self, attrname) -> Dict[str, Any]:
        with h5py.File(self.fname, 'r') as f:
            if self.pixels_groupname not in f:
                raise Exception
            pixels_group = f[self.pixels_groupname]
            data = {rna_name: rna_group.attrs.get(attrname) for rna_name, rna_group in pixels_group.items()}
        return data

    def write_scaling_splines(self, scaling_splines: Dict[str, SplineResult]) -> None:
        rna_spline_t = {rna_name: tck.t for rna_name, tck in scaling_splines.items()}
        rna_spline_c = {rna_name: tck.c for rna_name, tck in scaling_splines.items()}
        rna_spline_k = {rna_name: tck.k for rna_name, tck in scaling_splines.items()}
        self.write_attribute('scaling_spline_t', rna_spline_t)
        self.write_attribute('scaling_spline_c', rna_spline_c)
        self.write_attribute('scaling_spline_k', rna_spline_k)

    def read_scaling_splines(self) -> Dict[str, SplineResult]:
        rna_spline_t = self.read_attribute('scaling_spline_t')
        rna_spline_c = self.read_attribute('scaling_spline_c')
        rna_spline_k = self.read_attribute('scaling_spline_k')
        scaling_splines = {rna_name: SplineResult(t=rna_spline_t[rna_name], c=rna_spline_c[rna_name], k=rna_spline_k[rna_name])
                           for rna_name in rna_spline_t}
        return scaling_splines

    def read_scaling_single(self, rna_name: str) -> SplineResult:
        if not self.is_scaling_fitted:
            raise Exception
        with h5py.File(self.fname, 'r') as f:
            if self.pixels_groupname not in f:
                raise Exception
            pixels_group = f[self.pixels_groupname]
            if rna_name not in pixels_group:
                raise Exception
            rna_group = pixels_group[rna_name]
            spline_result = SplineResult(t=rna_group.attrs['scaling_spline_t'],
                                         c=rna_group.attrs['scaling_spline_c'],
                                         k=rna_group.attrs['scaling_spline_k'])
        return spline_result

    def _read_pixels_single(self,
                            f: h5py.File,
                            rna_name: str,
                            value_fields: Optional[List] = None,
                            chrom_type: Optional[str] = None) -> pd.DataFrame:
        if self.pixels_groupname not in f:
            raise Exception
        pixels_group = f[self.pixels_groupname]

        if rna_name not in pixels_group:
            raise Exception
        rna_group = pixels_group[rna_name]

        if chrom_type is None:
            valid_chroms = list(rna_group.keys())
        elif chrom_type == 'cis':
            valid_chroms = [rna_group.attrs['chrom']]
        elif chrom_type == 'trans':
            valid_chroms = [chrom for chrom in rna_group.keys() if chrom != rna_group.attrs['chrom']]
        else:
            raise ValueError

        mandatory_fields = ['start', 'end']
        all_fields = mandatory_fields
        chrom_dfs = list()

        for chrom_name, chrom_group in rna_group.items():
            if chrom_name not in valid_chroms:
                continue

            if value_fields is None:
                value_fields = [item for item in chrom_group.keys() if item not in mandatory_fields]
            all_fields = mandatory_fields + value_fields

            data = {field_name: field_data[()]
                    for field_name, field_data in chrom_group.items()
                    if field_name in all_fields}
            data['chrom'] = chrom_name
            chrom_dfs.append(pd.DataFrame(data)[['chrom'] + all_fields])

        result = pd.concat(chrom_dfs, ignore_index=True)
        return result

    def read_pixels_single(self,
                           rna_name: str,
                           value_fields: Optional[List] = None,
                           chrom_type: Optional[str] = None) -> pd.DataFrame:
        with h5py.File(self.fname, 'r') as f:
            return self._read_pixels_single(f, rna_name, value_fields=value_fields, chrom_type=chrom_type)

    def _get_annotation(self) -> Dict[str, GeneCoord]:
        with h5py.File(self.fname, 'r') as f:
            if self.pixels_groupname not in f:
                raise Exception
            pixels_group = f[self.pixels_groupname]
            annotation_dict = {rna_name: GeneCoord(chrom=rna_group.attrs['chrom'],
                                                   start=rna_group.attrs['start'],
                                                   end=rna_group.attrs['end'])
                               for rna_name, rna_group in pixels_group.items()}
        return annotation_dict

    @property
    def annotation(self):
        if self._annotation is None:
            self._annotation = self._get_annotation()
        return self._annotation
