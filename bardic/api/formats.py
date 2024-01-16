"""
Describes two types of data used in BaRDIC and corresponding file types:
1. DnaDataset and corresponding .dnah5 file type.
2. Rdc and corresponding .rdc file type.
"""


from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import bioframe as bf
import h5py
import numpy as np
import pandas as pd

from .schemas import GeneCoord, RnaAttrs, RnaPixelRecord, SplineResult


class StatusProperty:
    """
    A descriptor class that represents a status property.

    This class allows accessing and setting a boolean status property
    stored in an HDF5 file. The property is lazily loaded from the file
    when accessed for the first time.

    Attributes
    ----------
    private_name : str
        The name of the private attribute where the property value is stored.
    public_name : str
        The name of the public property.

    Methods
    -------
    __set_name__(self, owner, name: str) -> None:
        Sets the private and public names of the property.
    __get__(self, obj, objtype=None) -> bool:
        Retrieves the value of the property.
    __set__(self, obj, value: bool) -> None:
        Sets the value of the property.

    Usage
    -----
    Define a class attribute using the `StatusProperty` descriptor
    to represent a status property.

    Example
    -------
    class MyClass:
        status : StatusProperty
            Represents the status of something.
        is_active : StatusProperty
            Represents the active status of something.
    """

    def __set_name__(self, owner, name: str) -> None:
        """
        Set the name of the descriptor.

        Parameters
        ----------
        owner : type
            The class that owns the descriptor.
        name : str
            The name of the descriptor.

        Returns
        -------
        None
            This method does not return anything.
        """
        self.private_name = '_' + name
        self.public_name = name

    def __get__(self, obj, objtype=None) -> bool:
        """
        Retrieve the value of the attribute from the object.

        If the attribute value is `None`, it retrieves the value from the HDF5 file associated with the object.
        The attribute value is expected to be a boolean, and if it is not, an exception is raised.

        Parameters
        ----------
        obj : object
            The object on which the attribute is accessed.
        objtype : type, optional
            The type of the object.

        Returns
        -------
        value : bool
            The value of the attribute.

        Raises
        ------
        Exception
            If the attribute value is not a boolean.

        Examples
        --------
        >>> obj = MyClass()
        >>> obj.my_attribute
        True
        """
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
        """
        Set the value of the attribute.

        Parameters
        ----------
        obj : object
            The object on which the attribute is being set.
        value : bool
            The value to be set.

        Raises
        ------
        ValueError
            If the value is not a boolean.

        Returns
        -------
        None
        """
        if value not in (True, False):
            raise ValueError
        with h5py.File(obj.fname, 'a') as f:
            f.attrs[self.public_name] = value
        setattr(obj, self.private_name, value)


class VersionProperty:
    """
    A descriptor class that represents a version property.

    This class allows accessing and setting a version property of an object.
    The version property is stored in an HDF5 file and can only have values
    from a predefined list of supported versions.

    Parameters:
    -----------
    supported_versions : Tuple[str]
        A tuple of strings representing the supported versions.

    Attributes:
    -----------
    supported_versions : Tuple[str]
        A tuple of strings representing the supported versions.

    private_name : str
        The private name of the version property.

    public_name : str
        The public name of the version property.
    """

    def __init__(self, supported_versions: Tuple[str]) -> None:
        """
        Initialize the Formats class.

        Parameters
        ----------
        supported_versions : Tuple[str]
            A tuple of supported versions.

        Returns
        -------
        None
        """
        self.supported_versions = supported_versions

    def __set_name__(self, owner, name: str) -> None:
        """
        Set the name of the version property.

        This method is automatically called by the Python interpreter when
        the descriptor is assigned to a class attribute.

        Parameters:
        -----------
        owner : type
            The owner class of the version property.

        name : str
            The name of the version property.
        """
        self.private_name = '_' + name
        self.public_name = name

    def __get__(self, obj, objtype=None) -> str:
        """
        Get the value of the version property.

        This method is automatically called when the version property is accessed.

        Parameters:
        -----------
        obj : object
            The object instance that the version property belongs to.

        objtype : type, optional
            The type of the object instance.

        Returns:
        --------
        value : str
            The value of the version property.

        Raises:
        -------
        Exception
            If the version property is not found in the HDF5 file or if the value
            is not in the list of supported versions.
        """
        value = getattr(obj, self.private_name)
        if value is None:
            with h5py.File(obj.fname, 'r') as f:
                value = f.attrs[self.public_name]
                if value in self.supported_versions:
                    setattr(obj, self.private_name, value)
                else:
                    raise Exception("Unsupported version")
        return value

    def __set__(self, obj, value: str) -> None:
        """
        Set the value of the version property.

        This method is automatically called when the version property is set.

        Parameters:
        -----------
        obj : object
            The object instance that the version property belongs to.

        value : str
            The value to set for the version property.

        Raises:
        -------
        ValueError
            If the value is not in the list of supported versions.
        """
        if value not in self.supported_versions:
            raise ValueError("Unsupported version")
        with h5py.File(obj.fname, 'a') as f:
            f.attrs[self.public_name] = value
        setattr(obj, self.private_name, value)


class DnaDataset:
    """
    Handles DNA parts of contacts that are stored in an HDF5 .dnah5 file.

    Parameters
    ----------
    fname
        Filename of the corresponding .dnah5 file.
        If the file doesn't exist, it will be created.
    chromsizes
        Dictionary in a form of `{chromosome_name: chromosome_size}`,
        that contains chromosome sizes.
    annotation
        Dictionary in a form of `{gene_name: gene_coordinates}`,
        that contains gene coordinates.

    Attributes
    ----------
    chrom_groupname: str
        Name of the chromosome sizes group in the .dnah5 file.
    dna_groupname: str
        Name of the dna parts group in the .dnah5 file.
    are_binsizes_selected: bool
        Indicates whether bin sizes were selected.
    version: str
        Version of .dnah5 schema.
    """
    supported_versions = ("1", )
    chrom_groupname: str = "chrom_sizes"
    dna_groupname: str = "dna_parts"

    are_binsizes_selected = StatusProperty()
    version = VersionProperty(supported_versions)

    def __init__(self,
                 fname: str,
                 chromsizes: Optional[Dict[str, int]] = None,
                 annotation: Optional[Dict[str, GeneCoord]] = None) -> None:
        """Constructor

        Parameters
        ----------
        fname
            Filename of the corresponding .dnah5 file.
            If the file doesn't exist, it will be created.
        chromsizes
            Dictionary in a form of `{chromosome_name: chromosome_size}`,
            that contains chromosome sizes.
        annotation
            Dictionary in a form of `{gene_name: gene_coordinates}`,
            that contains gene coordinates.

        Raises
        ------
        ValueError
            In case the version of the existing .dnah5 file is not `"1"`.

        """
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
            self._write_chromsizes()
        if annotation is not None:
            self.validate_annotation(annotation, self.chromsizes)
        self._annotation: Optional[Dict[str, GeneCoord]] = annotation

    @property
    def chromsizes(self) -> Dict[str, int]:
        """Gets chromosome sizes from the .dnah5 file.

        Returns
        -------
        dict
            A dictionary of chromosome sizes in a form
            of `{chromosome_name: chromosome_size}`.

        Raises
        ------
        Exception
            In case there are no chromosome sizes in the file.

        """
        if self._chromsizes is None:
            try:
                self._chromsizes = self._read_chromsizes()
            except Exception:
                raise Exception
        return self._chromsizes

    @chromsizes.setter
    def chromsizes(self, chromsizes_dict: Dict[str, int]) -> None:
        """Writes new chromosome sizes to the .dnah5 file.

        Parameters
        ----------
        chromsizes_dict: dict
            A dictionary of chromosome sizes in a form
            of `{chromosome_name: chromosome_size}`.

        """
        self._chromsizes = chromsizes_dict
        self._write_chromsizes()

    @property
    def annotation(self) -> Dict[str, GeneCoord]:
        """Gets gene annotation from the .dnah5 file.

        Returns
        -------
        dict
            Dictionary in a form of `{gene_name: gene_coordinates}`,
            that contains gene coordinates.

        Raises
        ------
        Exception
            In case there is no gene annotation in the file.

        """
        if self._annotation is None:
            try:
                self._annotation = self._get_annotation()
            except Exception:
                raise Exception
        return self._annotation

    @annotation.setter
    def annotation(self, annotation_dict: Dict[str, GeneCoord]) -> None:
        """Writes new gene annotation in the .dnah5 file.

        Parameters
        ----------
        annotation_dict: dict
            Dictionary in a form of `{gene_name: gene_coordinates}`,
            that contains gene coordinates.
        """
        self.validate_annotation(annotation_dict, self.chromsizes)
        self._annotation = annotation_dict

    @staticmethod
    def validate_annotation(annotation_dict: Dict[str, GeneCoord],
                            chromsizes_dict: Dict[str, int]) -> bool:
        """Validates provided gene annotation.

        Checks that for each RNA:
        1. start does not exceed end.
        2. chromosome name is present in chromosome sizes.
        3. start is nonnegative.
        4. end does not exceed the chromosome size.

        Parameters
        ----------
        annotation_dict: dict
            Dictionary in a form of `{gene_name: gene_coordinates}`,
            that contains gene coordinates.
        chromsizes_dict: dict
            A dictionary of chromosome sizes in a form
            of `{chromosome_name: chromosome_size}`.

        Returns
        -------
        bool
            True if all validatation checks are successful.
            Otherwise, raises exceptions.

        Raises
        ------
        Exception
            In case one of the validation rules is broken
            for any RNA.

        """
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
        """"Validates a dataframe with DNA parts of contacts.

        Validation checks:
        1. The dataframe is a valid BED dataframe (according to bioframe).
        2. For every RNA:
            2.1. RNA name is present in the annotation dictionary.
            2.2. Every chromosome name is present in the chromosome
                 sizes dictionary.
        3. For every DNA part:
            2.1. End does not exceed the chromosome size.
            2.2. Start is nonnegative.

        Returns
        -------
        bool
            True if all validatation checks are successful.
            Otherwise, raises exceptions.

        Raises
        ------
        Exception
            In case one of the validation rules is broken.
        """
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

    def _read_chromsizes(self) -> Dict[str, int]:
        """
        """
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

    def _write_chromsizes(self) -> None:
        chromsizes_dict = self.chromsizes
        with h5py.File(self.fname, 'a') as f:
            names = np.array(list(chromsizes_dict.keys()), dtype='O')
            sizes = np.array(list(chromsizes_dict.values()), dtype='int64')
            if self.chrom_groupname in f:
                del f[self.chrom_groupname]
            chromsizes_group = f.create_group(self.chrom_groupname)
            chromsizes_group.create_dataset('chrom', data=names)
            chromsizes_group.create_dataset('size', data=sizes)

    def _read_dna_parts(self, f: h5py.File, rna_name: str) -> pd.DataFrame:
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

    def read_dna_parts(self, rna_name: str) -> pd.DataFrame:
        with h5py.File(self.fname, 'r') as f:
            dna_parts = self._read_dna_parts(f, rna_name)
        return dna_parts

    def _write_dna_parts(self,
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

    def write_dna_parts(self,
                        rna_name: str,
                        dna_parts: pd.DataFrame) -> None:
        with h5py.File(self.fname, 'a') as f:
            self._write_dna_parts_single(f, rna_name, dna_parts)

    def write_dna_parts_batch(self,
                              dna_frame: pd.DataFrame,
                              rna_col: str = 'name') -> None:
        annotation_dict = self.annotation
        self.validate_dna_frame(dna_frame, annotation_dict, self.chromsizes)
        self.validate_annotation(annotation_dict, self.chromsizes)
        with h5py.File(self.fname, 'a') as f:
            for rna_name, dna_parts in dna_frame.groupby(rna_col):
                self._write_dna_parts(f, rna_name, dna_parts)

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

    def read_rna_attribute_batch(self, attrname: str) -> Dict[str, Any]:
        with h5py.File(self.fname, 'r') as f:
            if self.dna_groupname not in f:
                raise Exception
            dna_group = f[self.dna_groupname]
            data = {rna_name: rna_group.attrs.get(attrname) for rna_name, rna_group in dna_group.items()}
        return data

    def write_rna_attribute_batch(self, attrname: str, data: Dict[str, Any]) -> None:
        with h5py.File(self.fname, 'a') as f:
            if self.dna_groupname not in f:
                raise Exception
            dna_group = f[self.dna_groupname]
            for rna_name, value in data.items():
                if rna_name in dna_group:
                    dna_group[rna_name].attrs[attrname] = value


class Rdc:
    """
    # add docstring in numpydoc style
    """
    supported_versions = ("1", "1.1")  # version 1 is deprecated and will be removed in the future.

    pixels_cols_by_version = {'1': {'start': 'int64',
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
                                    'qvalue': 'float'},
                              '1.1': {'start': 'int64',
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
                                      'qvalue_global': 'float',
                                      'qvalue_rna': 'float'}}
    # pixels_cols = {'start': 'int64',
    #                'end': 'int64',
    #                'signal_count': 'int64',
    #                'signal_prob': 'float',
    #                'impute': 'bool',
    #                'bg_count': 'float',
    #                'bg_prob': 'float',
    #                'raw_bg_prob': 'float',
    #                'scaling_factor': 'float',
    #                'fc': 'float',
    #                'pvalue': 'float',
    #                'qvalue': 'float'}

    chrom_groupname: str = "chrom_sizes"
    bg_groupname: str = "background"
    pixels_groupname: str = "pixels"

    is_scaling_fitted = StatusProperty()
    are_peaks_estimated = StatusProperty()
    version = VersionProperty(supported_versions)

    def __init__(self,
                 fname: str,
                 chromsizes: Optional[Dict[str, int]] = None) -> None:
        self.fname: Path = Path(fname)
        if not self.fname.exists():
            with h5py.File(self.fname, 'w'):
                pass
            self.is_scaling_fitted = False
            self.are_peaks_estimated = False
            self.version = "1.1"  # we write the newest version

        self._chromsizes: Optional[Dict[str, int]] = chromsizes
        if chromsizes is not None:
            self._write_chromsizes()

        self._annotation: Optional[Dict[str, GeneCoord]] = None
        self._is_scaling_fitted: Optional[bool] = None
        self._are_peaks_estimated: Optional[bool] = None

        self._version: Optional[str] = None
        if self.version not in self.supported_versions:
            raise Exception(f"The RDC file version {self.version} is not supported.")

        self.pixels_cols = self.pixels_cols_by_version[self.version]

    @property
    def chromsizes(self) -> Dict[str, int]:
        if self._chromsizes is None:
            try:
                self._chromsizes = self._read_chromsizes()
            except Exception:
                raise Exception
        return self._chromsizes

    @chromsizes.setter
    def chromsizes(self, chromdict: Dict[str, int]) -> None:
        self._chromsizes = chromdict
        self._write_chromsizes()

    def _read_chromsizes(self) -> Dict[str, int]:
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

    def _write_chromsizes(self) -> None:
        chromsizes_dict = self.chromsizes
        with h5py.File(self.fname, 'a') as f:
            names = np.array(list(chromsizes_dict.keys()), dtype='O')
            sizes = np.array(list(chromsizes_dict.values()), dtype='int64')
            if self.chrom_groupname in f:
                del f[self.chrom_groupname]
            chromsizes_group = f.create_group(self.chrom_groupname)
            chromsizes_group.create_dataset('chrom', data=names)
            chromsizes_group.create_dataset('size', data=sizes)

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

    def write_bg_track(self, bg_track: pd.DataFrame) -> None:
        with h5py.File(self.fname, 'a') as f:
            if self.bg_groupname in f:
                del f[self.bg_groupname]
            bg_group = f.create_group(self.bg_groupname)
            for chrom_name, chrom_df in bg_track.groupby('chrom'):
                chrom_group = bg_group.create_group(chrom_name)
                for key in ('start', 'end', 'count'):
                    chrom_group.create_dataset(key, data=chrom_df[key].values)

    def _read_pixels(self,
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

    def read_pixels(self,
                    rna_name: str,
                    value_fields: Optional[List] = None,
                    chrom_type: Optional[str] = None) -> pd.DataFrame:
        with h5py.File(self.fname, 'r') as f:
            return self._read_pixels(f, rna_name, value_fields=value_fields, chrom_type=chrom_type)

    def _write_pixels(self,
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

    def write_pixels(self, rna_name, pixels_df, rna_coords, rna_attrs) -> None:
        with h5py.File(self.fname, 'a') as f:
            self._write_pixels(f, rna_name, pixels_df, rna_coords, rna_attrs)

    def write_pixels_batch(self, data: Dict[str, RnaPixelRecord]) -> None:
        with h5py.File(self.fname, 'a') as f:
            for rna_name, rna_record in data.items():
                pixels_df = rna_record.pixels
                rna_coord = rna_record.gene_coord
                rna_attrs = rna_record.rna_attrs
                self._write_pixels(f, rna_name, pixels_df, rna_coord, rna_attrs)

    def _write_pixels_column(self, f: h5py.File, rna_name: str, col_name: str, col_frame: pd.DataFrame) -> None:
        if self.pixels_groupname not in f:
            raise Exception
        pixels_group = f[self.pixels_groupname]
        if rna_name not in pixels_group:
            raise Exception
        rna_group = pixels_group[rna_name]

        for chrom_name, chrom_df in col_frame.groupby('chrom'):
            if chrom_name not in rna_group:
                raise Exception
            chrom_group = rna_group[chrom_name]
            if col_name not in self.pixels_cols:
                raise Exception
            arr_dtype = self.pixels_cols[col_name]
            arr_data = chrom_df[col_name].values.astype(arr_dtype)
            if col_name in chrom_group:
                chrom_group[col_name][...] = arr_data
            else:
                chrom_group.create_dataset(col_name, data=arr_data)

    def write_pixels_column(self, rna_name: str, col_name: str, col_frame: pd.DataFrame) -> None:
        with h5py.File(self.fname, 'a') as f:
            self._write_pixels_column(f, rna_name, col_name, col_frame)

    def write_pixels_column_batch(self, col_name: str, col_frames: Dict[str, pd.DataFrame]) -> None:
        with h5py.File(self.fname, 'a') as f:
            for rna_name, col_frame in col_frames.items():
                self._write_pixels_column(f, rna_name, col_name, col_frame)

    def read_rna_attribute_batch(self, attrname) -> Dict[str, Any]:
        with h5py.File(self.fname, 'r') as f:
            if self.pixels_groupname not in f:
                raise Exception
            pixels_group = f[self.pixels_groupname]
            data = {rna_name: rna_group.attrs.get(attrname)
                    for rna_name, rna_group in pixels_group.items()}
        return data

    def write_rna_attribute_batch(self, attrname: str, data: Dict[str, Any]) -> None:
        with h5py.File(self.fname, 'a') as f:
            if self.pixels_groupname not in f:
                raise Exception
            pixels_group = f[self.pixels_groupname]
            for rna_name, value in data.items():
                if rna_name in pixels_group:
                    pixels_group[rna_name].attrs[attrname] = value

    def read_scaling(self, rna_name: str) -> SplineResult:
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

    def read_scaling_batch(self) -> Dict[str, SplineResult]:
        rna_spline_t = self.read_rna_attribute_batch('scaling_spline_t')
        rna_spline_c = self.read_rna_attribute_batch('scaling_spline_c')
        rna_spline_k = self.read_rna_attribute_batch('scaling_spline_k')
        scaling_splines = {rna_name: SplineResult(t=rna_spline_t[rna_name],
                                                  c=rna_spline_c[rna_name],
                                                  k=rna_spline_k[rna_name])
                           for rna_name in rna_spline_t}
        return scaling_splines

    def write_scaling_batch(self, scaling_splines: Dict[str, SplineResult]) -> None:
        rna_spline_t = {rna_name: tck.t for rna_name, tck in scaling_splines.items()}
        rna_spline_c = {rna_name: tck.c for rna_name, tck in scaling_splines.items()}
        rna_spline_k = {rna_name: tck.k for rna_name, tck in scaling_splines.items()}
        self.write_rna_attribute_batch('scaling_spline_t', rna_spline_t)
        self.write_rna_attribute_batch('scaling_spline_c', rna_spline_c)
        self.write_rna_attribute_batch('scaling_spline_k', rna_spline_k)

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
    def annotation(self) -> Dict[str, GeneCoord]:
        """
        Returns a dictionary of gene coordinates for the current sequence.

        If the annotation has not been loaded yet, it will be loaded from the appropriate file.

        Returns
        -------
        Dict[str, GeneCoord]
            A dictionary of gene coordinates for the current sequence.
        """
        if self._annotation is None:
            self._annotation = self._get_annotation()
        return self._annotation
