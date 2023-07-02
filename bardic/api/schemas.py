from dataclasses import astuple, dataclass
from typing import Dict, Optional, Union, Tuple, Any

import numpy as np
import pandas as pd


bed_schema = ['chrom',
              'start',
              'end',
              'name',
              'score',
              'strand']

bed_dtypes_noscore = {'chrom': 'category',
                      'start': 'int',
                      'end': 'int',
                      'name': 'category',
                      'score': 'category',
                      'strand': 'category'}

bedgraph_schema = ['chrom',
                   'start',
                   'end',
                   'count']

bedgraph_dtypes = {'chrom': 'category',
                   'start': 'int',
                   'end': 'int'}


@dataclass
class GeneCoord:
    chrom: str
    start: int
    end: int


@dataclass
class RnaAttrs:
    eligible: bool
    cis_factor: float
    cis_start: int
    trans_bin_size: int
    total_contacts: int
    genic_contacts: int
    cis_contacts: int
    trans_contacts: int
    spline_params: Optional[Dict] = None


@dataclass
class RnaPixelRecord:
    pixels: pd.DataFrame
    gene_coord: GeneCoord
    rna_attrs: RnaAttrs


@dataclass
class SplineResult:
    t: np.ndarray
    c: np.ndarray
    k: Union[int, float]

    def __iter__(self):
        return iter(astuple(self))


class BedValidator:

    valid_strands: Tuple = ('+', '-', '.')
    fail_reasons: Dict = {'chrom': 'chrom not in chromsizes',
                          'start': 'start < 0',
                          'end': 'end > chromsize',
                          'coord': 'start >= end',
                          'name': 'name contains "/" or starts with "."',
                          'strand': f'strand is not one of {valid_strands}',
                          'chrom_dtype': 'chrom is not str',
                          'start_dtype': 'start is not int',
                          'end_dtype': 'end is not int',
                          'name_dtype': 'name is not str',
                          'score_dtype': 'score is not int or "."',
                          'strand_dtype': 'strand is not str'}

    def __init__(self, chromsizes: Dict):
        self.chromsizes = chromsizes

    def _check_chroms(self, chroms: Union[str, pd.Series]):
        if isinstance(chroms, str):
            return chroms in self.chromsizes
        elif isinstance(chroms, pd.Series):
            return chroms.map(lambda s: s in self.chromsizes)
        else:
            raise Exception

    def _check_starts(self, starts: Union[int, pd.Series]):
        if isinstance(starts, int) or isinstance(starts, pd.Series):
            return starts >= 0
        else:
            raise Exception

    def _check_ends(self, ends: Union[int, pd.Series], chroms: Union[str, pd.Series]):
        if isinstance(ends, int) and isinstance(chroms, str):
            return ends <= self.chromsizes.get(chroms, float('nan'))
        elif isinstance(ends, pd.Series) and isinstance(chroms, pd.Series):
            vec_chromsizes = chroms.map(self.chromsizes)
            return ends <= vec_chromsizes
        else:
            raise Exception

    def _check_coords(self, starts: Union[int, pd.Series], ends: Union[int, pd.Series]):
        if (isinstance(starts, int) and isinstance(ends, int)) or (isinstance(starts, pd.Series) and isinstance(ends, pd.Series)):
            return starts < ends
        else:
            raise Exception

    def _check_names(self, names: Union[str, pd.Series]):
        if isinstance(names, str):
            return not (('/' in names) or names.startswith('.'))
        elif isinstance(names, pd.Series):
            return names.map(lambda s: not (('/' in s) or s.startswith('.')))
        else:
            raise Exception

    def _check_strands(self, strands: Union[str, pd.Series]):
        if isinstance(strands, str):
            return strands in self.valid_strands
        elif isinstance(strands, pd.Series):
            return strands.map(lambda s: s in self.valid_strands)

    def _validate_chroms_dtype(self, chroms: Any):
        if isinstance(chroms, pd.Series):
            return chroms.map(lambda s: isinstance(s, str))
        else:
            return isinstance(chroms, str)

    def _validate_starts_dtype(self, starts: Any):
        if isinstance(starts, pd.Series):
            return starts.map(lambda x: isinstance(x, int))
        else:
            return isinstance(starts, int)

    def _validate_ends_dtype(self, ends: Any):
        if isinstance(ends, pd.Series):
            return ends.map(lambda x: isinstance(x, int))
        else:
            return isinstance(ends, int)

    def _validate_names_dtype(self, names: Any):
        if isinstance(names, pd.Series):
            return names.map(lambda s: isinstance(s, str))
        else:
            return isinstance(names, str)

    def _validate_scores_dtype(self, scores: Any):
        if isinstance(scores, pd.Series):
            return scores.map(lambda x: isinstance(x, int) or x == '.')
        else:
            return isinstance(scores, int) or scores == '.'

    def _validate_strands_dtype(self, strands: Any):
        if isinstance(strands, pd.Series):
            return strands.map(lambda s: isinstance(s, str))
        else:
            return isinstance(strands, str)

    def _filter_df(self, df: pd.DataFrame, selections: Dict) -> Tuple[pd.DataFrame, pd.DataFrame]:
        valid_selection = np.logical_and.reduce(list(selections.values()))
        filtered_df = df[valid_selection].reset_index(drop=True)
        dropped_df = df[~valid_selection].copy()

        fail_reasons_handler = list()
        for fail_category, check_array in selections.items():
            fail_idx = check_array.index[~check_array]
            reason_series = pd.Series([self.fail_reasons[fail_category] + "; "
                                       for _ in range(len(fail_idx))],
                                      index=fail_idx,
                                      dtype='str')
            fail_reasons_handler.append(reason_series)

        fail_reasons_combined = fail_reasons_handler[0].str.cat(fail_reasons_handler[1:], join='outer', na_rep='')\
                                                       .str.rstrip('; ')
        dropped_df['reason'] = fail_reasons_combined
        dropped_df = dropped_df.reset_index(drop=True)

        return filtered_df, dropped_df

    def _validate_df_dtypes(self, df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
        dtype_selections = {'chrom_dtype': self._validate_chroms_dtype(df['chrom']),
                            'start_dtype': self._validate_starts_dtype(df['start']),
                            'end_dtype': self._validate_ends_dtype(df['end']),
                            'name_dtype': self._validate_names_dtype(df['name']),
                            'score_dtype': self._validate_scores_dtype(df['score']),
                            'strand_dtype': self._validate_strands_dtype(df['strand'])}
        return self._filter_df(df, dtype_selections)

    def _validate_df_values(self, df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
        check_selections = {'chrom': self._check_chroms(df['chrom']),
                            'start': self._check_starts(df['start']),
                            'end': self._check_ends(df['end'], df['chrom']),
                            'coord': self._check_coords(df['start'], df['end']),
                            'name': self._check_names(df['name']),
                            'strand': self._check_strands(df['strand'])}

        check_selections['end'] = pd.Series(np.where(check_selections['chrom'], check_selections['end'], True))
        return self._filter_df(df, check_selections)

    def validate_df(self, df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
        dtype_correct_df, dtype_failed_df = self._validate_df_dtypes(df)
        value_correct_df, value_failed_df = self._validate_df_values(dtype_correct_df)
        return value_correct_df, pd.concat([dtype_failed_df, value_failed_df], ignore_index=True)

    def validate_file(self, input_file, filtered_file, dropped_file):
        raise NotImplementedError


class ContactsValidator(BedValidator):

    def __init__(self, chromsizes: Dict, annotation: Dict):
        super().__init__(chromsizes)
        self.annotation = annotation
        self.fail_reasons['name_annot'] = "name is not found in the annotation"

    def _check_name_in_annotation(self, names: Union[str, pd.Series]):
        if isinstance(names, str):
            return names in self.annotation
        elif isinstance(names, pd.Series):
            return names.map(lambda s: s in self.annotation)
        else:
            raise Exception

    def _validate_name_in_annotation(self, df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
        name_annot_selection = {'name_annot': self._check_name_in_annotation(df['name'])}
        return self._filter_df(df, name_annot_selection)

    def validate_df(self, df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
        bed_correct_df, bed_failed_df = super().validate_df(df)
        name_annot_df, name_failed_df = self._validate_name_in_annotation(bed_correct_df)
        return name_annot_df, pd.concat([bed_failed_df, name_failed_df], ignore_index=True)

    def validate_file(self, input_file, filtered_file, dropped_file):
        raise NotImplementedError
