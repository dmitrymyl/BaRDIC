from .convert import encode_rna_track, encode_bg_track
import ray


def binnify_rnas(contacts_data, chromsizes, annot_dict, selection_dict, bg_track, impute=False, ifactor=0):
    selected_rnas = selection_dict.keys()
    rna_track_dicts = dict()
    for rna_name in selected_rnas:
        rna_track_dicts[rna_name] = encode_rna_track(contacts_data,
                                                     rna_name,
                                                     chromsizes,
                                                     selection_dict,
                                                     annot_dict, 
                                                     bg_track,
                                                     impute=impute,
                                                     imputation_bg=ifactor)
    return rna_track_dicts
