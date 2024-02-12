from pathlib import Path
from typing import Union, Tuple

import bioframe as bf
import numpy as np
import pandas as pd
import scipy.stats as ss
from numpy.typing import ArrayLike


Numeric = Union[int, float]


def _point_to_bed(point: int,
                  chrom: str,
                  name: str = '.',
                  score: Union[str, int] = '.',
                  strand: str = '.') -> str:
    return "\t".join([chrom, str(point), str(point + 1), name, str(score), strand])


def _get_B_C(D_high: Numeric, D_low: Numeric, F_high: float, F_low: float) -> Tuple[float, float]:
    B = (D_low + D_high) / (D_low - D_high) * np.log(F_low / F_high)
    C = 2 / (D_low - D_high) * np.log(F_high / F_low)
    return B, C


def _get_simulator_from_bins(bins: pd.DataFrame) -> ss.rv_histogram:
    bin_edges = np.concatenate([bins['start'].values, [bins['end'].iloc[-1]]])
    bin_freqs = bins['count'] / bins['count'].sum()
    return ss.rv_histogram((bin_freqs, bin_edges), density=True)


def _sample_peak_summits(simulator: ss.rv_histogram, peak_exp: Numeric, size: int) -> Tuple[ArrayLike, ArrayLike]:
    summits = simulator.rvs(size=size).round().astype('int')
    abundance_sampler = ss.poisson(peak_exp)
    peak_abundances = list()
    while len(peak_abundances) < len(summits):
        abundance = round(abundance_sampler.rvs())
        if abundance > 0:
            peak_abundances.append(abundance)
    peak_abundances = np.array(peak_abundances)
    return summits, peak_abundances


def _sample_specific_contacts(peak_summit: int, abundance: int, sd: int, L_chrom: int) -> ArrayLike:
    n_specif_contacts = abundance - 1  # -1 to account for the peak summit
    peak_sampler = ss.norm(loc=peak_summit, scale=sd)
    specif_coords = peak_sampler.rvs(size=n_specif_contacts).round().astype('int')  # if n_specif_contacts == 0, it will return an empty array
    specif_coords = np.append(specif_coords, peak_summit)
    specif_coords = np.where(specif_coords < 0, 0, specif_coords)
    specif_coords = np.where(specif_coords > L_chrom, L_chrom - 1, specif_coords)  # -1 to account for the edge and +1 point
    return specif_coords


def _sigmoid(d: ArrayLike, A: Numeric, B: float, C: float) -> ArrayLike:
    return A / (1 + np.e ** (B + C * d))


def simulate(outdir: Union[str, Path],
             rna_name: str = 'simulRNA',
             L_cis: int = 100_000_000,
             L_genome: int = 2_000_000_000,
             bin_size: int = 1_000,
             bg_exp: Numeric = 8,
             imputation_factor: float = 0.01,
             A: Numeric = 3,
             D_high: Numeric = 3,
             D_low: Numeric = 7,
             F_high: float = 0.9,
             F_low: float = 0.1,
             B_left: float = 0.05,
             B_right: float = 0.95,
             N_total: int = 20_000,
             frac_cis: float = 0.75,
             N_peaks: int = 5,
             frac_cis_peaks: float = 0.2,
             peak_exp: int = 100,
             sd: int = 5_000,
             ):

    outdir = Path(outdir)
    if not outdir.exists():
        outdir.mkdir()

    # Get chromsizes

    L_trans = L_genome - L_cis
    chromsizes = pd.Series([L_cis, L_trans], index=['cis', 'trans'])

    chromsizes.to_csv(outdir / "chrom.sizes", index=True, header=False, sep='\t')

    # Simulate bins and the background

    simulated_bins = bf.binnify(chromsizes, bin_size)
    bg_simulator = ss.poisson(bg_exp)
    bg_counts = bg_simulator.rvs(simulated_bins.shape[0], random_state=None)
    simulated_bins['count'] = bg_counts

    simulated_bins.to_csv(outdir / 'background.bedGraph', header=False, index=False, sep='\t')

    # Impute bins with zero counts

    simulated_bins.loc[simulated_bins['count'] == 0, 'count'] = imputation_factor * bg_exp

    # Drop the source gene point

    random_generator = np.random.default_rng()

    gene_point = random_generator.integers(L_cis * B_left, L_cis * B_right)
    gene_record = _point_to_bed(gene_point, 'cis', rna_name)

    with open(outdir / 'annot.bed', 'w') as outfile:
        outfile.write(gene_record + '\n')

    # Calculate separations between the gene and the bins

    cis_bins = simulated_bins.query('chrom == "cis"').copy()
    trans_bins = simulated_bins.query('chrom == "trans"').copy()

    dist_start = np.abs(cis_bins['start'] - gene_point)
    dist_end = np.abs(cis_bins['end'] - gene_point)
    dist_mid = np.abs((cis_bins['start'] + cis_bins['end']) // 2 - gene_point)

    # Calculate scalings

    B, C = _get_B_C(D_high, D_low, F_high, F_low)

    log_scaling_start = _sigmoid(np.log10(dist_start), A, B, C)
    log_scaling_end = _sigmoid(np.log10(dist_end), A, B, C)
    log_scaling_mid = _sigmoid(np.log10(dist_mid), A, B, C)

    log_scaling_avg = (log_scaling_start + log_scaling_end + log_scaling_mid) / 3
    scaling_factors = 10 ** log_scaling_avg
    cis_bins['count'] = cis_bins['count'] * scaling_factors

    # Normalize and fit backgrounds for simulations of non-specific contacts

    cis_simulator = _get_simulator_from_bins(cis_bins)
    trans_simulator = _get_simulator_from_bins(trans_bins)

    # Get non-specific contacts

    cis_non_specif = cis_simulator.rvs(size=int(N_total * frac_cis)).round().astype('int')
    trans_non_specif = trans_simulator.rvs(size=int(N_total * (1 - frac_cis))).round().astype('int')

    cis_non_specif_record = "\n".join((_point_to_bed(item, 'cis', rna_name)
                                       for item in cis_non_specif))
    trans_non_specif_record = "\n".join((_point_to_bed(item, 'trans', rna_name)
                                         for item in trans_non_specif))

    with open(outdir / 'non_specific_contacts.bed', 'w') as outfile:
        outfile.write('\n'.join((cis_non_specif_record, trans_non_specif_record)) + '\n')

    # Sample peak summits for peaks

    n_cis_peaks = round(N_peaks * frac_cis_peaks)
    n_trans_peaks = N_peaks - n_cis_peaks

    cis_peak_summits, cis_peak_abundances = _sample_peak_summits(cis_simulator, peak_exp, n_cis_peaks)
    trans_peak_summits, trans_peak_abundances = _sample_peak_summits(trans_simulator, peak_exp, n_trans_peaks)

    # Save peak summits

    summit_records = list()
    for summit, abundance in zip(cis_peak_summits, cis_peak_abundances):
        summit_records.append(_point_to_bed(summit, 'cis', rna_name, abundance))

    for summit, abundance in zip(trans_peak_summits, trans_peak_abundances):
        summit_records.append(_point_to_bed(summit, 'trans', rna_name, abundance))

    with open(outdir / 'peak_summits.bed', 'w') as outfile:
        outfile.write('\n'.join(summit_records) + '\n')

    # Sample coordinates of specific contacts

    cis_specif_coords = np.concatenate([_sample_specific_contacts(cis_summit, peak_exp, sd, L_cis)
                                        for cis_summit in cis_peak_summits])

    trans_specif_coords = np.concatenate([_sample_specific_contacts(trans_summit, peak_exp, sd, L_trans)
                                          for trans_summit in trans_peak_summits])

    # Save specific contacts

    cis_specif_record = "\n".join((_point_to_bed(item, 'cis', rna_name)
                                   for item in cis_specif_coords))
    trans_specif_record = "\n".join((_point_to_bed(item, 'trans', rna_name)
                                     for item in trans_specif_coords))

    with open(outdir / 'specific_contacts.bed', 'w') as outfile:
        outfile.write('\n'.join((cis_specif_record, trans_specif_record)) + '\n')

    # Save all contacts

    with open(outdir / 'contacts.bed', 'w') as outfile:
        outfile.write('\n'.join((cis_specif_record, trans_specif_record, cis_non_specif_record, trans_non_specif_record)) + '\n')

    # Voila!
