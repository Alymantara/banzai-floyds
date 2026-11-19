"""Quick-look checks of a reduction made with `banzai_floyds_offline`: what was reduced, the spectra, the accuracy
of the wavelength solutions, the flux accuracy of the standard stars, and the signal-to-noise.

Used by `notebooks/quickstart_reduction.ipynb`; every function works on the products in a processed folder:

    from banzai_floyds.diagnostics import quicklook as ql
    frames = ql.raw_frames('raw/')                 # what the raw folder contains
    spectra = ql.spectra('processed/')             # one row per extracted spectrum (-1d)
    arcs = ql.arc_accuracy('processed/')           # rms of the arc-line residuals per order
    std = ql.standard_flux_accuracy('processed/')  # standards against their reference spectra
"""
import glob
import importlib.resources
import os

import numpy as np
from astropy.io import fits
from astropy.table import Table

from banzai_floyds.utils.flux_utils import standard_at_resolution

ORDERS = {1: 'red', 2: 'blue'}
# Telluric bands (A) left out of the flux comparisons: the O2 A and B bands and the water bands
TELLURIC_BANDS = ((6270, 6330), (6860, 7090), (7150, 7400), (7550, 7720), (7870, 8050), (8080, 8420), (8900, 9900))
MAD_TO_SIGMA = 1.4826


def _robust_rms(values) -> float:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return np.nan
    return float(MAD_TO_SIGMA * np.median(np.abs(values - np.median(values))))


def _header(path: str) -> fits.Header:
    with fits.open(path) as hdul:
        return hdul[0].header if 'OBJECT' in hdul[0].header else hdul[1].header


def raw_frames(raw_dir: str) -> Table:
    """One row per raw frame: file, night, type (ARC, LAMPFLAT, SPECTRUM), object, slit width and exposure time."""
    rows = []
    for path in sorted(glob.glob(os.path.join(raw_dir, '*.fits*'))):
        header = fits.getheader(path, 'SCI') if 'SCI' in [h.name for h in fits.open(path)] else fits.getheader(path)
        rows.append({'file': os.path.basename(path), 'night': os.path.basename(path).split('-')[2],
                     'type': header.get('OBSTYPE'), 'object': header.get('OBJECT'),
                     'slit': float(header.get('APERWID', np.nan)), 'exptime': float(header.get('EXPTIME', np.nan))})
    return Table(rows=rows) if rows else Table()


def products(processed_dir: str, kind: str = '1d') -> list:
    """Paths of the processed products: kind '1d' or '2d' for spectra (e91), 'arc' for arcs (a91)."""
    pattern = {'1d': '*-e91-1d.fits*', '2d': '*-e91-2d.fits*', 'arc': '*-a91.fits*'}[kind]
    return sorted(glob.glob(os.path.join(processed_dir, '**', pattern), recursive=True))


def standard_names() -> dict:
    """The flux standards known to the pipeline: normalised object name -> reference spectrum file."""
    directory = os.path.join(importlib.resources.files('banzai_floyds'), 'data', 'standards')
    return {os.path.basename(p)[:-5].replace('+', '').replace('-', '').lower(): p
            for p in glob.glob(os.path.join(directory, '*.fits'))}


def _standard_file(name: str):
    key = name.replace('+', '').replace('-', '').replace(' ', '').lower()
    known = standard_names()
    if key in known:
        return known[key]
    # BD+28 4211 is stored as bdp284211
    return known.get(key.replace('bd', 'bdp', 1))


def spectra(processed_dir: str) -> Table:
    """One row per extracted spectrum, with its object, slit, airmass, whether it is a flux standard, and the
    quality keywords of the pipeline: median S/N per bin in the red (6000-8800 A) and blue (3800-5400 A) orders
    (L1QSN1/2), the seeing measured in the slit (L1QSEE1, arcsec), the object wavelength zero point applied
    (L1WZP1/2, A) and the quality flags (L1QFLAGS); for science frames, how the standards entered: the number of
    standards in the ensemble flux calibration (L1FLXNST) and the spread of their grey factors (max/min, L1FLXGRY),
    the number in the telluric correction (L1TELNST), the days to the nearest of them (L1QSTDDT), and whether the
    residual fringe template of the slit was applied (L1FRGTPL)."""
    rows = []
    for path in products(processed_dir, '1d'):
        header = _header(path)
        name = str(header.get('OBJECT', ''))
        rows.append({'file': os.path.basename(path), 'date': str(header.get('DATE-OBS', ''))[:16], 'object': name,
                     'standard': _standard_file(name) is not None, 'slit': header.get('APERWID'),
                     'airmass': header.get('AIRMASS'), 'exptime': header.get('EXPTIME'),
                     'snr_red': header.get('L1QSN1'), 'snr_blue': header.get('L1QSN2'),
                     'seeing': header.get('L1QSEE1'), 'zp_red': header.get('L1WZP1'),
                     'zp_blue': header.get('L1WZP2'), 'flags': header.get('L1QFLAGS', ''),
                     'flux_standards': header.get('L1FLXNST'), 'grey_spread': header.get('L1FLXGRY'),
                     'telluric_standards': header.get('L1TELNST'), 'days_to_standard': header.get('L1QSTDDT'),
                     'fringe_template': header.get('L1FRGTPL'), 'path': path})
    return Table(rows=rows) if rows else Table()


def extracted(path: str, order: int = None, good_only: bool = True) -> Table:
    """The extracted spectrum of a -1d product (optionally one order), sorted by wavelength; by default only the
    bins with mask == 0."""
    table = Table(fits.getdata(path, 'EXTRACTED'))
    if order is not None:
        table = table[table['order'] == order]
    if good_only:
        table = table[table['mask'] == 0]
    table.sort('wavelength')
    return table


def snr(path: str, order: int) -> tuple:
    """Wavelength and S/N per bin (optimally extracted flux over its error) of one order of a -1d product."""
    table = extracted(path, order)
    error = np.asarray(table['fluxerror_optimal'], dtype=float)
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = np.where(error > 0, np.asarray(table['flux_optimal'], dtype=float) / error, np.nan)
    return np.asarray(table['wavelength'], dtype=float), ratio


def _rejected_lines(residuals: np.ndarray) -> np.ndarray:
    """The lines the wavelength fit rejects as misidentified (wavelengths.fit_wavelength_solution): more than 5
    robust standard deviations, and more than 3 A, from the solution."""
    residuals = np.asarray(residuals, float)
    if residuals.size == 0:
        return np.zeros(0, bool)
    return np.abs(residuals - np.nanmedian(residuals)) > max(5 * _robust_rms(residuals), 3.0)


def arc_accuracy(processed_dir: str) -> Table:
    """For each arc, per order: the number of arc lines used in the wavelength solution, the number rejected as
    misidentified, and the robust rms and largest absolute value of the residuals of the lines used (measured minus
    catalogue wavelength, A): the typical statistical accuracy of the solution and its worst line."""
    rows = []
    for path in products(processed_dir, 'arc'):
        residuals = arc_residuals(path)
        for order, colour in ORDERS.items():
            these = residuals[(residuals['order'] == order) & residuals['used']]
            rejected = residuals[(residuals['order'] == order) & ~residuals['used']]
            rows.append({'arc': os.path.basename(path), 'order': colour, 'lines': len(these),
                         'rejected': ', '.join(f'{w:.0f}' for w in rejected['reference_wavelength']) or '-',
                         'rms_A': _robust_rms(these['residual']),
                         'max_abs_A': float(np.nanmax(np.abs(these['residual']))) if len(these) else np.nan})
    return Table(rows=rows) if rows else Table()


def arc_residuals(path: str) -> Table:
    """The residual of every identified arc line of one arc: order, catalogue wavelength, residual (A), and whether
    the wavelength fit used it (False: rejected as misidentified)."""
    table = Table(fits.getdata(path, 'RESIDUALS'))['order', 'reference_wavelength', 'residual']
    table['used'] = np.ones(len(table), bool)
    for order in ORDERS:
        these = table['order'] == order
        table['used'][these] = ~_rejected_lines(table['residual'][these])
    return table


def _lsf_sigma_angstrom(spectrum_path: str, order: int) -> float:
    """Gaussian sigma of the line spread function in A: the arc's sigma in pixels times the bin width."""
    arc_name = str(_header(spectrum_path).get('L1IDARC', ''))
    arc = glob.glob(os.path.join(os.path.dirname(spectrum_path), arc_name))
    table = extracted(spectrum_path, order, good_only=False)
    bin_width = float(np.nanmedian(table['binwidth']))
    sigma_pixels = fits.getheader(arc[0], 'LSF').get(f'SIGMA_{order}', 2.0) if arc else 2.0
    return sigma_pixels * bin_width


def flux_ratio_to_reference(path: str, order: int, sensitivity_from: str = None) -> tuple:
    """Wavelength and the ratio of a standard's calibrated flux to its reference spectrum (at the FLOYDS resolution)
    for one order. With `sensitivity_from`, the standard is calibrated with the sensitivity stored in another -1d
    product instead of its own (its own sensitivity reproduces the reference by construction, so this is the test
    that measures the flux accuracy a science frame gets)."""
    name = str(_header(path).get('OBJECT', ''))
    reference = Table(fits.getdata(_standard_file(name), 1))
    table = extracted(path, order)
    wavelength = np.asarray(table['wavelength'], dtype=float)
    flux = np.asarray(table['flux_optimal'], dtype=float)
    if sensitivity_from is not None:
        own = Table(fits.getdata(path, 'SENSITIVITY'))
        other = Table(fits.getdata(sensitivity_from, 'SENSITIVITY'))
        own, other = own[own['order'] == order], other[other['order'] == order]
        own.sort('wavelength')
        other.sort('wavelength')
        flux = flux * (np.interp(wavelength, other['wavelength'], other['sensitivity'])
                       / np.interp(wavelength, own['wavelength'], own['sensitivity']))
    expected = standard_at_resolution(wavelength, np.asarray(reference['wavelength'], dtype=float),
                                      np.asarray(reference['flux'], dtype=float), _lsf_sigma_angstrom(path, order))
    with np.errstate(divide='ignore', invalid='ignore'):
        return wavelength, flux / expected


def outside_telluric(wavelength) -> np.ndarray:
    keep = np.ones(np.shape(wavelength), bool)
    for low, high in TELLURIC_BANDS:
        keep &= ~((wavelength > low) & (wavelength < high))
    return keep


def standard_flux_accuracy(processed_dir: str) -> Table:
    """Each standard calibrated with the sensitivity of every other standard of the same slit width, compared with
    its reference spectrum outside the telluric bands, per order: median ratio (the grey offset, mostly clouds and
    slit losses) and robust rms of the ratio after removing that median (the error in the shape of the spectrum),
    both in per cent."""
    table = spectra(processed_dir)
    if len(table) == 0:
        return Table()
    standards = table[table['standard']]
    rows = []
    for row in standards:
        for other in standards:
            if other['path'] == row['path'] or other['slit'] != row['slit']:
                continue
            for order, colour in ORDERS.items():
                wavelength, ratio = flux_ratio_to_reference(row['path'], order, sensitivity_from=other['path'])
                good = outside_telluric(wavelength) & np.isfinite(ratio) & (ratio > 0)
                if good.sum() < 20:
                    continue
                median = float(np.median(ratio[good]))
                rows.append({'standard': row['file'], 'sensitivity_from': other['file'], 'order': colour,
                             'offset_pct': 100 * (median - 1), 'shape_rms_pct': 100 * _robust_rms(ratio[good] / median)})
    return Table(rows=rows) if rows else Table()


def recalibrated(science_path: str, standard_path: str, order: int) -> tuple:
    """Wavelength and the ratio of a science spectrum calibrated with one standard's own sensitivity to the same
    spectrum as the pipeline calibrated it (with the ensemble sensitivity stored in its SENSITIVITY extension).
    Everything else (extinction, telluric and fringe corrections) is the same, so the ratio is the effect of the
    choice of standard alone."""
    table = extracted(science_path, order)
    wavelength = np.asarray(table['wavelength'], dtype=float)
    used = Table(fits.getdata(science_path, 'SENSITIVITY'))
    single = Table(fits.getdata(standard_path, 'SENSITIVITY'))
    used, single = used[used['order'] == order], single[single['order'] == order]
    used.sort('wavelength')
    single.sort('wavelength')
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = (np.interp(wavelength, single['wavelength'], single['sensitivity'])
                 / np.interp(wavelength, used['wavelength'], used['sensitivity']))
    return wavelength, ratio


def science_vs_standards(processed_dir: str) -> Table:
    """For every science spectrum and every standard: how the science flux would change if it were calibrated with
    that standard alone instead of the ensemble the pipeline used, per order, outside the telluric bands. offset_pct
    is the overall change (a constant factor: clouds and slit losses of that standard's night, and the different
    slit width of the standards); shape_rms_pct the change in the shape of the spectrum after removing it;
    slope_pct the change from the blue to the red end of the order (difference of the median ratio in its first and
    last fifth)."""
    table = spectra(processed_dir)
    if len(table) == 0:
        return Table()
    rows = []
    for science in table[~table['standard']]:
        for standard in table[table['standard']]:
            for order, colour in ORDERS.items():
                wavelength, ratio = recalibrated(science['path'], standard['path'], order)
                good = outside_telluric(wavelength) & np.isfinite(ratio) & (ratio > 0)
                if good.sum() < 20:
                    continue
                median = float(np.median(ratio[good]))
                fifth = max(good.sum() // 5, 1)
                values = ratio[good] / median
                rows.append({'science': science['file'], 'object': science['object'], 'standard': standard['file'],
                             'order': colour, 'offset_pct': 100 * (median - 1),
                             'shape_rms_pct': 100 * _robust_rms(values),
                             'slope_pct': 100 * (np.median(values[-fifth:]) - np.median(values[:fifth]))})
    return Table(rows=rows) if rows else Table()
