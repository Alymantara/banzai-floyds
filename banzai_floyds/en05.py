"""COJ en05, the FLOYDS camera at Siding Spring until early 2019 (AVA B60).

There are no en05 sky flats in the LCO archive, so the order trace comes from a *hybrid* solution:

* order 1 (red) is traced on an en05 lamp flat processed as a sky flat (the lamp lights it well: it matches a real
  en12 sky flat to 0.35 px rms);
* order 2 (blue) has almost no lamp light below ~4600 A, so its shape is taken from an en12 sky-flat solution of the
  same slit width and shifted onto the en05 lamp light with banzai's matched filter (``fit_order_tweak``) over the
  lit columns (validated on en12 to ~1 px).

``register_en05`` puts the camera into a calibration database (instrument, order heights and locations, LSF
parameters copied from en12, and the shipped hybrid solutions in ``data/orders/skyflats_en05``); it is idempotent and
``setup_db`` calls it. ``build_hybrid`` makes a hybrid solution from one en05 lamp flat.
"""

import glob
import importlib.resources
import os
import shutil

import numpy as np

SITE, CAMERA, NAME, SOURCE_CAMERA = 'coj', 'en05', 'floyds02', 'en12'
#: columns of order 2 where an en05 / en12 lamp flat has light (lambda >~ 4700 A)
LIT_COLUMNS = (1300, 1900)


def shipped_solutions():
    directory = os.path.join(importlib.resources.files('banzai_floyds'), 'data', 'orders', 'skyflats_en05')
    return sorted(glob.glob(os.path.join(directory, '*-f91.fits*')))


def _instrument(db_address, camera):
    from banzai.dbs import query_for_instrument
    return query_for_instrument(db_address, SITE, camera)


def register_en05(db_address, solutions=None):
    """Add en05 to ``db_address`` (idempotent): the instrument, en12's order heights / locations and LSF parameters
    (same spectrograph and detector format), and the hybrid order solutions. Returns the number of solutions."""
    from banzai.dbs import add_instrument, get_session
    from banzai_floyds.dbs import LSFParams, OrderHeight, OrderLocation, register_skyflat
    add_instrument({'site': SITE, 'camera': CAMERA, 'name': NAME, 'type': '2m0-FLOYDS-SciCam',
                    'nx': 2079, 'ny': 512}, db_address)
    source, target = _instrument(db_address, SOURCE_CAMERA), _instrument(db_address, CAMERA)
    with get_session(db_address=db_address) as session:
        for table in (OrderHeight, OrderLocation, LSFParams):
            session.query(table).filter(table.instrument_id == target.id).delete()
            columns = [c.name for c in table.__table__.columns if c.name not in ('id', 'instrument_id')]
            for row in session.query(table).filter(table.instrument_id == source.id).all():
                session.add(table(instrument_id=target.id, **{c: getattr(row, c) for c in columns}))
        session.commit()
    paths = shipped_solutions() if solutions is None else list(solutions)
    for path in paths:
        register_skyflat(path, db_address)
    return len(paths)


def hybrid_coeffs(own, reference, shift):
    """ORDER_COEFFS table: order 1 from ``own`` (the lamp-flat trace), order 2 from ``reference`` shifted by
    ``shift`` px (added to the constant Legendre term)."""
    new = own.copy()
    i, j = list(own['order']).index(2), list(reference['order']).index(2)
    for column in [c for c in own.colnames if c.startswith('c')] + ['domainmin', 'domainmax', 'height']:
        new[column][i] = reference[column][j]
    new['c0'][i] += shift
    return new


def fit_heights(table, shape, margin=2):
    """Order heights reduced, if needed, so each order band (centre +- height / 2) stays on the detector with
    ``margin`` px to spare. en05's blue order sits high (centre y ~ 470 of 512 in 2016-2018), and banzai needs a
    full rectangle per order (``get_order_2d_region``); the rows lost are sky at the end of the slit."""
    from numpy.polynomial.legendre import Legendre
    new = table.copy()
    for i, row in enumerate(table):
        model = Legendre([row[c] for c in table.colnames if c.startswith('c')],
                         domain=(row['domainmin'], row['domainmax']))
        y = model(np.linspace(row['domainmin'], row['domainmax'], 200))
        room = min(np.min(y), shape[0] - 1 - np.max(y)) - margin
        new['height'][i] = int(min(row['height'], 2 * np.floor(room)))
    new['height'] = np.asarray(new['height']).astype(np.int64)          # as banzai writes it (K)
    return new


def order2_shift(data, error, mask, reference):
    """Shift (px) of the reference's order 2 that best matches the lamp light, over ``LIT_COLUMNS``."""
    from numpy.polynomial.legendre import Legendre
    from banzai_floyds.orders import fit_order_tweak
    row = reference[list(reference['order']).index(2)]
    coeffs = [row[c] for c in reference.colnames if c.startswith('c')]
    domain = (float(row['domainmin']), float(row['domainmax']))
    model = Legendre(coeffs, domain=domain)
    x2d, y2d = np.meshgrid(np.arange(data.shape[1]), np.arange(data.shape[0]))
    lo, hi = max(LIT_COLUMNS[0], domain[0]), min(LIT_COLUMNS[1], domain[1])
    region = (x2d >= lo) & (x2d <= hi) & (np.abs(y2d - model(x2d)) <= row['height']) & (mask == 0)
    shift, = fit_order_tweak(data[region], error[region], int(row['height']), coeffs,
                             (x2d[region], y2d[region]), domain)
    return float(shift)


def make_hybrid_file(lamp_f91, reference_f91, out_path):
    """Write the hybrid solution of a processed lamp flat (``lamp_f91``, with ORDER_COEFFS) to ``out_path``."""
    from astropy.io import fits
    from astropy.table import Table
    from numpy.polynomial.legendre import Legendre
    from banzai_floyds.orders import Orders
    reference = Table(fits.getdata(reference_f91, 'ORDER_COEFFS'))
    with fits.open(lamp_f91) as hdul:
        data, error = hdul['SCI'].data.astype(float), hdul['ERR'].data.astype(float)
        mask = hdul['BPM'].data if 'BPM' in hdul else np.zeros(data.shape, dtype=np.uint8)
        shift = order2_shift(data, error, mask, reference)
        new = fit_heights(hybrid_coeffs(Table(hdul['ORDER_COEFFS'].data), reference, shift), data.shape)
        models = [Legendre([r[c] for c in new.colnames if c.startswith('c')], domain=(r['domainmin'], r['domainmax']))
                  for r in new]
        hdul['ORDER_COEFFS'] = fits.BinTableHDU(new, header=hdul['ORDER_COEFFS'].header, name='ORDER_COEFFS')
        hdul['ORDERS'].data = Orders(models, data.shape, [int(r['height']) for r in new]).data
        header = hdul['SCI'].header
        header['L1ORDHYB'] = (os.path.basename(reference_f91)[:40], 'en05 hybrid: order 2 shape from this solution')
        header['L1ORDSH2'] = (round(shift, 3), '[px] shift of that order 2 onto the lamp light')
        header['OBSTYPE'] = 'SKYFLAT'
        tmp = out_path + '.tmp.fits'
        hdul.writeto(tmp, overwrite=True)
    shutil.move(tmp, out_path)
    return shift


def build_hybrid(lamp_raw, reference_f91, db_address, processed_dir, out_dir):
    """Process an en05 raw lamp flat as a sky flat (order 1 traced by OrderSolver), then write the hybrid solution
    to ``out_dir`` as ``<lamp stem>-f91.fits.fz``. Returns (path, order-2 shift)."""
    from astropy.io import fits
    from banzai.dbs import get_session
    from banzai.utils.stage_utils import run_pipeline_stages
    from banzai_floyds import offline
    from banzai_floyds.dbs import FLOYDSCalibrationImage
    context = offline._context_for(db_address, processed_dir)
    staging = os.path.join(processed_dir, 'input')
    os.makedirs(staging, exist_ok=True)
    stem = os.path.basename(lamp_raw).split('.fits')[0].replace('-w00', '-f00')
    staged = os.path.join(staging, stem + '.fits.fz')
    with fits.open(lamp_raw) as hdul:
        for hdu in hdul:
            if 'OBSTYPE' in hdu.header:
                hdu.header['OBSTYPE'] = 'SKYFLAT'
        hdul.writeto(staged, overwrite=True)
    run_pipeline_stages([{'path': staged}], context)
    f91 = stem.replace('-f00', '-f91') + '.fits.fz'
    found = glob.glob(os.path.join(processed_dir, '**', f91), recursive=True)
    if not found:
        raise RuntimeError(f'{os.path.basename(lamp_raw)}: no order solution from the lamp flat (see the log)')
    os.makedirs(out_dir, exist_ok=True)
    out = os.path.join(out_dir, f91)
    shift = make_hybrid_file(found[0], reference_f91, out)
    with get_session(db_address=db_address) as session:      # the lamp-only solution must not be used
        session.query(FLOYDSCalibrationImage).filter(FLOYDSCalibrationImage.filename == f91).delete()
        session.commit()
    return out, shift
