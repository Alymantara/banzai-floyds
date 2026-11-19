"""Run BANZAI-FLOYDS offline on a directory of raw frames, with a local (SQLite) calibration database.

    banzai_floyds_offline setup-db  --db-address sqlite:///floyds.db
    banzai_floyds_offline fetch-skyflats --raw-dir skyflats/ [--site coj]
    banzai_floyds_offline skyflats  --db-address ... --raw-dir skyflats/ --processed-dir processed/  (optional)
    banzai_floyds_offline run       --db-address ... --raw-dir raw/ --processed-dir processed/ --kind arcs
    banzai_floyds_offline chain     --db-address ... --raw-dir raw/ --processed-dir processed/

``chain`` runs the whole reduction of the improved pipeline in the order it needs: arcs -> lamp flats -> fringe
phase masters (one per slit width and ~60-day window) -> standards (pass 1) -> residual fringe templates ->
standards (pass 2) -> science. Frames are classified from their headers: ARC, LAMPFLAT, and SPECTRUM split into
standards (a position in the flux-standard table) and science.
"""
import argparse
import datetime
import glob
import os
import time
import traceback
from multiprocessing import Pool

os.environ.setdefault('OPENTSDB_PYTHON_METRICS_TEST_MODE', 'True')

SITES = {'coj': ('-31.272932', '149.070648', '1116', '10'), 'ogg': ('20.7069444444', '-156.258055556', '3065', '-10')}
INSTRUMENTS = {'coj': ('en12', 'floyds02'), 'ogg': ('en06', 'floyds01')}
KINDS = ('arcs', 'flats', 'standards', 'science')
_context = None


def setup_db(db_address):
    """Create the calibration database with both FLOYDS sites and instruments, order locations, LSF parameters, the
    flux standards and the sky-flat order solutions shipped with the package."""
    import subprocess
    import sys
    from banzai_floyds import dbs
    bindir = os.path.dirname(sys.executable)
    os.environ['DB_ADDRESS'] = db_address

    def run(*args):
        subprocess.run([os.path.join(bindir, args[0]), *args[1:], f'--db-address={db_address}'], check=True)
    run('banzai_floyds_create_db')
    for site, (lat, lon, elev, tz) in SITES.items():
        run('banzai_add_site', '--site', site, '--latitude', lat, '--longitude', lon, '--elevation', elev,
            '--timezone', tz)
    for site, (camera, name) in INSTRUMENTS.items():
        run('banzai_add_instrument', '--site', site, '--camera', camera, '--name', name,
            '--instrument-type', '2m0-FLOYDS-SciCam', '--nx', '2079', '--ny', '512')
    dbs.populate_order_heights_locations(db_address)
    dbs.populate_lsf_params(db_address)
    dbs.ingest_standards(db_address)
    # the order solutions of all the sky flats in data/orders/skyflats.dat ship with the package
    print('sky flats registered:', dbs.ingest_skyflats(db_address))


def db_is_set_up(db_address):
    """True if the database exists and setup_db has filled it (an empty file, left by a crashed run, is not)."""
    from sqlalchemy import create_engine, inspect
    if db_address.startswith('sqlite:///') and not os.path.exists(db_address[len('sqlite:///'):]):
        return False
    engine = create_engine(db_address)
    try:
        return {'fluxstandards', 'sites'} <= set(inspect(engine).get_table_names())
    finally:
        engine.dispose()


def _context_for(db_address, processed_dir):
    from banzai_floyds import settings
    import banzai.main
    os.environ['DB_ADDRESS'] = db_address
    settings.processed_path = processed_dir
    settings.fpack = True
    settings.db_address = db_address
    return banzai.main.parse_args(settings, parse_system_args=False)


def _init(db_address, processed_dir):
    global _context
    _context = _context_for(db_address, processed_dir)


def _process(path):
    from banzai.utils.stage_utils import run_pipeline_stages
    start = time.time()
    # where SameVisitCosmicRays looks for the other exposures of the visit (banzai keeps only the file name)
    from banzai_floyds import visit_cosmics
    visit_cosmics.RAW_DIRECTORY = os.path.dirname(os.path.abspath(path))
    for attempt in range(3):
        try:
            run_pipeline_stages([{'path': path}], _context)
            return path, 'ok', time.time() - start
        except Exception as error:  # SQLite lock contention between workers: retry
            if 'locked' in str(error) and attempt < 2:
                time.sleep(5)
                continue
            return path, 'ERROR ' + ''.join(traceback.format_exception_only(error)).strip()[:300], time.time() - start


def _register_processed_skyflat(path, db_address, processed_dir):
    """Copy an archive-processed sky flat (-f91, with ORDERS and ORDER_COEFFS) into the processed tree and register
    it as an order solution, as if it had been processed here."""
    import shutil
    from banzai_floyds.dbs import register_skyflat
    from astropy.io import fits
    header = fits.getheader(path, 'SCI')
    out_dir = os.path.join(processed_dir, header['SITEID'].strip(), header['INSTRUME'].strip(),
                           os.path.basename(path).split('-')[2], 'processed')
    os.makedirs(out_dir, exist_ok=True)
    out = os.path.join(out_dir, os.path.basename(path))
    shutil.copy(path, out)
    register_skyflat(out, db_address)


ARCHIVE_API = 'https://archive-api.lco.global/frames/'


def _archive_get(url, max_wait=120.0):
    """GET from the LCO archive, waiting out short 429 throttles (Retry-After) and giving up on long ones."""
    import urllib.error
    import urllib.request
    request = urllib.request.Request(url, headers={'User-Agent': 'banzai-floyds-offline'})
    for _ in range(5):
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                return response.read()
        except urllib.error.HTTPError as error:
            if error.code != 429:
                raise
            wait = float(error.headers.get('Retry-After', 60))
            if wait > max_wait:
                raise RuntimeError(f'LCO archive throttled for {wait:.0f} s; try again later') from error
            time.sleep(wait)
    raise RuntimeError('LCO archive kept throttling')


def fetch_skyflats(raw_dir, site=None, pause=1.0):
    """Download the raw sky flats listed in data/orders/skyflats.dat (public frames, by archive frame id) that are
    not in raw_dir yet. One request a second: the archive throttles bursts for hours."""
    import importlib.resources
    import json
    from astropy.io import ascii
    table = ascii.read(os.path.join(importlib.resources.files('banzai_floyds'), 'data', 'orders', 'skyflats.dat'))
    os.makedirs(raw_dir, exist_ok=True)
    present = {os.path.basename(f).split('-x00')[0].split('-f91')[0] for f in glob.glob(os.path.join(raw_dir, '*'))}
    fetched = 0
    for row in table:
        if site and row['site'] != site:
            continue
        root = row['filename'].split('-x00')[0]
        if root in present:
            continue
        frame = json.loads(_archive_get(f'{ARCHIVE_API}{row["frameid"]}/'))
        time.sleep(pause)
        if frame.get('basename', '').split('-x00')[0] != root:
            raise RuntimeError(f'archive frame {row["frameid"]} is {frame.get("basename")}, not {row["filename"]}')
        data = _archive_get(frame['url'])
        out = os.path.join(raw_dir, row['filename'])
        with open(out + '.part', 'wb') as handle:
            handle.write(data)
        os.replace(out + '.part', out)
        print('fetched', out, flush=True)
        fetched += 1
        time.sleep(pause)
    return fetched


def process_skyflats(db_address, raw_dir, processed_dir):
    """Relabel raw sky flats as SKYFLAT, process them, and bound their validity windows (data/orders/skyflats.dat).
    Sky flats already processed by the archive (-f91 files with an ORDER_COEFFS extension) are registered as they
    are."""
    from astropy.io import fits
    from banzai.utils.stage_utils import run_pipeline_stages
    from banzai_floyds.dbs import bound_skyflats_to_windows
    context = _context_for(db_address, processed_dir)
    staging = os.path.join(processed_dir, 'input')
    os.makedirs(staging, exist_ok=True)
    for path in sorted(glob.glob(os.path.join(raw_dir, '*.fits*'))):
        with fits.open(path) as hdul:
            processed = 'ORDER_COEFFS' in hdul
        if processed:
            _register_processed_skyflat(path, db_address, processed_dir)
            print('registered', path, flush=True)
            continue
        with fits.open(path) as hdul:
            hdul['SCI'].header['OBSTYPE'] = 'SKYFLAT'
            out = os.path.join(staging, os.path.basename(path).replace('x00.fits', 'f00.fits'))
            hdul.writeto(out, overwrite=True)
        run_pipeline_stages([{'path': out}], context)
        print('processed', out, flush=True)
    print('sky-flat windows bounded:', bound_skyflats_to_windows(db_address))


def classify(raw_dir, db_address):
    """{kind: [paths sorted by DATE-OBS]} from the headers of the raw frames."""
    from astropy.io import fits
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    from banzai_floyds.dbs import find_flux_standard
    groups = {kind: [] for kind in KINDS}
    for path in sorted(glob.glob(os.path.join(raw_dir, '*.fits*'))):
        header = fits.getheader(path, 'SCI') if 'SCI' in [h.name for h in fits.open(path)] else fits.getheader(path)
        obstype = header.get('OBSTYPE', '')
        if obstype == 'ARC':
            kind = 'arcs'
        elif obstype == 'LAMPFLAT':
            kind = 'flats'
        elif obstype == 'SPECTRUM':
            # the catalogue position, as BANZAI uses (the pointing can be several arcsec off)
            ra, dec = header.get('CAT-RA', header['RA']), header.get('CAT-DEC', header['DEC'])
            position = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
            standard = find_flux_standard(position.ra.deg, position.dec.deg, db_address)
            kind = 'standards' if standard is not None else 'science'
        else:
            continue
        groups[kind].append((header.get('DATE-OBS', ''), path))
    return {kind: [p for _, p in sorted(items)] for kind, items in groups.items()}


def run(paths, db_address, processed_dir, workers=4, record=None):
    """Reduce `paths` with `workers` processes; return the number of failures. With `record` (a list), append
    (file, status, seconds) for every frame."""
    failures = 0
    with Pool(workers, initializer=_init, initargs=(db_address, processed_dir)) as pool:
        for path, status, seconds in pool.imap_unordered(_process, paths):
            print(f'{os.path.basename(path)} {status} {seconds:.0f}s', flush=True)
            failures += status != 'ok'
            if record is not None:
                record.append((os.path.basename(path), status, float(seconds)))
    return failures


def _windows(db_address, days=60):
    """Consecutive date windows of about `days` covering the processed lamp flats."""
    from banzai.dbs import get_session
    from banzai_floyds.dbs import FLOYDSCalibrationImage
    with get_session(db_address=db_address) as session:
        dates = [r.dateobs for r in session.query(FLOYDSCalibrationImage).filter(
            FLOYDSCalibrationImage.type == 'LAMPFLAT', FLOYDSCalibrationImage.is_master.is_(False)).all()]
    if not dates:
        return []
    start, end = min(dates), max(dates) + datetime.timedelta(days=1)
    n = max(int(round((end - start).days / days)), 1)
    edges = [start + (end - start) * i / n for i in range(n + 1)]
    return [(edges[i].isoformat(), edges[i + 1].isoformat()) for i in range(n)]


def chain(raw_dir, db_address, processed_dir, workers=4, window_days=60, timings=None):
    """Reduce everything in `raw_dir` in the order the improved pipeline needs; return the number of science
    frames that failed. With `timings` (a dict), record the wall-clock seconds and frame count of each step in
    timings['steps'] and the time of every frame in timings['frames'] (step, file, status, seconds)."""
    from banzai_floyds.fringe_masters import keep_template_free, make_phase_masters, make_residual_templates
    steps, frames = [], []
    if timings is not None:
        timings['steps'], timings['frames'] = steps, frames

    def timed(name, paths=None, function=None):
        start, record = time.time(), []
        result = run(paths, db_address, processed_dir, workers, record) if paths is not None else function()
        steps.append({'step': name, 'frames': len(paths) if paths is not None else 0,
                      'seconds': time.time() - start})
        frames.extend({'step': name, 'file': f, 'status': st, 'seconds': sec} for f, st, sec in record)
        return result

    groups = classify(raw_dir, db_address)
    print({kind: len(paths) for kind, paths in groups.items()}, flush=True)
    timed('arcs', groups['arcs'])
    timed('lamp flats', groups['flats'])

    def phase_masters():
        for low, high in _windows(db_address, window_days):
            for path in make_phase_masters(db_address, processed_dir, low, high):
                print('phase master', path, flush=True)
        make_residual_templates(db_address, strip=True)
    timed('fringe phase masters', function=phase_masters)
    timed('standards, first pass', groups['standards'])

    def templates():
        print('template-free standards kept', keep_template_free(db_address), flush=True)
        print('residual templates', make_residual_templates(db_address), flush=True)
    timed('residual fringe templates', function=templates)
    timed('standards, second pass', groups['standards'])
    return timed('science', groups['science'])


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('command', choices=['setup-db', 'fetch-skyflats', 'skyflats', 'run', 'chain'])
    parser.add_argument('--db-address', help='e.g. sqlite:///floyds.db (all commands but fetch-skyflats)')
    parser.add_argument('--site', choices=sorted(SITES), help='for fetch-skyflats: one site only')
    parser.add_argument('--raw-dir')
    parser.add_argument('--processed-dir')
    parser.add_argument('--kind', choices=KINDS, help='for run: which frames')
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--window-days', type=float, default=60.0, help='for chain: phase-master window')
    args = parser.parse_args(argv)
    if args.command == 'fetch-skyflats':
        print('sky flats fetched:', fetch_skyflats(args.raw_dir or 'skyflats', args.site))
        return
    if not args.db_address:
        parser.error('--db-address is required')
    if args.command == 'setup-db':
        setup_db(args.db_address)
    elif args.command == 'skyflats':
        process_skyflats(args.db_address, args.raw_dir, args.processed_dir)
    elif args.command == 'run':
        paths = classify(args.raw_dir, args.db_address)[args.kind]
        raise SystemExit(run(paths, args.db_address, args.processed_dir, args.workers) > 0)
    else:
        raise SystemExit(chain(args.raw_dir, args.db_address, args.processed_dir, args.workers,
                               args.window_days) > 0)


if __name__ == '__main__':
    main()
