# %% [markdown]
# # FLOYDS quick start: from raw frames to checked spectra
#
# This notebook reduces a folder of raw FLOYDS frames with the improved BANZAI-FLOYDS pipeline (branch
# `rm-improvements`) on your own computer, and then checks the result:
#
# 1. what the raw folder contains;
# 2. the calibration database and the whole reduction (arcs, lamp flats, fringe phase masters, standards, science),
#    and how long each step took;
# 3. the extracted spectra;
# 4. how accurate the wavelength solutions are;
# 5. how accurate the flux calibration of the standard stars is;
# 6. the signal-to-noise of each spectrum;
# 7. how the standard stars affect the calibrated science spectra.
#
# **What you need.** The pipeline installed as in `docs/INSTALL_RM.md`, plus Jupyter in the same environment
# (`pip install jupyterlab`). A folder of raw frames (`*-e00`, `*-a00`, `*-w00` `.fits.fz` files from the LCO
# archive) with, for each night, the science exposures and the arc and lamp flat of their block, and some flux
# standards (stars in `banzai_floyds/data/standards`) taken within a few weeks. Sky flats are not needed: their
# order solutions come with the package.
#
# **Time.** About 4 minutes on an 8-core laptop for the 21-frame example (three nights of Fairall 9 and Feige 110
# at COJ); about 1 minute per visit and 4 workers in general. See `docs/QUICKSTART.md` for the command-line version.

# %%
import os

os.environ.setdefault('OPENTSDB_PYTHON_METRICS_TEST_MODE', 'True')  # BANZAI must not try to send metrics

import matplotlib.pyplot as plt
import numpy as np

from banzai_floyds import offline
from banzai_floyds.diagnostics import quicklook as ql

# Edit these lines for your data (or set the environment variables FLOYDS_RAW and FLOYDS_WORK)
RAW_DIR = os.path.expanduser(os.environ.get('FLOYDS_RAW', 'raw'))    # the raw frames
WORK_DIR = os.path.expanduser(os.environ.get('FLOYDS_WORK', '.'))    # the database and processed products go here
WORKERS = 4                                                         # frames reduced in parallel (half your cores)
REDUCE = True                                               # False: only look at products already reduced

DB = f"sqlite:///{os.path.join(WORK_DIR, 'floyds.db')}"
PROCESSED = os.path.join(WORK_DIR, 'processed')

# %% [markdown]
# ## 1. What is in the raw folder
# Each block (one telescope request) should have an `ARC`, a `LAMPFLAT` and its `SPECTRUM` exposures with the
# same slit width.

# %%
frames = ql.raw_frames(RAW_DIR)
print(f'{len(frames)} raw frames')
frames['night', 'type', 'object', 'slit', 'exptime', 'file'].pprint(max_lines=-1, max_width=-1)

# %% [markdown]
# ## 2. Reduce
# `setup_db` creates a local calibration database (an SQLite file) with both FLOYDS sites, the order positions,
# the line-spread-function seeds, the flux standards and the sky-flat order solutions shipped with the package.
# `chain` then reduces every frame in the order the pipeline needs. It reduces all the frames each time it runs
# (about 4 minutes for the example), so set `REDUCE = False` above to look at the results again without reducing.
# Progress is logged below, one line per frame with `ok` or the error.

# %%
import json
import time

timings = {}
if REDUCE:
    start = time.time()
    if not offline.db_is_set_up(DB):
        offline.setup_db(DB)
    timings['setup_seconds'] = time.time() - start
    n_failed = offline.chain(RAW_DIR, DB, PROCESSED, workers=WORKERS, timings=timings)
    timings['total_seconds'] = time.time() - start
    timings['workers'] = WORKERS
    json.dump(timings, open(os.path.join(WORK_DIR, 'timings.json'), 'w'), indent=1)
    print('science frames that failed:', n_failed)
elif os.path.exists(os.path.join(WORK_DIR, 'timings.json')):
    timings = json.load(open(os.path.join(WORK_DIR, 'timings.json')))

# %% [markdown]
# ### How long it took
# Wall-clock time of each step of the chain (with `WORKERS` frames in parallel), and the median time one frame of
# each kind takes on one worker. Arcs take longest: every arc line is fitted in every detector row along the slit
# to measure the tilt and the shape of the lines. Standards are reduced twice, before and after the residual fringe
# template of their slit is measured from them.

# %%
from astropy.table import Table

if timings:
    steps = Table(rows=timings['steps'])
    steps['seconds'].format = '.0f'
    steps.pprint(max_width=-1)
    print(f"set-up {timings['setup_seconds']:.0f} s; total {timings['total_seconds']:.0f} s "
          f"({timings['total_seconds'] / 60:.1f} min) with {timings['workers']} workers")
    per_frame = Table(rows=timings['frames'])
    for step in dict.fromkeys(per_frame['step']):
        these = per_frame[per_frame['step'] == step]
        print(f"{step:24s} {len(these):3d} frames, {np.median(these['seconds']):5.1f} s per frame (median, one worker)")
    n_visits = int(np.sum(ql.raw_frames(RAW_DIR)['type'] == 'ARC'))
    print(f"{n_visits} blocks (one arc each): {timings['total_seconds'] / n_visits:.0f} s of wall-clock time per block")

# %% [markdown]
# ## 3. The spectra
# One row per extracted spectrum (`-1d` product). `snr_red`/`snr_blue` are the median signal-to-noise per bin in
# the red (6000–8800 Å) and blue (3800–5400 Å) orders, `seeing` the width of the object along the slit (arcsec),
# `zp_red`/`zp_blue` the wavelength zero point measured on the telluric bands in the object's own light (Å), and
# `flags` the quality flags (`docs/QC_KEYWORDS.md`).

# %%
spec = ql.spectra(PROCESSED)
for column in ('airmass', 'snr_red', 'snr_blue', 'seeing', 'zp_red', 'zp_blue'):
    spec[column].format = '.2f'
spec['date', 'object', 'standard', 'slit', 'airmass', 'snr_red', 'snr_blue', 'seeing', 'zp_red', 'zp_blue',
     'flags'].pprint(max_lines=-1, max_width=-1)

# %%
# Flux-calibrated spectra, both orders, one panel per object; masked bins (cosmic rays, the O2 A band) left out
objects = sorted(set(spec['object']))
fig, axes = plt.subplots(len(objects), 1, figsize=(10, 3.2 * len(objects)), squeeze=False)
for ax, name in zip(axes[:, 0], objects):
    for row in spec[spec['object'] == name]:
        for order, colour in ((2, 'tab:blue'), (1, 'tab:red')):
            table = ql.extracted(row['path'], order)
            ax.plot(table['wavelength'], table['flux_optimal'] * 1e15, lw=0.6, color=colour, alpha=0.6)
    ax.set_title(f"{name} ({np.sum(spec['object'] == name)} spectra; blue and red orders)", loc='left')
    ax.set_ylabel('flux (10$^{-15}$ erg s$^{-1}$ cm$^{-2}$ Å$^{-1}$)')
    ax.set_xlim(3200, 10000)
    top = np.nanpercentile(np.concatenate([ql.extracted(p, 1)['flux_optimal'] for p in
                                           spec[spec['object'] == name]['path']]), 99.5) * 1e15
    ax.set_ylim(0, 1.6 * top)
axes[-1, 0].set_xlabel('observed wavelength (Å)')
fig.tight_layout()

# %% [markdown]
# ## 4. Wavelength accuracy
# Each arc gives a wavelength solution per order from the arc lines it detects. The robust rms of the residuals
# (measured minus catalogue wavelength) is the statistical accuracy of the solution, and `max_abs_A` its worst line.
# A line matched to the wrong feature (in 2\" arcs, at times the Ar 6753/6767 blend or Ar 6416) lands several Å off;
# the fit rejects it, and `rejected` lists it (catalogue wavelength). The object's own zero point (`zp_red`/`zp_blue` above) then corrects for
# where the object sits in a wide slit, which the arc cannot know.

# %%
arcs = ql.arc_accuracy(PROCESSED)
arcs['rms_A'].format = '.3f'
arcs['max_abs_A'].format = '.2f'
arcs.pprint(max_lines=-1, max_width=-1)
for colour in ('red', 'blue'):
    these = arcs[arcs['order'] == colour]
    print(f'{colour} order: median rms {np.nanmedian(these["rms_A"]):.3f} Å over {len(these)} arcs, '
          f'{int(np.median(these["lines"]))} lines per arc')

# %%
# Residuals of every line of every arc against wavelength
fig, ax = plt.subplots(figsize=(10, 3.5))
for path in ql.products(PROCESSED, 'arc'):
    residuals = ql.arc_residuals(path)
    for order, colour in ((2, 'tab:blue'), (1, 'tab:red')):
        these = residuals[residuals['order'] == order]
        used = these[these['used']]
        ax.plot(used['reference_wavelength'], used['residual'], 'o', ms=3, color=colour, alpha=0.6)
        rejected = these[~these['used']]
        ax.plot(rejected['reference_wavelength'], np.clip(rejected['residual'], -1.9, 1.9), 'x', ms=6, color=colour)
ax.axhline(0, color='k', lw=0.6)
ax.set_ylim(-2, 2)
ax.set_xlabel('catalogue wavelength of the arc line (Å)')
ax.set_ylabel('measured − catalogue (Å)')
ax.set_title('Arc-line residuals, all arcs (blue and red orders; × rejected lines, drawn at the edge)', loc='left')
fig.tight_layout()

# %% [markdown]
# ## 5. Flux accuracy of the standard stars
# A standard calibrated with its own sensitivity reproduces its reference spectrum by construction. The honest
# test is to calibrate each standard with the sensitivity of another standard of the same slit width (another
# night), as a science frame is, and compare with the reference spectrum outside the telluric bands. `offset_pct`
# is the overall difference (mostly clouds and slit losses, which change from night to night); `shape_rms_pct` is
# the scatter of the ratio after removing that offset, the error in the shape of the calibrated spectrum.

# %%
std = ql.standard_flux_accuracy(PROCESSED)
if len(std):
    std['offset_pct'].format = '+.1f'
    std['shape_rms_pct'].format = '.1f'
    std.pprint(max_lines=-1, max_width=-1)
    for colour in ('red', 'blue'):
        these = std[std['order'] == colour]
        print(f'{colour} order: shape error {np.median(these["shape_rms_pct"]):.1f}% (median), '
              f'offsets {np.min(these["offset_pct"]):+.1f} to {np.max(these["offset_pct"]):+.1f}%')
else:
    print('Fewer than two standards of the same slit width: no transfer test possible.')

# %%
# The ratio to the reference for each standard calibrated with another night's sensitivity
standards = spec[spec['standard']]
if len(standards) > 1:
    fig, ax = plt.subplots(figsize=(10, 3.5))
    for row in standards:
        other = [o for o in standards if o['path'] != row['path'] and o['slit'] == row['slit']][0]
        for order, colour in ((2, 'tab:blue'), (1, 'tab:red')):
            wavelength, ratio = ql.flux_ratio_to_reference(row['path'], order, sensitivity_from=other['path'])
            ratio = np.where(ql.outside_telluric(wavelength), ratio, np.nan)
            ax.plot(wavelength, ratio, lw=0.7, color=colour, alpha=0.7)
    ax.axhline(1, color='k', lw=0.6)
    ax.set_ylim(0.8, 1.2)
    ax.set_xlim(3200, 10000)
    ax.set_xlabel('observed wavelength (Å)')
    ax.set_ylabel('calibrated / reference')
    ax.set_title('Standards calibrated with another night\'s sensitivity (telluric bands left out)', loc='left')
    fig.tight_layout()

# %% [markdown]
# ## 6. Signal-to-noise
# The signal-to-noise per bin (about one pixel: 1.7 Å blue, 3.5 Å red) of each spectrum, from the propagated
# errors of the optimal extraction. The table in section 3 has the medians (`snr_red`, `snr_blue`).

# %%
fig, ax = plt.subplots(figsize=(10, 3.5))
for row in spec:
    style = '-' if row['standard'] else '--'
    for order, colour in ((2, 'tab:blue'), (1, 'tab:red')):
        wavelength, ratio = ql.snr(row['path'], order)
        # a running median over 25 bins keeps the curve readable
        smooth = np.array([np.nanmedian(ratio[max(i - 12, 0):i + 13]) for i in range(len(ratio))])
        ax.plot(wavelength, smooth, style, lw=0.8, color=colour, alpha=0.6)
ax.set_yscale('log')
ax.set_xlim(3200, 10000)
ax.set_xlabel('observed wavelength (Å)')
ax.set_ylabel('S/N per bin')
ax.set_title('Signal-to-noise per bin (solid: standards, dashed: science; blue and red orders)', loc='left')
fig.tight_layout()

# %% [markdown]
# ## 7. How the standard stars affect the science spectra
# The standards enter a science spectrum in three ways:
#
# * **Flux calibration.** The science frame is calibrated with the *ensemble* sensitivity: the median of the
#   sensitivities of the (up to nine) standards nearest in time, each scaled by a constant (grey) factor, with the
#   overall level taken from the clearest of them. Below, every science spectrum is calibrated again with each
#   standard alone, everything else unchanged, to show how much the choice of standard matters.
# * **Telluric correction.** The depth of the atmospheric absorption bands is measured on the nearest standards,
#   scaled to the airmass of the science frame; it changes the spectrum only inside the bands.
# * **Residual fringe template.** The small fringe residual left by the fringe model is measured on the standards
#   of each slit width and removed from the science spectra of the same slit. Here the standards are taken with the
#   2 arcsec slit and the science with the 6 arcsec slit, so no template is applied to the science.

# %%
# How the standards were used for each science frame (header keywords of the pipeline)
science_rows = spec[~spec['standard']]
science_rows['date', 'object', 'flux_standards', 'grey_spread', 'telluric_standards', 'days_to_standard',
             'fringe_template'].pprint(max_width=-1)

# %%
influence = ql.science_vs_standards(PROCESSED)
if len(influence):
    for column in ('offset_pct', 'shape_rms_pct', 'slope_pct'):
        influence[column].format = '+.1f'
    influence['science', 'standard', 'order', 'offset_pct', 'shape_rms_pct', 'slope_pct'].pprint(max_lines=20,
                                                                                                max_width=-1)
    for colour in ('red', 'blue'):
        these = influence[influence['order'] == colour]
        print(f"{colour} order, one standard instead of the ensemble: overall level "
              f"{np.min(these['offset_pct']):+.1f} to {np.max(these['offset_pct']):+.1f}%, shape "
              f"{np.median(these['shape_rms_pct']):.1f}% rms (median; max {np.max(these['shape_rms_pct']):.1f}%), "
              f"blue-to-red tilt {np.min(these['slope_pct']):+.1f} to {np.max(these['slope_pct']):+.1f}%")

# %%
# One science spectrum calibrated with each standard alone, divided by the pipeline's calibration
science = spec[~spec['standard']]
if len(science) and np.any(spec['standard']):
    fig, ax = plt.subplots(figsize=(10, 3.5))
    for number, standard in enumerate(spec[spec['standard']]):
        for order in (2, 1):
            wavelength, ratio = ql.recalibrated(science[0]['path'], standard['path'], order)
            ratio = np.where(ql.outside_telluric(wavelength), ratio, np.nan)
            ax.plot(wavelength, ratio, lw=0.8, color=f'C{number}',
                    label=standard['date'][:10] if order == 2 else None)
    ax.axhline(1, color='k', lw=0.6)
    ax.set_xlim(3200, 10000)
    ax.set_xlabel('observed wavelength (Å)')
    ax.set_ylabel('one standard / ensemble')
    ax.set_title(f"{science[0]['object']} {science[0]['date'][:10]} calibrated with each standard alone", loc='left')
    ax.legend(title='standard taken on', fontsize=8)
    fig.tight_layout()

# %% [markdown]
# For reverberation mapping the overall level does not matter: the spectra are scaled to a common [O III] flux
# afterwards. What matters is the shape, which the choice of standard changes by much less.

# %% [markdown]
# ## Next steps
# * Every product carries quality keywords and flags; `banzai_floyds_qc` writes them for a whole folder to one
#   table (`docs/QC_KEYWORDS.md`).
# * The `-1d` files hold the spectra: use the rows of the `EXTRACTED` table with `mask == 0`, columns
#   `wavelength`, `flux_optimal`, `fluxerror_optimal`, `order`.
# * The report `floyds_rm_report.pdf` explains each step and how much it improves on BANZAI 1.0.1 and 1.1.0.
