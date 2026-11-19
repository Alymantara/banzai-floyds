# Quick start: reduce FLOYDS raw data and check the result

This guide takes you from a folder of raw FLOYDS frames to extracted, flux-calibrated spectra and a short set of
quality checks. Install the pipeline first ([INSTALL_RM.md](INSTALL_RM.md), sections 1-4).

## What you need

One folder of raw frames (`*-e00.fits.fz`, `*-a00.fits.fz`, `*-w00.fits.fz`, as downloaded from the LCO archive)
with, for every block (one telescope request):

- the `ARC` and the `LAMPFLAT` of the block, and
- its `SPECTRUM` exposures,

plus at least one block of a flux standard (Feige 110, Feige 34, Feige 67, GD 153, BD+28 4211, ...) at the same
site within the ~60 days around the science. Sky flats are not needed: their order solutions ship with the package.

## Option 1: the notebook

`notebooks/quickstart_reduction.ipynb` reduces the folder and shows:

1. the raw frames, grouped by night, type and slit;
2. the reduction, with the time of each step;
3. the extracted spectra of every object;
4. the accuracy of the wavelength solution (residuals of the arc lines);
5. the accuracy of the flux calibration (each standard calibrated with another standard of the same slit);
6. the signal-to-noise of every spectrum;
7. how the standard stars enter each science spectrum, and how much its flux changes with the standard used.

Jupyter is not installed with the pipeline. In the pipeline's environment:

```bash
pip install jupyterlab
```

Then point the notebook to your data, either by editing `RAW_DIR` and `WORK_DIR` in its first code cell or with
environment variables:

```bash
FLOYDS_RAW=/path/to/raw FLOYDS_WORK=/path/to/work jupyter lab notebooks/quickstart_reduction.ipynb
```

`WORK_DIR` receives the calibration database (`floyds.db`), the products (`processed/`) and the timings
(`timings.json`). Run all the cells. The reduction reduces every frame each time it is run; set `REDUCE = False` to
look at products already reduced.

The notebook is generated from `notebooks/quickstart_reduction.py` (same cells, plain Python) by
`python notebooks/make_notebook.py`; edit the `.py` and regenerate.

## Option 2: the command line

```bash
banzai_floyds_offline setup-db --db-address sqlite:///work/floyds.db
```

```bash
banzai_floyds_offline chain --db-address sqlite:///work/floyds.db --raw-dir raw/ --processed-dir work/processed/ --workers 4
```

```bash
banzai_floyds_qc work/processed/ --table work/qc.csv
```

`setup-db` is needed once per database. `chain` reduces, in order: arcs, lamp flats, fringe phase masters,
standards (first pass), residual fringe templates, standards (second pass) and science. `banzai_floyds_qc` collects the
quality keywords of every product into one table ([QC_KEYWORDS.md](QC_KEYWORDS.md)).

From Python, the same functions and the checks of the notebook:

```python
from banzai_floyds import offline
from banzai_floyds.diagnostics import quicklook as ql

if __name__ == '__main__':          # needed in a script (not in a notebook): the workers are separate processes
    offline.setup_db('sqlite:///work/floyds.db')
    timings = {}
    offline.chain('raw', 'sqlite:///work/floyds.db', 'work/processed', workers=4, timings=timings)
    print(ql.spectra('work/processed'))                   # one row per extracted spectrum
    print(ql.arc_accuracy('work/processed'))              # arc-line residuals per arc and order
    print(ql.standard_flux_accuracy('work/processed'))    # standard calibrated with another standard
    print(ql.science_vs_standards('work/processed'))      # science calibrated with each standard
```

## The products

For every science or standard exposure, in `processed/<site>/<instrument>/<night>/processed/`:

- `...-e91-1d.fits.fz`: the extracted spectrum. Extension `EXTRACTED` has, per order (1 red, 2 blue),
  `wavelength`, `flux_optimal` and `fluxerror_optimal` (erg s⁻¹ cm⁻² Å⁻¹) and `mask` (non-zero: do not use);
  `SENSITIVITY` and `TELLURIC` the calibrations applied. The header has the quality keywords (`L1Q...`), the
  standards used (`L1FLXNST`, `L1TELNST`) and the arc (`L1IDARC`).
- `...-e91-2d.fits.fz`: the calibrated 2-d image.
- `...-a91.fits.fz`: the reduced arc, with its line residuals (`RESIDUALS`).

## The example

The checks below were run on 21 raw frames from COJ, 2026 September 12-14: 6 science exposures of Fairall 9 (6″
slit), 3 of the standard Feige 110 (2″ slit), and their 6 arcs and 6 lamp flats. With 4 workers on an Apple M-series
laptop:

| step | frames | wall-clock (s) | one frame, one worker (s, median) |
|---|---|---|---|
| set up the database | - | 4 | - |
| arcs | 6 | 75 | 39 |
| lamp flats | 6 | 11 | 6.5 |
| fringe phase masters | - | 1 | - |
| standards, first pass | 3 | 14 | 12 |
| residual fringe templates | - | 1 | - |
| standards, second pass | 3 | 12 | 10 |
| science | 6 | 48 | 31 |
| **total** | 21 | **165 (2.7 min)** | |

About 27 s of wall-clock time per block (arc, lamp flat and its spectra). Arcs take longest: every line is fitted in
every detector row along the slit to measure the tilt and the shape of the lines.

Typical results of the checks on this set:

| check | red order | blue order |
|---|---|---|
| arc residuals, robust rms (median of 6 arcs) | 0.14 Å, 16-24 lines | 0.53 Å, 5-6 lines |
| arc lines rejected as misidentified | Ar 6753 in two of the three 2″ arcs | none |
| flux of a standard calibrated with another standard: shape error (rms) | 1.2 % | 1.8 % |
| the same: overall level (clouds, slit losses; 2″ slit in 2.4-2.8″ seeing) | −9 to +10 % | −9 to +10 % |
| science calibrated with one standard instead of the ensemble (here all three): level | −3.9 to +5.9 % | −4.2 to +5.0 % |
| the same: shape (rms; max) | 0.4 %; 0.6 % | 0.5 %; 1.5 % |
| S/N per bin, Fairall 9, 600 s | 87-94 | 59-68 |

The level of a science spectrum follows the grey level of its standards; in reverberation mapping it is rescaled to
the narrow [O III] flux anyway. The shape, set by the ensemble, changes by less than 1 %.
