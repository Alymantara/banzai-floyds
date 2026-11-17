# Installing and running the improved BANZAI-FLOYDS (rm-improvements branch)

This branch adds the improved reduction for reverberation mapping (improved v3 in the FLOYDS RM report): harmonic
fringe model fitted jointly with a telluric template, residual fringe template per slit, physical object-trace
model, band-limited ensemble telluric correction, ensemble flux calibration and the O2 A band flagged in science
spectra. All of it is on by default.

## 1. Get the code

```bash
git clone https://github.com/Alymantara/banzai-floyds.git
cd banzai-floyds
git checkout rm-improvements
```

## 2. A Python environment native to the machine

BANZAI needs PyTorch >= 2.3. On Apple-silicon Macs the environment **must be arm64**: an Intel (x86_64) Anaconda
running under Rosetta has no PyTorch >= 2.3 wheels and the installation fails. Check with

```bash
python -c "import platform; print(platform.machine())"     # must print arm64 on Apple silicon
```

Create a native Python 3.12 environment, for example with conda (works even from an Intel Anaconda):

```bash
CONDA_SUBDIR=osx-arm64 conda create -y -p ./env-banzai python=3.12 pip
```

(or Miniforge for arm64, or `uv venv --python 3.12`). On Linux any Python 3.10-3.12 virtual environment works.

## 3. Install with a clean compiler environment

`lco-banzai` builds C extensions. If the login profile exports compiler settings (`CC`, `CFLAGS`, `LDFLAGS`, ...)
from an Intel conda or oneAPI installation, the build fails with errors such as

    error: This is not the Intel oneAPI compiler
    ld: unsupported tapi file type '!tapi-tbd' ... for architecture x86_64

Run the installation (and, simplest, everything else) from a clean environment with the system compiler. A small
wrapper script avoids retyping it; in zsh do not put the command in a variable, since zsh does not split words:

```bash
cat > benv <<'EOF'
#!/bin/sh
# run a command in the BANZAI environment with a clean shell environment
exec env -i HOME="$HOME" PATH="$(cd "$(dirname "$0")" && pwd)/env-banzai/bin:/usr/bin:/bin:/usr/sbin:/sbin" \
  OPENTSDB_PYTHON_METRICS_TEST_MODE=True CC=/usr/bin/clang CXX=/usr/bin/clang++ \
  ${DB_ADDRESS:+DB_ADDRESS="$DB_ADDRESS"} "$@"
EOF
chmod +x benv
./benv python -m pip install -e ./banzai-floyds
```

(The `benv` file and `env-banzai/` sit next to the `banzai-floyds` checkout.) On Linux a plain
`pip install -e .` inside the virtual environment is usually enough.

## 4. Check the installation

```bash
cd banzai-floyds
../benv python -m pytest banzai_floyds/tests --deselect banzai_floyds/tests/test_e2e.py
```

Expected: 176 passed, 1 xfailed. The end-to-end tests (`test_e2e.py`) need LCO's test database and are skipped.

## 5. Reduce a set of raw frames offline

`OPENTSDB_PYTHON_METRICS_TEST_MODE=True` (set by `benv`) stops BANZAI from trying to send metrics.

```bash
# calibration database with both FLOYDS sites, order locations, LSF parameters, flux standards and the sky-flat
# order solutions (they ship with the package, one per frame of banzai_floyds/data/orders/skyflats.dat, and are
# registered from where the package is installed: nothing to download)
./benv banzai_floyds_offline setup-db --db-address sqlite:///$PWD/floyds.db

# the whole reduction: arcs -> lamp flats -> fringe phase masters -> standards -> residual fringe templates ->
# standards again -> science
./benv banzai_floyds_offline chain --db-address sqlite:///$PWD/floyds.db \
    --raw-dir raw/ --processed-dir processed/ --workers 4
```

`raw/` holds the raw frames of the science targets, the flux standards, and the arcs and lamp flats of their blocks.
Standards are recognised by position (`banzai_floyds/data/standards`). Single steps can be run with
`banzai_floyds_offline run --kind {arcs,flats,standards,science}`, and the fringe products with
`banzai_floyds_make_phase_master` and `banzai_floyds_make_residual_template`.

Products: `processed/<site>/<camera>/<date>/processed/*-1d.fits.fz` (extracted, flux-calibrated spectra; in the
EXTRACTED table use `mask == 0`; bit 2 marks the O2 A band in science spectra) and `*-2d.fits.fz`.

## 6. Practical notes

* **Disk.** A processed 2-d frame is ~50 MB (most of it the BINNED2D table) and BANZAI memory-maps arrays in
  the temporary directory; keep > 10 GB free.
* **Time.** On an 8-core laptop: arc ~34 s, lamp flat ~4 s, standard ~12 s (processed twice), science ~10 s per
  frame; a science visit (arc, lamp flat, two exposures) ~58 s and a standard visit ~62 s with one worker, ~15 s of
  wall-clock time per science visit with four workers. (The arc is ~48 s before the faster line fits of 327e941; the
  visit times are scaled from that measured 29 % saving.)
* **Sites.** The object-trace template and its priors exist for COJ and OGG (built from each site's standards with
  `work/scripts/build_trace_template.py`; without a template the trace falls back to the 1.1.0 polynomial). The
  residual fringe template is built per slit from the standards of that slit; without 6" standards there is none for
  the 6" slit.
* **Sky flats.** The order solutions of every sky flat in `data/orders/skyflats.dat` ship in
  `banzai_floyds/data/orders/skyflats/` (images zeroed, ~170 kB each) and `setup-db` registers them with their
  validity windows. Re-run `setup-db` if the package moves. To reprocess them, or after adding a new set to
  `skyflats.dat`: `banzai_floyds_offline fetch-skyflats --raw-dir skyflats/` downloads the listed raw frames from the
  LCO archive (public, one request a second), `banzai_floyds_offline skyflats ...` reduces them (it also takes the
  archive-processed `-f91` files), and `python scripts/pack_skyflats.py processed/` writes the shipped copies.
  For OGG the only listed set is 2019-03-29; it works for the 2024-25 data.
* **Flux standards.** Only stars in `banzai_floyds/data/standards` are used as standards (BD+28 4211, EG 274,
  Feige 34, Feige 67, Feige 110, GD 108, GD 153); other standards are reduced as science. Make sure every campaign has
  standards before and after any change of the telescope throughput.
* **Disk.** After a run, `work/banzai_run/lean_2d.py --drop-binned <processed>` drops the diagnostic images and the
  BINNED2D table from finished -2d spectra (examples kept): ~2 GB per site-season.
* **Stock behaviour.** Setting `FRINGE_METHOD='2d'`, `PROFILE_TRACE_MODEL='polynomial'`, `TELLURIC_MODEL_BANDS=None`,
  `TELLURIC_ENSEMBLE_DAYS=None`, `FLUX_ENSEMBLE_DAYS=None`, `TELLURIC_FLAG_BANDS=None` in `banzai_floyds/settings.py`
  restores 1.1.0.
