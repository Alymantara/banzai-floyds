"""A failed arc fit (flat red order, negative dispersion) must be caught at once, never hang the sky fit (AVA B100)."""
import time

import numpy as np
import pytest
from astropy.table import Table
from numpy.polynomial.legendre import Legendre
from scipy.interpolate import BSpline

from banzai_floyds.background import adaptive_knots, fit_background, uniform_knots
from banzai_floyds.utils.wavelength_utils import wavelength_solution_problems


def good_solution():
    return [Legendre.fit(np.arange(0, 1500), 4800 + 3.485 * np.arange(0, 1500), 1, domain=[0, 1500]),
            Legendre.fit(np.arange(530, 2000), 3300 + 1.723 * np.arange(0, 1470), 1, domain=[530, 2000])]


def test_good_solution_has_no_problems():
    assert wavelength_solution_problems(good_solution()) == []


@pytest.mark.parametrize('start, stop', [(2662.0, 2618.0), (16696.0, 16693.0)])
def test_flat_decreasing_red_order_is_caught(start, stop):
    # coj en12 20260721-0019 and 20250119-0009
    solution = good_solution()
    solution[0] = Legendre([0.5 * (start + stop), 0.5 * (stop - start)], domain=[0, 1500])
    problems = wavelength_solution_problems(solution)
    assert len(problems) == 1 and problems[0].startswith('order 1')


def test_wrong_scale_dispersion_is_caught():
    solution = good_solution()
    solution[1] = Legendre([4500.0, 0.1 * 735], domain=[530, 2000])   # 0.1 A/pixel instead of 1.7
    assert wavelength_solution_problems(solution)[0].startswith('order 2')


def test_adaptive_knots_refuses_negative_dispersion():
    wavelength = np.linspace(2615.5, 2662.4, 1500)
    knots = uniform_knots(wavelength, 1.0, 3)
    sky = BSpline(knots, np.ones(len(knots) - 4), 3)
    with pytest.raises(ValueError):
        adaptive_knots(wavelength, sky, -0.0303, 3)


def test_fit_background_skips_order_with_broken_wavelengths_quickly():
    rng = np.random.default_rng(1)
    rows = []
    for order, wavelengths in [(1, np.linspace(2662.4, 2615.5, 300)), (2, np.linspace(3300, 5700, 300))]:
        bin_width = np.median(np.diff(wavelengths))
        for y_order in np.arange(-20, 21):
            rows.append(Table({'x': np.arange(300), 'y': np.full(300, y_order + 25 + 50 * order),
                               'order': np.full(300, order), 'y_order': np.full(300, float(y_order)),
                               'wavelength': wavelengths, 'order_wavelength_bin': np.arange(1, 301),
                               'order_wavelength_bin_width': np.full(300, bin_width),
                               'y_profile': np.full(300, float(y_order)), 'profile_sigma': np.full(300, 2.0),
                               'data': 100.0 + rng.normal(0, 10, 300), 'uncertainty': np.full(300, 10.0),
                               'mask': np.zeros(300, dtype=int), 'extraction_window': np.full(300, abs(y_order) < 4)}))
    from astropy.table import vstack
    data = vstack(rows)
    start = time.time()
    background, fits = fit_background(data, read_noise=5.0)
    assert time.time() - start < 30
    assert 1 not in fits and 2 in fits
    assert np.all(np.isin(background['y'], data['y'][data['order'] == 2]))
