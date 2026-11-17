import importlib.resources
import os

from astropy.io import ascii, fits


def test_every_listed_sky_flat_ships_with_its_order_solution():
    # setup-db registers the order solutions from the package, so every frame of skyflats.dat must have one, with the
    # extensions OrderLoader reads, for the site, slit and date of the list
    data = importlib.resources.files('banzai_floyds') / 'data' / 'orders'
    listed = ascii.read(str(data / 'skyflats.dat'))
    for row in listed:
        path = str(data / 'skyflats' / row['filename'].replace('-x00', '-f91'))
        assert os.path.exists(path), path
        with fits.open(path) as hdul:
            assert 'ORDERS' in hdul and 'ORDER_COEFFS' in hdul
            assert len(hdul['ORDER_COEFFS'].data) == 2
            header = hdul['SCI'].header
            assert header['SITEID'].strip() == row['site']
            assert float(header['APERWID']) == float(row['slit_width'])
