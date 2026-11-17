"""Write the sky-flat order solutions shipped with the package (banzai_floyds/data/orders/skyflats/).

Only the ORDERS map and the ORDER_COEFFS table of a processed sky flat (-f91) are used by the pipeline
(orders.OrderLoader); the images are set to zero so that each file compresses to ~170 kB, with every header kept.
Run it on the processed sky flats of all the frames in data/orders/skyflats.dat, e.g. after
    banzai_floyds_offline fetch-skyflats --raw-dir skyflats/
    banzai_floyds_offline skyflats --db-address sqlite:///sf.db --raw-dir skyflats/ --processed-dir processed/

usage: python scripts/pack_skyflats.py PROCESSED_DIR"""
import glob
import importlib.resources
import os
import sys

import numpy as np
from astropy.io import ascii, fits


def pack(processed_dir: str) -> int:
    data = importlib.resources.files('banzai_floyds') / 'data' / 'orders'
    listed = ascii.read(str(data / 'skyflats.dat'))
    out_dir = str(data / 'skyflats')
    os.makedirs(out_dir, exist_ok=True)
    found = {os.path.basename(f): f for f in glob.glob(os.path.join(processed_dir, '**', '*-f91.fits.fz'),
                                                        recursive=True)}
    written = 0
    for name in listed['filename']:
        name = name.replace('-x00', '-f91')
        if name not in found:
            print('missing', name)
            continue
        with fits.open(found[name]) as hdul:
            for extension in ('SCI', 'ERR', 'BPM'):
                hdul[extension].data = np.zeros_like(hdul[extension].data)
            hdul.writeto(os.path.join(out_dir, name), overwrite=True)
        written += 1
    return written


if __name__ == '__main__':
    print('written', pack(sys.argv[1]))
