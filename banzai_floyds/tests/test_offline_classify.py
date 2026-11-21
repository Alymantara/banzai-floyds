"""AVA B50: raw frames whose compressed image extension is not named SCI must still be classified."""
import numpy as np
from astropy.io import fits

from banzai_floyds import offline


def _raw(path, obstype, extname):
    header = fits.Header({'OBSTYPE': obstype, 'DATE-OBS': '2014-08-04T06:00:00'})
    hdu = fits.CompImageHDU(np.zeros((8, 8), np.float32), header=header, name=extname)
    fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(path)


def test_raw_header_without_sci_extension(tmp_path):
    new, old = tmp_path / 'new-a00.fits.fz', tmp_path / 'old-a00.fits.fz'
    _raw(new, 'ARC', 'SCI')
    _raw(old, 'ARC', 'COMPRESSED_IMAGE')
    assert offline.raw_header(new)['OBSTYPE'] == 'ARC'
    assert offline.raw_header(old)['OBSTYPE'] == 'ARC'          # was '' (primary header), so the frame was skipped


def test_classify_counts_old_style_frames(tmp_path):
    _raw(tmp_path / 'x-0001-a00.fits.fz', 'ARC', 'COMPRESSED_IMAGE')
    _raw(tmp_path / 'x-0002-w00.fits.fz', 'LAMPFLAT', 'COMPRESSED_IMAGE')
    groups = offline.classify(str(tmp_path), 'sqlite:///unused.db')   # no SPECTRUM frame: no DB lookup
    assert len(groups['arcs']) == 1 and len(groups['flats']) == 1


def test_floyds_default_saturate_for_old_headers():
    """AVA B50: 2014-2016 en05/en06 raw headers have SATURATE = 0; banzai raised MissingSaturate (no FLOYDS default)."""
    from types import SimpleNamespace
    from banzai_floyds.frames import FLOYDSFrameFactory, FLOYDS_DEFAULT_SATURATE

    old = SimpleNamespace(meta=fits.Header({'SATURATE': 0.0}), binning=[1, 1], gain=2.0)
    image = SimpleNamespace(ccd_hdus=[old], meta=fits.Header({'SATURATE': 0.0}),
                            instrument=SimpleNamespace(type='2m0-FLOYDS-SciCam'))
    FLOYDSFrameFactory._init_saturate(image)
    assert old.meta['SATURATE'] == FLOYDS_DEFAULT_SATURATE and old.meta['MAXLIN'] == 56000.0

    new = SimpleNamespace(meta=fits.Header({'SATURATE': 38400.0}), binning=[1, 1], gain=2.3)
    image = SimpleNamespace(ccd_hdus=[new], meta=fits.Header({'SATURATE': 38400.0}),
                            instrument=SimpleNamespace(type='2m0-FLOYDS-SciCam'))
    FLOYDSFrameFactory._init_saturate(image)
    assert new.meta['SATURATE'] == 38400.0                       # a header value is never replaced here


def test_unnamed_image_extension_becomes_sci():
    """AVA B50: a 2014 raw image extension without EXTNAME is named SCI, so the product can be written."""
    from types import SimpleNamespace
    from banzai_floyds.frames import FLOYDSFrameFactory

    old = SimpleNamespace(name=None, meta=fits.Header())
    raw2016 = SimpleNamespace(name='RAW', meta=fits.Header({'EXTNAME': 'RAW'}))
    new = SimpleNamespace(name='SCI', meta=fits.Header({'EXTNAME': 'SCI'}))
    FLOYDSFrameFactory._name_unnamed_ccd_hdus(SimpleNamespace(ccd_hdus=[old, raw2016, new]))
    assert old.name == 'SCI' and old.meta['EXTNAME'] == 'SCI' and new.name == 'SCI'
    assert raw2016.name == 'SCI' and raw2016.meta['EXTNAME'] == 'SCI'          # 2016 raw files (AVA B60)
