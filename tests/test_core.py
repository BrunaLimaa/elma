import pytest
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from elma.core import calculate_kpc, load_and_process

def test_calculate_kpc():
    # Setup a mock WCS with 1 arcsec per pixel (1/3600 deg)
    w = WCS(naxis=2)
    w.wcs.cdelt = [1/3600, 1/3600]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    
    # Test parameters
    bar_length_pixels = 100
    redshift = 0.1
    
    size_kpc = calculate_kpc(bar_length_pixels, redshift, w)
    
    assert size_kpc > 0
    assert isinstance(size_kpc, float)
    # Precise check based on astropy's calculation
    assert pytest.approx(size_kpc, rel=0.01) == 184.4  # Value at z=0.1 with standard cosmo

def test_load_and_process_2d(tmp_path):
    # Create a dummy 2D FITS file
    data = np.ones((50, 50))
    header = fits.Header()
    header['CTYPE1'] = 'RA---TAN'
    header['CTYPE2'] = 'DEC--TAN'
    fits_file = tmp_path / "test_2d.fits"
    fits.writeto(fits_file, data, header)
    
    img, wcs_info = load_and_process(str(fits_file))
    
    assert img.shape == (50, 50)
    assert isinstance(wcs_info, WCS)

def test_load_and_process_3d(tmp_path):
    # Create a dummy 3D FITS file
    data = np.ones((10, 50, 50))
    header = fits.Header()
    fits_file = tmp_path / "test_3d.fits"
    fits.writeto(fits_file, data, header)
    
    img, wcs_info = load_and_process(str(fits_file))
    
    # Should be collapsed to 2D
    assert img.shape == (50, 50)
    assert np.all(img == 10)  # Sum of 10 layers of 1s

def test_find_bar_geometry():
    # Create a synthetic 2D galaxy (Gaussian profile)
    y, x = np.mgrid[0:100, 0:100]
    # An elongated Gaussian to simulate a bar
    z = np.exp(-(((x-50)**2 / (2*20**2)) + ((y-50)**2 / (2*10**2))))
    
    from elma.core import find_bar_geometry
    isolist, radius_px = find_bar_geometry(z)
    
    assert len(isolist) > 0
    assert radius_px > 0
    assert radius_px < 50


