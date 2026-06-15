# Elma: Automated Galaxy Bar Detection Pipeline

`elma` is a Python package for the automated measurement of galactic bars in FITS imaging data. Utilizing iterative elliptical-isophote fitting, the code measures the central galactic component on the assumption that this region is bar-dominated, which holds true for galaxies with prominent bars. It combines this geometric approach with manual cosmological calculations to convert raw telescope data into physical bar-length measurements in kiloparsecs.

## Features

* **Flexible input:** handles 2D single-band images and 3D data cubes (collapsed spectrally via sum); automatically finds the first populated HDU.
* **Isophote fitting:** fits concentric ellipses outward from the galaxy centre following the Jedrzejewski (1987) method as implemented in `photutils`.
* **Bar detection:** identifies the bar as the isophote at the peak ellipticity, excluding the PSF-dominated nucleus (sma < 3 px).
* **Physical units:** converts pixels to kpc using the WCS pixel scale and angular diameter distance integrated from a flat ΛCDM cosmology (H₀ = 70, Ωm = 0.3, ΩΛ = 0.7).
* **Adaptive zoom:** output plots are centred on the galaxy and zoomed to the bar scale, with a minimum field floor to avoid over-cropping compact galaxies.
* **High-quality visualisation:** 300 DPI RGB output using Lanczos interpolation and `make_lupton_rgb` stretch.

## Installation

```bash
git clone https://github.com/BrunaLimaa/elma.git
cd elma
pip install .
```

## Usage

### Single galaxy

```python
from elma import run_pipeline

bar_size_kpc = run_pipeline(filename="barred1.fits", redshift=0.42)
print(f"Detected bar length: {bar_size_kpc:.2f} kpc")
```

### Batch processing

A ready-to-use batch script is provided at [`example_batch_run.py`](example_batch_run.py). It processes a directory of FITS files, each with its own redshift, and prints a summary table of bar lengths and filter metadata.

```bash
python example_batch_run.py
```

The script iterates over a dictionary of `filename: redshift` pairs, calls `run_pipeline` for each, collects results, and reports them alongside any filter keywords found in the FITS headers.

After all files are processed, a summary table is printed to the terminal with the detected bar length and FITS filter metadata for each galaxy:

```
Summary of Results:
------------------------------------------------------------
galaxy_a.fits  :    10.96 kpc  |  filters: FILTER1=F090W, FILTER2=F115W, FILTER3=F150W
galaxy_b.fits  :    13.91 kpc  |  filters: FILTER1=F090W, FILTER2=F115W, FILTER3=F150W
galaxy_c.fits  :    10.28 kpc  |  filters: FILTER1=F090W, FILTER2=F150W, FILTER3=F182M
------------------------------------------------------------
```

## Outputs

For each processed galaxy `elma` saves three files next to the input FITS:

| File | Contents |
|---|---|
| `*_analysis_all.png` | Galaxy + all fitted isophotes (cyan, uniform weight) |
| `*_analysis_bar_only.png` | Galaxy + the detected bar ellipse only (red), with kpc label |
| `*_DEBUG_INPUT.png` | Raw RGB input check before any processing |

## Methodology

### Centre and initialisation
The galaxy centre is defined as the position of the brightest pixel. The isophote fitter is seeded at a semi-major axis of 5 pixels with ellipticity ε = 0.2, and the centre is held fixed throughout.

### Bar length
Because the pipeline relies strictly on geometric surface brightness rather than statistical structural decomposition, it assumes the central component being fit is the bar. Under this assumption, the bar radius is extracted as the semi-major axis of the valid isophote exhibiting the highest ellipticity. The full bar length is twice this radius.

### Physical scaling
The pixel scale is read from the WCS header. The angular diameter distance D_A is integrated numerically from the flat ΛCDM Friedmann equation. Bar length in kpc = 2 × radius_px × (pixel scale in rad) × D_A × 1000.

## Testing

```bash
python -m pytest tests/test_core.py
```

## Dependencies

* `numpy`
* `matplotlib`
* `astropy`
* `photutils`
* `scipy`
* `pytest`

## Author

**Bruna Lima**  
Computer Science Undergraduate, UFRGS
