# Elma: Automated Galaxy Bar Detection Pipeline

`Elma` is a Python package designed to automate the detection and measurement of galactic bars in FITS imaging data. It combines photometic isophote fitting with cosmological calculations to transform raw telescope data into physical scientific measurements.

## Features

* **Adaptive Input:** Automatically handles 2D images and 3D IFU data cubes (via spectral collapse).
* **Isophote Fitting:** fits galaxy light profiles using Iterative Ellipse Fitting.
* **Automated Detection:** Identifies bar length via peak ellipticity, masking the noisy nucleus.
* **Physical Units:** Converts pixels to Kiloparsecs (kpc) using WCS and LambdaCDM Cosmology.
* **Visualization:** Exports diagnostic plots of the fitted isophotes and detected bar vectors.

## Installation

### Prerequisites
* Python 3.8+

### Install via pip
You can install the package and its dependencies directly from the source:

```bash
git clone https://github.com/BrunaLimaa/elma.git
cd elma
pip install .

```

Or install the dependencies manually:

```bash
pip install -r requirements.txt

```

## Usage

The package is designed to be imported and run with a single high-level command.

```python
from elma import run_pipeline

# Path to your FITS file (2D image or 3D cube)
fits_file = "data/galaxy_sample.fits"

# You MUST provide the redshift (z) for accurate physical sizing
galaxy_redshift = 0.42

# Run the pipeline
# Returns the physical size of the bar in kpc
bar_size_kpc = run_pipeline(filename=fits_file, redshift=galaxy_redshift)

print(f"Detected Bar Length: {bar_size_kpc:.2f} kpc")

```

## Outputs

For every processed galaxy, `elma` generates two artifacts in the working directory:

* **`*_analysis.png`**: The primary result.
* **Background:** Log-scaled flux map of the galaxy.
* **Blue Ellipses:** The isophote fits tracing the galaxy's potential.
* **Red Ellipse:** The specific isophote identified as the "Bar" based on the peak ellipticity.


* **`*_DEBUG_INPUT.png`**: A quality control plot showing the raw data (or the collapsed cube) before processing, useful for verifying signal-to-noise ratios.

## Methodology

### 1. The Algorithm

The pipeline treats galaxy morphology as a signal processing problem. It fits concentric ellipses (isophotes) expanding from the galaxy center.

* **Step Size:** 0.05 (Fixed high-resolution sampling).
* **Safety Mask:** The central region (pixels) is ignored to prevent the Point Spread Function (PSF) or nuclear noise from triggering false positives.

### 2. The Physics

To enable fair comparisons across different redshifts, `elma` calculates the size at the epoch of emission.

## 🛠 Dependencies

* `numpy`
* `matplotlib`
* `astropy` (FITS handling, WCS, Cosmology)
* `photutils` (Isophote fitting engine)
* `scipy`

## 👤 Author

**Bruna Lima**

Computer Science Undergraduate, UFRGS


