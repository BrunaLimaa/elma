# Elma: Automated Galaxy Bar Detection Pipeline

`Elma` is a Python package designed to automate the detection and measurement of galactic bars in FITS imaging data. It combines photometric isophote fitting with manual cosmological calculations to transform raw telescope data into physical scientific measurements.

## Features

* **Adaptive Input:** Automatically handles 2D images and 3D IFU data cubes (via spectral collapse).
* **Isophote Fitting:** Fits galaxy light profiles using Iterative Ellipse Fitting.
* **Automated Detection:** Identifies bar length via peak ellipticity, masking the noisy nucleus.
* **Physical Units:** Converts pixels to Kiloparsecs (kpc) using WCS and manual angular diameter distance integration.
* **High-Quality Visualization:** Exports high-resolution (300 DPI) RGB diagnostic plots with Lanczos interpolation for maximum clarity.

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

The package is designed to be run with a single high-level command.

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

### Batch Processing
An example script is provided to process multiple files with specific redshifts:

```bash
python example_batch_run.py
```

## Outputs

For every processed galaxy, `elma` generates three artifacts:

* **`*_analysis_all.png`**: Full-field RGB image showing all fitted isophotes (Cyan) and the detected Bar (Red).
* **`*_analysis_bar_only.png`**: Full-field RGB image highlighting only the detected Bar ellipse.
* **`*_DEBUG_INPUT.png`**: A quality control plot showing the RGB input before processing.

All plots are generated at **300 DPI** using **Lanczos interpolation** to match the quality of the original input data.

## Methodology

### 1. The Algorithm
The pipeline fits concentric ellipses expanding from the galaxy center using the `photutils` isophote engine. It identifies the bar by locating the peak in ellipticity while ignoring the PSF-dominated central pixels.

### 2. The Physics
To enable fair comparisons across different redshifts, `elma` calculates the angular diameter distance using manual trapezoidal integration of the LambdaCDM expansion history ($H_0=70, \Omega_m=0.3, \Omega_\Lambda=0.7$).

## Testing
The project includes a suite of unit tests using `pytest`.

```bash
python -m pytest tests/test_core.py
```

## 🛠 Dependencies

* `numpy`
* `matplotlib`
* `astropy` (FITS, WCS)
* `photutils` (Isophote fitting)
* `scipy`
* `pytest`

## 👤 Author

**Bruna Lima**
Computer Science Undergraduate, UFRGS
