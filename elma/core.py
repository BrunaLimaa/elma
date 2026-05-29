import os
import logging
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.visualization import simple_norm, make_lupton_rgb
from astropy.wcs import WCS
from astropy import wcs
import astropy.units as u
from photutils.isophote import EllipseGeometry, Ellipse, IsophoteList
from matplotlib.patches import Ellipse as MplEllipse
from typing import Tuple, Optional, Any

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

_FILTER_KEYS = ("FILTER", "FILTER1", "FILTER2", "FILTER3", "FILTNAM1", "FILTNAM2", "FILTNAM3", "PHOTFILT", "BAND")

def get_fits_filters(filename: str) -> dict:
    """
    Returns filter/band header keywords found across all HDUs in a FITS file.
    Keys are header keyword names, values are the corresponding header values.
    """
    found = {}
    with fits.open(filename) as hdul:
        for hdu in hdul:
            for key in _FILTER_KEYS:
                if key in hdu.header and key not in found:
                    found[key] = hdu.header[key]
    return found


def calculate_angular_diameter_distance(redshift: float) -> float:
    """
    Calculates the angular diameter distance locally.
    """
    H0 = 70 # hubble constant
    c = 299792.458 # speed of light
    D_H = c / H0

    omega_m = 0.30
    omega_l = 0.70

    def E(z_prime):
        return np.sqrt(omega_m * (1 + z_prime)**3 + omega_l)

    num_steps = 10000
    dz = redshift / num_steps
    
    integral_sum = 0.0
    current_z = 0.0
    
    for _ in range(num_steps):
        left_height = 1.0 / E(current_z)
        right_height = 1.0 / E(current_z + dz)
        slice_area = ((left_height + right_height) / 2.0 ) * dz
        integral_sum += slice_area
        current_z += dz
        
    D_C = D_H * integral_sum
    D_M = D_C
    D_A = D_M / (1 + redshift)
    
    return D_A

def load_image_rgb(filename: str) -> np.ndarray:
    """
    Loads 3D FITS data and converts to an RGB image using make_lupton_rgb.
    """
    with fits.open(filename) as hdul:
        data = next(hdu.data for hdu in hdul if hdu.data is not None)

        if data.ndim == 3:
            r, g, b = data[0], data[1], data[2]
        else:
            r = g = b = data

        r = r / np.nanpercentile(r, 99)
        g = g / np.nanpercentile(g, 99)
        b = b / np.nanpercentile(b, 99)

    rgb_img = make_lupton_rgb(r, g, b, stretch=1.5, Q=5, minimum=0.01)
    
    return rgb_img

def load_and_process(filename: str) -> Tuple[np.ndarray, WCS]:
    """
    Processes the FITS file. Extracts the image and the header info.
    
    Finds the first HDU that contains data. If the data is 3D, it collapses it to 2D.

    Args:
        filename: Path to the FITS file.

    Returns:
        A tuple containing the 2D image data (float) and the WCS information.

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If no valid image data is found in the HDU list.
    """
    with fits.open(filename) as hdul:
        # Find the first HDU with data
        data_hdu = None
        for hdu in hdul:
            if hdu.data is not None:
                data_hdu = hdu
                break
        
        if data_hdu is None:
            raise ValueError(f"No data found in FITS file: {filename}")

        header = data_hdu.header
        wcs_info = WCS(header)
        data = data_hdu.data
        
        if data.ndim == 3:
            logger.info("[%s] Detected 3D cube. Collapsing to 2D via sum.", filename)
            image_2d = np.sum(data, axis=0)
        elif data.ndim == 2:
            logger.info("[%s] Detected 2D Image.", filename)
            image_2d = data
        else:
            raise ValueError(f"Unsupported data dimensions: {data.ndim} in {filename}")
        
        return image_2d.astype(float), wcs_info

def calculate_kpc(
    bar_length_pixels: float, 
    redshift: float, 
    wcs_info: WCS, 
    h0: float = 70.0, 
    om0: float = 0.3
) -> float:
    """
    Converts pixel length to physical distance in kpc.
    """
    scales = wcs.utils.proj_plane_pixel_scales(wcs_info)
    deg_per_pixel = scales[0]
    
    bar_size_deg = bar_length_pixels * deg_per_pixel
    bar_size_rad = np.deg2rad(bar_size_deg)
    
    distance_mpc = calculate_angular_diameter_distance(redshift)
    distance_kpc = distance_mpc * 1000.0
    
    return bar_size_rad * distance_kpc

def find_bar_geometry(image_data: np.ndarray) -> Tuple[IsophoteList, float]:
    """
    Detects the bar geometry using isophote fitting.

    Args:
        image_data: 2D galaxy image.

    Returns:
        A tuple containing the list of fitted isophotes and the detected bar radius in pixels.
    """
    y_center, x_center = np.unravel_index(np.argmax(image_data), image_data.shape)
    
    g = EllipseGeometry(x0=x_center, y0=y_center, sma=5.0, eps=0.2, pa=0.0)
    ellipse = Ellipse(image_data, geometry=g)
    
    isolist = ellipse.fit_image(maxsma=450, fix_center=True, step=0.05)
    
    # Filter for valid fits (avoiding central pixels)
    valid_fits = isolist.sma > 3
    
    if np.any(valid_fits):
        idx_peak = np.argmax(isolist.eps[valid_fits])
        bar_radius_pixels = isolist.sma[valid_fits][idx_peak]
    else:
        bar_radius_pixels = 0.0
        logger.warning("No valid isophotes found for bar detection.")
        
    return isolist, bar_radius_pixels

def generate_plot(
    image_data: np.ndarray, 
    isolist: IsophoteList, 
    bar_radius_pixels: float, 
    kpc_size: float, 
    filename_out: str,
    wcs_info: Optional[WCS] = None
) -> None:
    """
    Generates a diagnostic plot and saves it to a file.
    """
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.imshow(image_data, origin='lower', interpolation='lanczos')

    # Determine galaxy center from the bar isophote (or image peak as fallback)
    cx, cy = None, None
    for iso in isolist:
        if np.isclose(float(iso.sma), bar_radius_pixels, atol=0.5):
            cx, cy = float(iso.x0), float(iso.y0)
            break
    if cx is None:
        h, w = image_data.shape[:2]
        cy, cx = h / 2, w / 2

    h, w = image_data.shape[:2]
    pad = max(bar_radius_pixels * 2, min(h, w) * 0.25)
    ax.set_xlim(cx - pad, cx + pad)
    ax.set_ylim(cy - pad, cy + pad)

    for iso in isolist:
        sma = float(iso.sma)
        if sma < 3: continue

        x0 = float(iso.x0)
        y0 = float(iso.y0)
        pa = float(iso.pa)
        eps = float(iso.eps)
        smi = sma * (1 - eps)

        e = MplEllipse(xy=(x0, y0),
                       width=2*sma,
                       height=2*smi,
                       angle=np.degrees(pa),
                       edgecolor='cyan',
                       facecolor='none',
                       linewidth=1.0,
                       alpha=0.3)
        ax.add_patch(e)

    length_bar = 2 * bar_radius_pixels


    plt.title(f"Analysis: {os.path.basename(filename_out)}")
    plt.savefig(filename_out, dpi=300, bbox_inches='tight')
    plt.close(fig)
    logger.info("Plot saved to: %s", filename_out)

def generate_bar_plot(
    image_data: np.ndarray, 
    isolist: IsophoteList, 
    bar_radius_pixels: float, 
    kpc_size: float, 
    filename_out: str,
    wcs_info: Optional[WCS] = None
) -> None:
    """
    Creates a focused plot with ONLY the galaxy + the chosen red ellipse + size label.
    """
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.imshow(image_data, origin='lower', interpolation='lanczos')

    cx, cy = None, None
    for iso in isolist:
        sma = float(iso.sma)
        if np.isclose(sma, bar_radius_pixels, atol=0.5):
            cx, cy = float(iso.x0), float(iso.y0)
            x0 = cx
            y0 = cy
            pa = float(iso.pa)
            eps = float(iso.eps)
            smi = sma * (1 - eps)

            e = MplEllipse(xy=(x0, y0),
                           width=2*sma,
                           height=2*smi,
                           angle=np.degrees(pa),
                           edgecolor='red',
                           facecolor='none',
                           linewidth=3.0)
            ax.add_patch(e)
            break

    if cx is None:
        h, w = image_data.shape[:2]
        cy, cx = h / 2, w / 2

    h, w = image_data.shape[:2]
    pad = max(bar_radius_pixels * 2, min(h, w) * 0.25)
    ax.set_xlim(cx - pad, cx + pad)
    ax.set_ylim(cy - pad, cy + pad)

    length_bar = 2 * bar_radius_pixels
    info_text = f"Bar Length: {kpc_size:.2f} kpc"
    ax.text(0.05, 0.95, info_text, transform=ax.transAxes, color='white',
            fontsize=24, va='top', ha='left',
            bbox=dict(facecolor='black', alpha=0.8, edgecolor='none', boxstyle='square,pad=0.8'))

    plt.title(f"Analysis (Bar Only): {os.path.basename(filename_out)}")
    plt.savefig(filename_out, dpi=300, bbox_inches='tight', facecolor='black')
    plt.close(fig)
    logger.info("Bar plot saved to: %s", filename_out)

def run_pipeline(filename: str, redshift: float) -> float:
    """
    Runs the full bar detection pipeline.
    """
    img_2d, wcs_info = load_and_process(filename)
    isolist, radius_px = find_bar_geometry(img_2d)
    img_rgb = load_image_rgb(filename)
    
    save_debug_image(img_rgb, filename)
    
    length_px = radius_px * 2
    size_kpc = calculate_kpc(length_px, redshift, wcs_info)

    output_name_all = filename.replace(".fits", "_analysis_all.png")
    output_name_bar = filename.replace(".fits", "_analysis_bar_only.png")
    
    generate_plot(img_rgb, isolist, radius_px, size_kpc, output_name_all)
    generate_bar_plot(img_rgb, isolist, radius_px, size_kpc, output_name_bar)
    
    return size_kpc

def save_debug_image(image_data: np.ndarray, filename_original: str) -> None:
    """
    Saves a debug image of the input data.
    """
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.imshow(image_data, origin='lower', interpolation='lanczos')
    ax.set_title("RGB INPUT CHECK:", fontsize=10)
    
    debug_name = filename_original.replace(".fits", "_DEBUG_INPUT.png")
    plt.savefig(debug_name, dpi=300, bbox_inches='tight')
    plt.close()
    logger.info("Debug image saved to: %s", debug_name)
