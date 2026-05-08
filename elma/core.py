import logging
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.visualization import simple_norm
from astropy.wcs import WCS
from astropy import wcs
from astropy.cosmology import LambdaCDM
import astropy.units as u
from photutils.isophote import EllipseGeometry, Ellipse, IsophoteList
from matplotlib.patches import Ellipse as MplEllipse
from typing import Tuple, Optional, Any

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

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

    Args:
        bar_length_pixels: Length of the bar in pixels.
        redshift: Redshift of the galaxy.
        wcs_info: WCS information from the FITS header.
        h0: Hubble constant (km/s/Mpc). Defaults to 70.0.
        om0: Omega matter. Defaults to 0.3.

    Returns:
        Physical size in kpc.
    """
    cosmo = LambdaCDM(H0=h0, Om0=om0, Ode0=1.0 - om0)
    
    scales = wcs.utils.proj_plane_pixel_scales(wcs_info)
    # Use average scale if non-square pixels
    deg_per_pixel = np.mean(scales)
    
    bar_size_deg = bar_length_pixels * deg_per_pixel
    bar_size_rad = np.deg2rad(bar_size_deg)
    
    distance_mpc = cosmo.angular_diameter_distance(redshift)
    distance_kpc = distance_mpc.to(u.kpc).value
    
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
    filename_out: str
) -> None:
    """
    Generates a diagnostic plot and saves it to a file.

    Args:
        image_data: 2D galaxy image.
        isolist: List of fitted isophotes.
        bar_radius_pixels: Detected bar radius in pixels.
        kpc_size: Physical size of the bar in kpc.
        filename_out: Output path for the plot.
    """
    good_mask = ~np.logical_or(np.isnan(image_data), image_data <= 0)
    if not np.any(good_mask):
        logger.error("Cannot generate plot: all data is invalid (NaN or <= 0).")
        return

    # Using a higher percentile for better dynamic range in the log stretch
    norm = simple_norm(image_data[good_mask], stretch='log', percent=99.5)
    
    fig, ax = plt.subplots(figsize=(12, 12))
    # 'magma' is a perceptually uniform colormap, excellent for astrophysical data
    im = ax.imshow(image_data, cmap='magma', origin='lower', norm=norm)
    
    # Add a refined colorbar
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('Flux (Log Scale)', fontsize=12)

    for iso in isolist:
        if iso.sma < 3: 
            continue 

        is_bar = np.isclose(iso.sma, bar_radius_pixels, atol=0.5)
        
        # Cyan isophotes and Yellow bar provide high contrast against 'magma' background
        color = 'yellow' if is_bar else 'cyan'
        alpha = 1.0 if is_bar else 0.4
        linewidth = 3.5 if is_bar else 1.2
        
        angle_deg = np.degrees(iso.pa)
        smi = iso.sma * (1 - iso.eps)
        
        e = MplEllipse(xy=(iso.x0, iso.y0),
                       width=2*iso.sma,
                       height=2*smi,
                       angle=angle_deg,
                       edgecolor=color,
                       facecolor='none',
                       linewidth=linewidth,
                       alpha=alpha)
        ax.add_patch(e)
    
    length_px = bar_radius_pixels * 2
    info_text = f"Bar Length: {kpc_size:.2f} kpc\n({length_px:.1f} pixels)"
    ax.text(0.02, 0.98, info_text, transform=ax.transAxes, color='yellow', 
            fontsize=14, verticalalignment='top', fontweight='bold',
            bbox=dict(facecolor='black', alpha=0.6, edgecolor='yellow'))
    
    plt.title(f"Galaxy Bar Analysis: {os.path.basename(filename_out)}", fontsize=16, pad=20)
    
    plt.savefig(filename_out, dpi=200, bbox_inches='tight')
    plt.close() 
    logger.info("Enhanced plot saved to: %s", filename_out)

def run_pipeline(
    filename: str, 
    redshift: float, 
    h0: float = 70.0, 
    om0: float = 0.3
) -> float:
    """
    Runs the full bar detection pipeline.

    Args:
        filename: Path to the FITS file.
        redshift: Galaxy redshift.
        h0: Hubble constant. Defaults to 70.0.
        om0: Omega matter. Defaults to 0.3.

    Returns:
        The detected bar size in kpc.
    """
    img, wcs_info = load_and_process(filename)
    
    isolist, radius_px = find_bar_geometry(img)
    
    save_debug_image(img, filename)
    
    length_px = radius_px * 2
    size_kpc = calculate_kpc(length_px, redshift, wcs_info, h0=h0, om0=om0)
    
    output_name = filename.replace(".fits", "_analysis.png")
    generate_plot(img, isolist, radius_px, size_kpc, output_name)
    
    return size_kpc

def save_debug_image(image_data: np.ndarray, filename_original: str) -> None:
    """
    Saves a debug image of the input data.

    Args:
        image_data: 2D galaxy image.
        filename_original: Original filename for naming the debug plot.
    """
    fig, ax = plt.subplots(figsize=(8, 8))
    
    good_mask = ~np.logical_or(np.isnan(image_data), image_data <= 0)
    norm = simple_norm(image_data[good_mask], stretch='log', percent=99.0) if np.any(good_mask) else None
    
    im = ax.imshow(image_data, cmap='magma', origin='lower', norm=norm)
    
    cbar = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.04)
    cbar.set_label('Flux (Log Scale)', fontsize=8)
    
    ax.set_title(f"INPUT DATA CHECK:\n{os.path.basename(filename_original)}", fontsize=14, fontweight='bold')
    ax.set_xlabel("Pixels", fontsize=12)
    ax.set_ylabel("Pixels", fontsize=12)

    debug_name = filename_original.replace(".fits", "_DEBUG_INPUT.png")
    plt.savefig(debug_name, dpi=150, bbox_inches='tight')
    plt.close()
    logger.info("Debug image saved to: %s", debug_name)
