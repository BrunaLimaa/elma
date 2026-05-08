import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.visualization import simple_norm
from astropy.wcs import WCS
from astropy import wcs
from astropy.cosmology import LambdaCDM
import astropy.units as u
from photutils.isophote import EllipseGeometry, Ellipse
from matplotlib.patches import Ellipse as MplEllipse


"""
Processes the FITS file. Extracts the image and the header info. 
Return 2D image data and the WCS. 
"""
def load_and_process(filename):
    
    with fits.open(filename) as hdul:
        
        header = hdul[0].header
        wcs_info = WCS(header)
        data = hdul[0].data
        
        if data.ndim == 3:
            print(f"[{filename}] Detected 3D cube. Collapsing to 2D")
            image_2d = np.sum(data, axis=0)
        else:
            print(f"[{filename}] Detected 2D Image.")
            image_2d = data
        
    return image_2d.astype(float), wcs_info



"""
Converts pixel length to physical distance in Kpc.
"""
def calculate_kpc(bar_length_pixels, redshift, wcs_info):
    
    cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
    
    scales = wcs.utils.proj_plane_pixel_scales(wcs_info)
    deg_per_pixel = scales[0]
    
    bar_size_deg = bar_length_pixels * deg_per_pixel
    bar_size_rad = np.deg2rad(bar_size_deg)
    
    distance_mpc = cosmo.angular_diameter_distance(redshift)
    distance_kpc = distance_mpc.to(u.kpc).value
    
    return bar_size_rad * distance_kpc



"""
Find bar logic. Returns everything needed for the plot.
"""
def find_bar_geometry(image_data):
    y_center, x_center = np.unravel_index(np.argmax(image_data), image_data.shape)
    
    g = EllipseGeometry(x0=x_center, y0=y_center, sma=5.0, eps=0.2, pa=0.0)
    ellipse = Ellipse(image_data, geometry=g)
    
    isolist = ellipse.fit_image(maxsma=450, fix_center=True, step=0.05)
    
    valid_fits = isolist.sma > 3
    
    if np.any(valid_fits):
        idx_peak = np.argmax(isolist.eps[valid_fits])
        bar_radius_pixels = isolist.sma[valid_fits][idx_peak]
    else:
        bar_radius_pixels = 0.0
        print("Aviso: Nenhuma elipse válida encontrada.")
        
    return isolist, bar_radius_pixels


"""
Creates the plot with the galaxy + red ellipse + size label.
Saves it to a file.
"""
def generate_plot(image_data, isolist, bar_radius_pixels, kpc_size, filename_out):
    good_mask = ~np.logical_or(np.isnan(image_data), image_data <= 0)
    norm = simple_norm(image_data[good_mask], stretch='log', percent=98.5)
    
    fig, ax = plt.subplots(figsize=(10, 10))
    im = ax.imshow(image_data, cmap='gray_r', origin='lower', norm=norm)
    
    for iso in isolist:
        if iso.sma < 3: continue 

        is_bar = np.isclose(iso.sma, bar_radius_pixels, atol=0.5)
        
        color = 'red' if is_bar else 'blue'
        alpha = 1.0 if is_bar else 0.3
        linewidth = 2.5 if is_bar else 1.0
        
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
    ax.text(0.02, 0.98, info_text, transform=ax.transAxes, color='red', 
            fontsize=12, verticalalignment='top', fontweight='bold',
            bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
    
    plt.title(f"Analysis: {filename_out}")
    
    plt.savefig(filename_out, dpi=150)
    plt.close() 
    print(f"Plot saved to: {filename_out}")


"""
Runs full pipeline.
"""
def run_pipeline(filename, redshift):
    
    img, wcs = load_and_process(filename)
    
    isolist, radius_px = find_bar_geometry(img)
    
    save_debug_image(img, filename)
    
    length_px = radius_px * 2
    size_kpc = calculate_kpc(length_px, redshift, wcs)
    
    output_name = filename.replace(".fits", "_analysis.png")
    generate_plot(img, isolist, radius_px, size_kpc, output_name)
    
    return size_kpc


"""
Save original image.
"""

def save_debug_image(image_data, filename_original):
    
    fig, ax = plt.subplots(figsize=(8, 8))
    
    good_mask = ~np.logical_or(np.isnan(image_data), image_data <= 0)
    norm = simple_norm(image_data[good_mask], stretch='log', percent=98.5) if np.any(good_mask) else None
    
    im = ax.imshow(image_data, cmap='gray_r', origin='lower', norm=norm)
    
    cbar = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.04)
    cbar.set_label('Flux (Log Scale)', fontsize=8)
    
    ax.set_title(f"INPUT DATA CHECK:\n{filename_original}", fontsize=14, fontweight='bold')
    ax.set_xlabel("Pixels", fontsize=12)
    ax.set_ylabel("Pixels", fontsize=12)

    debug_name = filename_original.replace(".fits", "_DEBUG_INPUT.png")
    plt.savefig(debug_name, dpi=150, bbox_inches='tight') # bbox_inches='tight' cuts extra whitespace
    plt.close()
    print(f"Debug image saved to: {debug_name}")