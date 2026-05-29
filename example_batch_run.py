import os
import logging
from astropy.io import fits
from elma import run_pipeline

# Configure logging to show output in the console
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

def inspect_fits_filters(filepath):
    with fits.open(filepath) as hdul:
        print(f"\nFITS structure for: {os.path.basename(filepath)}")
        for i, hdu in enumerate(hdul):
            header = hdu.header
            print(f"  Extension {i} ({hdu.name}): shape={getattr(hdu.data, 'shape', None)}")
            for key in ("FILTER", "FILTER1", "FILTER2", "FILTER3", "FILTNAM1", "FILTNAM2", "FILTNAM3", "PHOTFILT", "BAND"):
                if key in header:
                    print(f"    {key} = {header[key]}")


def main():
    # Directory containing the FITS files
    data_dir = r"/Users/I774546/personal/barred"
    
    # Dictionary of files and their respective redshift (z) values
    fits_targets = {
        "barred1.fits": 0.42,
        "barred2.fits": 0.62,
        "barred3.fits": 0.25,
        "barred4.fits": 0.74,
        "barred5.fits": 0.67,
        "barred6.fits": 0.46
    }
    
    logger.info("Starting batch processing of %d files...", len(fits_targets))
    
    results = {}
    
    for filename, z in fits_targets.items():
        full_path = os.path.join(data_dir, filename)
        
        if not os.path.exists(full_path):
            logger.error("File not found: %s", full_path)
            continue
            
        try:
            logger.info("--- Processing %s (z=%.2f) ---", filename, z)
            inspect_fits_filters(full_path)
            # Run the pipeline
            bar_size_kpc = run_pipeline(filename=full_path, redshift=z)
            
            results[filename] = bar_size_kpc
            logger.info("Detected Bar Length for %s: %.2f kpc", filename, bar_size_kpc)
            
        except Exception as e:
            logger.error("Failed to process %s: %s", filename, e)

    logger.info("Batch processing complete.")
    print("\nSummary of Results:")
    print("-" * 40)
    for name, size in results.items():
        print(f"{name:15}: {size:8.2f} kpc")
    print("-" * 40)

if __name__ == "__main__":
    main()
