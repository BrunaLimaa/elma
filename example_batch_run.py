import os
import logging
from elma import run_pipeline, get_fits_filters

# Configure logging to show output in the console
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

def main():
    # Directory containing the FITS files
    data_dir = r"/Users/barred"

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
            filters = get_fits_filters(full_path)
            bar_size_kpc = run_pipeline(filename=full_path, redshift=z)

            results[filename] = {"bar_kpc": bar_size_kpc, "filters": filters}
            logger.info("Detected Bar Length for %s: %.2f kpc", filename, bar_size_kpc)

        except Exception as e:
            logger.error("Failed to process %s: %s", filename, e)

    logger.info("Batch processing complete.")
    print("\nSummary of Results:")
    print("-" * 60)
    for name, info in results.items():
        filter_str = ", ".join(f"{k}={v}" for k, v in info["filters"].items()) or "N/A"
        print(f"{name:15}: {info['bar_kpc']:8.2f} kpc  |  filters: {filter_str}")
    print("-" * 60)

if __name__ == "__main__":
    main()
