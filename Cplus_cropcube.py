import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.colors as colors
import matplotlib as mpl

import astropy.units as u
from astropy.io import fits
from spectral_cube import SpectralCube

from astroquery.esasky import ESASky
from astropy.wcs import WCS
from reproject import reproject_interp

from matplotlib import rc
from regions import Regions

# Font configuration for plots
rc('font', **{'family': 'serif', 'serif': ['Times']})
rc('text', usetex=True)
plt.rcParams.update({'font.size': 20})

# Input file names
hi_datafile = 'NGC7538_CII_merged_fullmap.fits'
region_file = 'region_ngc7538_2.reg'
sub_cube_file = 'NGC7538_CII_subcube.fits'  # New sub-cube output file

try:
    # Open the FITS file for reading
    hi_data = fits.open(hi_datafile)
    cube = SpectralCube.read(hi_data)
except FileNotFoundError:
    print(f"Error: The file {hi_datafile} was not found.")
    raise
except Exception as e:
    print(f"Error opening FITS file: {e}")
    raise

try:
    # Read regions and create sub-cube
    region_list = Regions.read(region_file, format='ds9')

    # Convert the regions to pixel coordinates if they are in fk5 coordinates
    wcs = cube.wcs.sub(['longitude', 'latitude'])
    pixel_regions = []
    for region in region_list:
        if hasattr(region, 'to_pixel'):
            pixel_region = region.to_pixel(wcs)
            pixel_regions.append(pixel_region)
        else:
            raise ValueError("Region could not be converted to pixel coordinates")

    # Define a mask from the pixel regions
    mask = np.zeros(cube.shape[1:], dtype=bool)
    for region in pixel_regions:
        mask |= region.to_mask().to_image(cube.shape[1:]) > 0

    # Apply the mask to the cube to create the sub-cube
    sub_cube = cube.with_mask(mask)
except FileNotFoundError:
    print(f"Error: The region file {region_file} was not found.")
    raise
except TypeError as e:
    print(f"Error reading region file: {e}")
    raise
except Exception as e:
    print(f"Unexpected error: {e}")
    raise

# Write the sub-cube to a new FITS file
sub_cube.write(sub_cube_file, overwrite=True)

# Define spectral slab range
spectral_range_min = -130.0 * u.km / u.s
spectral_range_max = 25.0 * u.km / u.s
sub_cube_slab = sub_cube.spectral_slab(spectral_range_min, spectral_range_max)

# Calculate moments
moment_0 = sub_cube_slab.with_spectral_unit(u.km / u.s).moment(order=0)  # Integrated intensity
moment_1 = sub_cube_slab.with_spectral_unit(u.km / u.s).moment(order=1)  # Velocity
moment_2 = sub_cube_slab.with_spectral_unit(u.km / u.s).moment(order=2)  # FWHM

# Output file names
moment_0_file = 'NGC7538_CII_moment0.fits'
moment_1_file = 'NGC7538_CII_moment1.fits'
moment_2_file = 'NGC7538_CII_moment2.fits'

# Write the moments as FITS images
moment_0.write(moment_0_file, overwrite=True)
moment_1.write(moment_1_file, overwrite=True)
moment_2.write(moment_2_file, overwrite=True)

print(f"Sub-cube written to: {sub_cube_file}")
print(f"Moment files written: {moment_0_file}, {moment_1_file}, {moment_2_file}")
