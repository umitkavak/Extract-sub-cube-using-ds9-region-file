
# Sub-cube and Moment Calculation for NGC 7538

This repository contains a Python script to create a sub-cube from a specified region in a FITS file and calculate the moments (integrated intensity, velocity, and FWHM) for the NGC 7538 region. The script handles region extraction, sub-cube creation, and moment calculations.

## Requirements

- Python 3.x
- Astropy
- Numpy
- Spectral-cube
- Astroquery
- Matplotlib
- Regions

You can install the necessary Python packages using:

```sh
pip install astropy numpy spectral-cube astroquery matplotlib regions
```

## Script Overview

The script `subcube_moment_calculation.py` performs the following steps:

1. **Open the FITS File**: Reads the FITS file and extracts the data and header.
2. **Read Regions and Create Sub-cube**: Reads regions from a DS9 region file, converts them to pixel coordinates, and creates a sub-cube.
3. **Write the Sub-cube to a New FITS File**: Saves the created sub-cube to a new FITS file.
4. **Define Spectral Slab Range**: Specifies the range of the spectral slab.
5. **Calculate Moments**: Calculates the integrated intensity, velocity, and FWHM moments.
6. **Write Moments to FITS Files**: Saves the calculated moments to separate FITS files.

## Detailed Steps

### 1. Open the FITS File

```python
from astropy.io import fits
from spectral_cube import SpectralCube

# Open the FITS file for reading
hi_datafile = 'NGC7538_CII_merged_fullmap.fits'
try:
    hi_data = fits.open(hi_datafile)
    cube = SpectralCube.read(hi_data)
except FileNotFoundError:
    print(f"Error: The file {hi_datafile} was not found.")
    raise
except Exception as e:
    print(f"Error opening FITS file: {e}")
    raise
```

### 2. Read Regions and Create Sub-cube

```python
from regions import Regions

# Read regions and create sub-cube
region_file = 'region_ngc7538_2.reg'
try:
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
```

### 3. Write the Sub-cube to a New FITS File

```python
sub_cube_file = 'NGC7538_CII_subcube.fits'  # New sub-cube output file

# Write the sub-cube to a new FITS file
sub_cube.write(sub_cube_file, overwrite=True)
```

### 4. Define Spectral Slab Range

```python
import astropy.units as u

# Define spectral slab range
spectral_range_min = -130.0 * u.km / u.s
spectral_range_max = 25.0 * u.km / u.s
sub_cube_slab = sub_cube.spectral_slab(spectral_range_min, spectral_range_max)
```

### 5. Calculate Moments

```python
# Calculate moments
moment_0 = sub_cube_slab.with_spectral_unit(u.km / u.s).moment(order=0)  # Integrated intensity
moment_1 = sub_cube_slab.with_spectral_unit(u.km / u.s).moment(order=1)  # Velocity
moment_2 = sub_cube_slab.with_spectral_unit(u.km / u.s).moment(order=2)  # FWHM
```

### 6. Write Moments to FITS Files

```python
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
```

## Running the Script

Ensure that the FITS file `NGC7538_CII_merged_fullmap.fits` and the region file `region_ngc7538_2.reg` are in the same directory as the script. Run the script using:

```sh
python subcube_moment_calculation.py
```

## Results
sub-cube
<img width="921" alt="Screenshot 2024-06-12 at 22 55 00" src="https://github.com/umitkavak/Extract-sub-cube-using-ds9-region-file/assets/26542534/a49815a8-46be-4c7e-8373-c137f65c6878">

Moment maps
![Screenshot 2024-06-12 at 23 01 30](https://github.com/umitkavak/Extract-sub-cube-using-ds9-region-file/assets/26542534/f355695c-3a02-4c73-b35e-4534b6f85e47)

## References

- Astropy Documentation: [https://docs.astropy.org/]

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
