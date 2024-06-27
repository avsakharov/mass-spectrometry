# Mass Spectrometry Data Processing

This repository contains tools and scripts for processing mass spectrometry data. The main features include loading and preprocessing spectra, background subtraction, resampling, and identifying substances based on spectral data.

## Table of Contents
+ [Overview](##Overview)
+ [Installation](##Installation)
+ [Usage](##Usage)
+ [Example](##Example)
+ [Contributing](##Contributing)
+ [License](##License)

## Overview
The main goal of this project is to provide tools for efficient handling and analysis of mass spectrometry data. The repository includes Spectrum and SpectrumLoader classes for download and representing individual spectra, methods for normalization and thresholding, and functions for plotting spectra.

## Installation

To get started with the project, clone the repository and install the required dependencies:
```bash
git clone https://github.com/avsakharov/mass-spectrometry.git
cd mass-spectrometry
```

## Usage

```
# Loading spectra. You can load spectra from `.mzml` files or `.txt` files using the `SpectrumLoader` class.
import pandas as pd
from spectrum import SpectrumLoader

# Load metadata
metadata_df = pd.read_csv('path_to_metadata.csv', index_col=0)

# Initialize the loader
loader = SpectrumLoader(step=1)

# Load spectra from the database
spectra_db = loader.load_spectra_from_mzml('data/database_trap/database', metadata_df)

# Load background spectrum
spectrum_bg = loader.load_background_spectrum('data/database_trap/background/Background.mzml')

# Load test spectra
spectra_test = loader.load_spectra_from_mzml('data/database_trap/test')

# Load spectra from text files
spectra_txt_db = loader.load_spectra_from_txt('results/database_trap')

# Background Subtraction.
# You can subtract background spectra from your database spectra using the `background_subtraction` method in the `Spectrum` class.
spectrum = spectra_db[0]
spectrum.background_subtraction(spectrum_bg)

# Spectrum Matching. To compare test spectra with the database, use the `calculate_metrics` method in teh `Spectrum` clss.
cos_measure, pearson_coefficient, euclid_distance, manhattan_distance = spectrum.calculate_metrics(spectrum_test[0])
```

## Example
```
# Create a Spectrum instance
mz = [3.0, 3.5, 4.0, 4.5, 5.0]
intensities = [10, 20, 15, 5, 10]
substance_name = "Sample Substance"
metadata = {
    "id": "1",
    "molecular_formula": "H2O",
    "molecular_weight": 18.01528
}
spectrum = Spectrum(mz, intensities, substance_name, **metadata)

# Normalize the spectrum
spectrum.normalize()

# Plot the spectrum
spectrum.plot()

# Plot the spectrum as a bar chart
spectrum.barplot()

# Comparative plot of two spectra
# To create a comparative plot of two spectra, where one spectrum is displayed with positive intensities and
# the other with negative, use the compare_plot method.
mz2 = [3.0, 3.5, 4.0, 4.5, 5.0]
intensities2 = [8, 18, 12, 4, 9]
substance_name2 = "Another Substance"
metadata2 = {"id": "2", "molecular_formula": "H2O2", "molecular_weight": 34.0147}

spectrum2 = Spectrum(mz2, intensities2, substance_name2, **metadata2)
spectrum2.normalize()
spectrum.compare_plot(spectrum2)

# Comparative bar plot of two spectra. To create a comparative bar plot of two spectra, use the compare_barplot method.
spectrum.compare_barplot(spectrum2)

# Apply a threshold to the intensities
threshold = 0.05
spectrum.apply_threshold(threshold)

# Change the sampling step
step = 1
spectrum.resample(step)

# Plot cosine measures matrix
spectra = [spectrum, spectrum2]
Spectrum.cosine_measures_matrix_plot(spectra)

# Cosine measures list. Get cosine measures between test_spectrum and spectrum in spectra.
test_spectrum.calculate_cosine_measures(spectra)
```

## Contributing
Contributions are welcome! Please feel free to submit a Pull Request or open an Issue to discuss any changes or improvements.

## License
This project is licensed under the MIT License. See the LICENSE file for details.
