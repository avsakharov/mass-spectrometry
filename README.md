# Mass Spectrometry Data Analysis

This repository contains code and scripts for processing and analyzing mass spectrometry data. It includes classes and methods for handling spectrum data, normalizing intensities, plotting spectra, and more.

## Table of Contents
+ [Overview](##Overview)
+ [Installation](##Installation)
+ [Usage](##Usage)
+ [Class Documentation](##Class-Documentation)
+ [Contributing](##Contributing)
+ [License](##License)

## Overview
The main goal of this project is to provide tools for efficient handling and analysis of mass spectrometry data. The repository includes a Spectrum class for representing individual spectra, methods for normalization and thresholding, and functions for plotting spectra.

## Installation
To get started, clone the repository and install the necessary dependencies. You can use `pip` to install the required packages.
```
git clone https://github.com/avsakharov/mass-spectrometry.git
cd mass-spectrometry
```

## Usage
Here is an example of how to use the Spectrum class:
```
from spectrum import Spectrum
```

### Example data
```
mz = [3.0, 3.5, 4.0, 4.5, 5.0]
intensities = [10, 20, 15, 5, 10]
substance_name = "Sample Substance"
metadata = {
    "id": "1",
    "molecular_formula": "H2O",
    "molecular_weight": 18.01528
}
```

### Create a Spectrum instance
```
spectrum = Spectrum(mz, intensities, substance_name, **metadata)
```

### Normalize the spectrum
```
spectrum.normalize()
```

### Plot the spectrum
```
spectrum.plot()
```

### Plot the spectrum as a bar chart
```
spectrum.barplot()
```

### Apply a threshold to the intensities
```
threshold = 0.05
spectrum.apply_threshold(threshold)
```

### Change the sampling step
```
step = 1
spectrum.resample(step)
```

## Class Documentation

### Spectrum
The Spectrum class represents a mass spectrum and includes various methods for processing and analyzing the data.

### Initialization
```
Spectrum(mz, intensities, substance_name, **metadata)
```
+ `mz`: Array of m/z values
+ `intensities`: Array of intensity values
+ `substance_name`: Name of the substance
+ `metadata`: Optional metadata such as id, molecular_formula, molecular_weight, etc.

### Methods
+ `normalize()`: Normalize the intensities to a range between 0 and 1.
+ `plot()`: Plot the spectrum.
+ `barplot()`: Plot the spectrum as a bar chart.
+ `apply_threshold(threshold)`: Set intensities below the threshold to zero.
+ `resample(step)`: Resample the spectrum with a new m/z step.

## Contributing
Contributions are welcome! Please feel free to submit a Pull Request or open an Issue to discuss any changes or improvements.

## License
This project is licensed under the MIT License. See the LICENSE file for details.
