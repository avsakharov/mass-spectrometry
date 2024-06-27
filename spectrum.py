# spectrum.py
import os
import glob
import numpy as np
import pandas as pd
from pyteomics import mzml  # Library used to open mzml format files (spectra format)

import matplotlib.pyplot as plt
import seaborn as sns
sns.set_theme()

class Spectrum:
    def __init__(self, mz, intensities, substance_name, **metadata):
        self.mz = np.array(mz)
        self.intensities = np.array(intensities)
        self.substance_name = substance_name
        self.metadata = metadata

    def normalize(self):
        # Min-max normalization
        min_intensity = np.min(self.intensities)
        max_intensity = np.max(self.intensities)
        self.intensities = (self.intensities - min_intensity) / (max_intensity - min_intensity)
    
    def resample(self, step=1):
        # Changing the sampling step
        mz_min = np.floor(self.mz.min())
        mz_max = np.ceil(self.mz.max())
        new_mz = np.arange(mz_min, mz_max + step, step)
        new_intensities = np.zeros(len(new_mz))
        for i in range(1, len(new_mz)-1):  # Group the intensities into a new array mz with step and sum up
            mask = (self.mz >= new_mz[i-1]) & (self.mz < new_mz[i])
            new_intensities[i] = self.intensities[mask].sum()
        self.mz = new_mz
        self.intensities = new_intensities
    
    def apply_threshold(self, threshold=0.001):  # Cutting off values below the threshold
        non_zero_mask = self.intensities >= threshold
        self.mz = self.mz[non_zero_mask]
        self.intensities = self.intensities[non_zero_mask]
    
    @staticmethod
    def merge_and_align_spectra(spectrum1, spectrum2):
        # Combine arrays spectrum1.mz and spectrum2.mz into array common_mz and get spectrum intensities
        common_mz = np.unique(np.concatenate((spectrum1.mz, spectrum2.mz)))
        intensities1 = np.zeros(len(common_mz))
        intensities2 = np.zeros(len(common_mz))
        mask1 = np.isin(common_mz, spectrum1.mz)
        mask2 = np.isin(common_mz, spectrum2.mz)
        intensities1[mask1] = spectrum1.intensities[np.where(spectrum1.mz == common_mz[mask1])[0]]
        intensities2[mask2] = spectrum2.intensities[np.where(spectrum2.mz == common_mz[mask2])[0]]
        spectrum1.mz = common_mz
        spectrum1.intensities = intensities1
        spectrum2.mz = common_mz
        spectrum2.intensities = intensities2
        return spectrum1, spectrum2

    def background_subtraction(self, spectrum_bg):  # Subtract the background spectrum from the database spectrum
        self, spectrum_bg = Spectrum.merge_and_align_spectra(self, spectrum_bg)
        self.intensities = np.maximum(self.intensities - spectrum_bg.intensities, 0)
        
    def to_dict(self):
        data = {
            'cas': self.metadata.get('cas', ''),
            'substance_name': self.substance_name,
            'other_names': self.metadata.get('other_names', ''),
            'molecular_formula': self.metadata.get('molecular_formula', ''),
            'monoisotopic_mass': self.metadata.get('monoisotopic_mass', ''),
            'array_length': len(self.mz),
            'mz': list(self.mz),
            'intensities': list(self.intensities)
        }
        return data
    
    def plot(self):
        plt.figure(figsize=(12, 4))
        plt.plot(self.mz, self.intensities, label=self.substance_name)
        plt.xlabel('m/z')
        plt.ylabel('Intensity')
        plt.title(f'Mass Spectrum {self.substance_name}')
        plt.legend()
        plt.show()

    def barplot(self):
        plt.figure(figsize=(12, 4))
        plt.bar(self.mz, self.intensities, width=0.5, edgecolor='none', label=self.substance_name)
        plt.xlabel('m/z')
        plt.ylabel('Intensity')
        plt.title(f'Mass Spectrum (Bar Plot) {self.substance_name}')
        plt.legend()
        plt.show()
    
    def compare_plot(self, other_spectrum):
        if not isinstance(other_spectrum, Spectrum):
            raise ValueError("The object to compare must be an instance of Spectrum class.")
        plt.figure(figsize=(12, 4))
        # Plotting the first spectrum with positive intensities
        plt.plot(self.mz, self.intensities, color='blue', label=f'{self.substance_name} (Positive)')
        # Plotting the second spectrum with negative intensities
        plt.plot(other_spectrum.mz, -other_spectrum.intensities, color='red', label=f'{other_spectrum.substance_name} (Negative)')
        plt.xlabel('m/z')
        plt.ylabel('Intensity')
        plt.title(f'Comparison of {self.substance_name} and {other_spectrum.substance_name}')
        plt.legend()
        plt.show()
    
    def compare_barplot(self, other_spectrum):
        if not isinstance(other_spectrum, Spectrum):
            raise ValueError("The object to compare must be an instance of Spectrum class.")
        plt.figure(figsize=(12, 4))
        # Plotting the first spectrum with positive intensities
        plt.bar(self.mz, self.intensities, width=0.5, color='blue', edgecolor='none', label=f'{self.substance_name} (Positive)')
        # Plotting the second spectrum with negative intensities
        plt.bar(other_spectrum.mz, -other_spectrum.intensities, width=0.5, color='red', edgecolor='none',
                label=f'{other_spectrum.substance_name} (Negative)'
               )
        plt.xlabel('m/z')
        plt.ylabel('Intensity')
        plt.title(f'Comparison of {self.substance_name} and {other_spectrum.substance_name}')
        plt.legend()
        plt.show()

    @staticmethod
    def cosine_measure(spectrum1, spectrum2):
        spectrum1, spectrum2 = Spectrum.merge_and_align_spectra(spectrum1, spectrum2)
        cos_measure = np.dot(spectrum1.intensities, spectrum2.intensities) / (
            np.linalg.norm(spectrum1.intensities) * np.linalg.norm(spectrum2.intensities)
        )
        return cos_measure

    def cosine_measures_list(self, spectra):
        cos_measures_list = []
        for spectrum in spectra:
            cos_measure = Spectrum.cosine_measure(self, spectrum)
            cos_measures_list.append({
                'cm': cos_measure,
                'id': spectrum.metadata.get('id'),
                'substance_name': spectrum.substance_name
            })

        return cos_measures_list

    @staticmethod
    def cosine_measures_matrix_plot(spectra):
        n = len(spectra)
        cosine_measures_matrix = np.zeros((n, n))
        names = [spectrum.substance_name for spectrum in spectra]  # Vertical and horizontal captions
        for i in range(n):
            for k in range(n):
                cos_measure = Spectrum.cosine_measure(spectra[i], spectra[k])
                cosine_measures_matrix[i, k] = cos_measure
        fig, ax = plt.subplots(figsize=(12, 9))
        ax.grid(False)
        vmin=0
        vmax=1
        plt.imshow(cosine_measures_matrix, cmap='copper', aspect='equal', vmin=vmin, vmax=vmax)
        for i in range(n):
            for k in range(n):
                cell_value = round(cosine_measures_matrix[i, k], 2)
                cell_color = 'white' if cell_value < (vmax - vmin)/2 else 'black'
                text = ax.text(i, k, cell_value, ha="center", va="center", color=cell_color, fontsize=8)  # Cell signatures
        plt.colorbar(label='Cosine measure')
        plt.title('Cosine measures matrix of spectra')
        plt.xticks(range(n), names, rotation=45, ha='right')
        plt.yticks(range(n), names)
        plt.tight_layout()
        plt.show()


class SpectrumLoader:
    def __init__(self, step=1):
        self.step = step

    def load_spectra_from_mzml(self, folder_path, metadata_df=None):
        mzml_files = glob.glob(os.path.join(folder_path, '*.mzml'))
        spectra = []

        for file_path in mzml_files:
            file_name = os.path.splitext(os.path.basename(file_path))[0]
            ms_data = mzml.read(file_path)
            spectrum_data = list(ms_data)[0]

            mz = spectrum_data.pop('m/z array')
            intensities = spectrum_data.pop('intensity array')

            if metadata_df is not None:
                metadata = metadata_df.loc[file_name].to_dict()
                substance_name = metadata.pop('substance_name')
            else:
                metadata = {'file_name': file_name}
                substance_name = 'Unknown substance'

            spectrum = Spectrum(mz, intensities, substance_name, **metadata)
            spectrum.resample(self.step)
            spectra.append(spectrum)

        return spectra

    def load_background_spectrum(self, file_path):
        ms_data = mzml.read(file_path)
        spectrum_data = list(ms_data)[0]

        mz = spectrum_data.pop('m/z array')
        intensities = spectrum_data.pop('intensity array')
        metadata = {}
        substance_name = 'Air'

        spectrum = Spectrum(mz, intensities, substance_name, **metadata)
        spectrum.resample(self.step)

        return spectrum

    def load_spectra_from_txt(self, folder_path):
        file_list = os.listdir(folder_path)
        spectra = []

        for file in file_list:
            if file.endswith('.txt'):
                file_path = os.path.join(folder_path, file)
                with open(file_path, 'r') as f:
                    lines = f.readlines()
                    file_data = {}
                    for line in lines:
                        key, value = line.strip().split(': ')
                        if key in ['mz', 'intensities']:
                            value = [float(v) for v in value[1:-1].split(', ')]
                        elif key == 'array_length':
                            value = int(value)
                        file_data[key] = value

                    mz = file_data.pop('mz')
                    intensities = file_data.pop('intensities')
                    substance_name = file_data.pop('substance_name')

                    spectrum = Spectrum(mz, intensities, substance_name, **file_data)
                    spectra.append(spectrum)

        return spectra
