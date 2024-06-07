# spectrum.py
import numpy as np
import matplotlib.pyplot as plt

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
    
    def resample(self, step):
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
    
    def apply_threshold(self, threshold):  # Cutting off values below the threshold
        self.intensities[self.intensities < threshold] = 0

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
                label=f'{other_spectrum.substance_name} (Negative)')
        plt.xlabel('m/z')
        plt.ylabel('Intensity')
        plt.title(f'Comparison of {self.substance_name} and {other_spectrum.substance_name}')
        plt.legend()
        plt.show()