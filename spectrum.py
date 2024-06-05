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

    def plot(self):
        plt.figure(figsize=(10, 6))
        plt.plot(self.mz, self.intensities, label=self.substance_name)
        plt.xlabel('m/z')
        plt.ylabel('Intensity')
        plt.title('Mass Spectrum')
        plt.legend()
        plt.show()

    def barplot(self):
        plt.figure(figsize=(10, 6))
        plt.bar(self.mz, self.intensities, width=0.01, label=self.substance_name)
        plt.xlabel('m/z')
        plt.ylabel('Intensity')
        plt.title('Mass Spectrum (Bar Plot)')
        plt.legend()
        plt.show()