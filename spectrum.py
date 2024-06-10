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
                label=f'{other_spectrum.substance_name} (Negative)'
               )
        plt.xlabel('m/z')
        plt.ylabel('Intensity')
        plt.title(f'Comparison of {self.substance_name} and {other_spectrum.substance_name}')
        plt.legend()
        plt.show()

    @staticmethod
    def cosine_measure(spectrum1, spectrum2):
        common_mz = np.unique(np.concatenate((spectrum1.mz, spectrum2.mz)))
        common_intensities1 = np.zeros(len(common_mz))
        common_intensities2 = np.zeros(len(common_mz))
        mask1 = np.isin(common_mz, spectrum1.mz)
        mask2 = np.isin(common_mz, spectrum2.mz)
        common_intensities1[mask1] = spectrum1.intensities[np.where(np.isin(spectrum1.mz, common_mz[mask1]))[0]]
        common_intensities2[mask2] = spectrum2.intensities[np.where(np.isin(spectrum2.mz, common_mz[mask2]))[0]]
        cos_measure = np.dot(common_intensities1, common_intensities2) / (
            np.linalg.norm(common_intensities1) * np.linalg.norm(common_intensities2)
        )
        return cos_measure

    def cosine_measures_list(self, spectra):
        cos_measures = []
        for spectrum in spectra:
            cos_measure = self.cosine_measure(self, spectrum)
            cos_measures.append({
                'cm': cos_measure,
                'id': spectrum.metadata.get('id'),
                'substance_name': spectrum.substance_name
            })

        return cos_measures

    @staticmethod
    def cosine_measures_matrix_plot(spectra):
        n = len(spectra)
        cosine_measures_matrix = np.zeros((n, n))
        names = [spectrum.substance_name for spectrum in spectra]  # Vertical and horizontal captions
        for i in range(n):
            for k in range(n):
                cos_measure = Spectrum.cosine_measure(spectra[i], spectra[k])
                cosine_measures_matrix[i, k] = cos_measure
        fig, ax = plt.subplots(figsize=(6.4, 4.8))
        ax.grid(False)
        vmin=0
        vmax=1
        plt.imshow(cosine_measures_matrix, cmap='gray', aspect='auto', vmin=vmin, vmax=vmax)
        for i in range(n):
            for k in range(n):
                cell_value = round(cosine_measures_matrix[i, k], 2)
                cell_color = 'white' if cell_value < (vmax - vmin)/2 else 'black'
                text = ax.text(i, k, cell_value, ha="center", va="center", color=cell_color, fontsize=8)  # Подписи ячеек
        plt.colorbar(label='Cosine measure')
        plt.title('Cosine measures matrix of spectra')
        plt.xticks(range(n), names, rotation=45, ha='right')
        plt.yticks(range(n), names)
        plt.tight_layout()
        plt.show()