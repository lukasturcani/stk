#!/usr/bin/env python3
""" Module for plotting. Under construction. """

import json
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as stats


class Plot(object):
    def __init__(self, save_figures=True):
        a = 1

    @classmethod
    def load_data_dict(cls, data):
        obj = cls()
        obj.data = data
        return obj

    @classmethod
    def load_data_json(cls, filepath):
        obj = cls()
        with open(filepath, 'r') as file_src:
            obj.data = json.load(file_src)
        return obj

    def plot_PLE(self):
        self.all_windows = np.array([], dtype=float)
        for frame in self.data:
            for mol in self.data[frame].keys():
                self.all_windows = np.concatenate(
                    (
                        self.all_windows,
                        np.array(self.data[frame][mol]['windows'][
                            'windows_diameters'], dtype=float)
                    )
                )
        min_x = self.all_windows[np.argmin(self.all_windows)]
        max_x = self.all_windows[np.argmax(self.all_windows)]
        linsp = np.linspace(min_x, max_x, 1000)
        kernel = stats.gaussian_kde(self.all_windows)
        kde = kernel(linsp)

        fig, ax = plt.subplots(figsize=(7, 3.5))
        plt.plot(linsp, kde)
        plt.show()
