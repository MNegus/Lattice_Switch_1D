import matplotlib.pyplot as plt
import numpy as np

def plot_bins(bins_filename):
    'Takes in an output file specifiying the bin locations and the number of particles in each'
    raw_bin_data = np.loadtxt(bins_filename, delimiter=",")
    x_pos = []
    bar_heights = []
    for line in raw_bin_data:
        x_pos.append(line[0])
        bar_heights.append(line[1])

    plt.bar(x_pos, bar_heights, width=0.01)
    plt.show()

if __name__ == "__main__":
    plot_bins("my_output_bins.csv")
