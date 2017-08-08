import matplotlib.pyplot as plt
import numpy as np
import exact_sol

def plot_bins(bins_filename, poten_name, kT):
    'Takes in an output file specifiying the bin locations and the number of particles in each'
    raw_bin_data = np.loadtxt(bins_filename, delimiter=",")
    poten_selection = exact_sol.potential_selector(poten_name, kT)
    poten_func = poten_selection[0]
    poten_shift = poten_selection[1]

    x_pos = []
    bar_heights = []
    for line in raw_bin_data:
        x_pos.append(line[0])
        bar_heights.append(line[1])

    nobins = len(x_pos)
    bar_width = x_pos[1] - x_pos[0]

    tot_area = 0

    for i in range(nobins):
        if x_pos[i] < 0:
            bar_heights[i] *= np.exp(-poten_shift / kT)

        tot_area += bar_width * bar_heights[i]

    for i in range(nobins):
        bar_heights[i] = bar_heights[i] / tot_area


    fig = plt.figure()

    ax1 = fig.add_subplot(111, label="1")
    ax2 = fig.add_subplot(111, label="2", frame_on=False)
    # ax3 = fig.add_subplot(111, label="3", frame_on=False)

    l1 = ax1.bar(x_pos, bar_heights, width=0.01, color="black", label="Bins")
    ax1.set_ylabel("Distribution of particles")

    potential_y = [poten_func(x) for x in x_pos]
    l2 = ax2.plot(x_pos, potential_y, label="Potential energy")
    ax2.yaxis.set_label_position('right')
    ax2.yaxis.tick_right()
    ax2.set_ylabel("Potential energy")

    equi_func = exact_sol.equilibrium_sol(poten_func, kT)
    equilibrium_y = [equi_func(x) for x in x_pos]
    l3 = ax1.plot(x_pos, equilibrium_y, color="red", label="Boltzmann distribution")
    ax1.legend(loc="upper left")
    ax2.legend()
    plt.show()

if __name__ == "__main__":
    plot_bins("my_output_bins.csv", "DIFF_WIDTH", 0.5)
