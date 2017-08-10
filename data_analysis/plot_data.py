import matplotlib.pyplot as plt
import numpy as np
import exact_sol
import sys
import csv


def plot_solutions(poten_name):
    'Plots the solutions of the simulation against their exact answers'
    raw_data_arr = []
    data_reader = csv.reader(open(poten_name + "_combined.csv"), delimiter=",")
    next(data_reader, None)
    for row in data_reader:
        raw_data_arr.append(row)

    dynamics_type = raw_data_arr[0][1]

    useful_data_arr = []
    temps = set()
    for row in raw_data_arr:
        timestep = float(row[3])
        kT = float(row[4])
        temps.add(kT)
        energy_diff = float(row[5])
        std_error = float(row[6])
        useful_data_arr.append([kT, timestep, energy_diff, std_error])

    for kT in temps:
        timesteps = []
        energy_diffs = []
        std_error = []
        for row in useful_data_arr:
            if row[0] == kT:
                timesteps.append(row[1])
                energy_diffs.append(row[2])
                std_error.append(row[3])

        exact_answer = exact_sol.exact_energy_diff(poten_name, kT)

        fig, ax = plt.subplots()
        exact_line = ax.plot(timesteps, [exact_answer for p in timesteps], linestyle="dashed", color="green", label="Exact")
        calc_line = ax.plot(timesteps, energy_diffs, color="blue")
        calc_scatter = ax.scatter(timesteps, energy_diffs, color="blue", label="Lattice switch")
        calc_errs = ax.errorbar(timesteps, energy_diffs, yerr=std_error, color="blue", fmt="none")

        plt.xlabel("Timestep size")
        plt.ylabel("Free energy difference")
        plt.title(poten_name + ", " + dynamics_type + ". kT = " + str(kT))
        ax.legend(loc="best")
        plt.savefig(poten_name + "_" + str(kT) +".png")



if __name__ == "__main__":
    plot_solutions("DIFF_WIDTH")
    plot_solutions("QUARTIC")
    plot_solutions("KT")
