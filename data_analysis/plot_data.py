import matplotlib.pyplot as plt
import numpy as np
import exact_sol
import sys


def plot_solutions(poten_name):
    raw_data = np.loadtxt(poten_name + "_combined.csv", delimiter=",", skiprows=1)
    print raw_data

if __name__ == "__main__":
    plot_solutions(sys.argv[1])
