import matplotlib.pyplot as plt
import numpy as np


x_vals = []
# gauss_vals = []
# bias_vals = []
mu_vals = []
total_data = np.loadtxt("pot_diff.csv", delimiter=",")




for row in total_data:
    mu_vals.append(row[0])
    x_vals.append(row[1])
    # gauss_vals.append(row[1])
    # bias_vals.append(row[2])
    # poten_vals.append(row[1])

# plt.plot(x_vals, gauss_vals)
# plt.plot(x_vals, bias_vals)
# plt.plot(x_vals)


for val in mu_vals:
    if val > 0:
        print val
plt.plot(mu_vals)

plt.show()
# for i in range(100):
#     gauss_data = np.loadtxt("gauss_"+str(i)+".csv", delimiter=",")
#     x_vals = []
#     gauss_vals = []
#     for row in gauss_data:
#         x_vals.append(row[0])
#         gauss_vals.append(row[1])
#
#     plt.plot(x_vals, gauss_vals)


# bin_locs = []
# bin_vals = []
# bin_data = np.loadtxt("bins.csv", delimiter=",")
# for row in bin_data:
#     bin_locs.append(row[0])
#     bin_vals.append(row[1])
#
# plt.hist(bin_vals, bins=bin_locs)
#
# plt.show()