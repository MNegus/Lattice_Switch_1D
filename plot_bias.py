import matplotlib.pyplot as plt
import numpy as np

raw_data = np.loadtxt("biased_data.csv", delimiter=",")

x_vals = []
y_vals = []
for row in raw_data:
    x_vals.append(row[0])
    y_vals.append(row[1])
    print row[1]
plt.plot(x_vals, y_vals)
plt.show()