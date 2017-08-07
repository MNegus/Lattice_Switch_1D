import scipy.integrate as integrate
import numpy as np
import matplotlib.pyplot as plt

def KT_Poten(x):
	'Potential function from kinetic theory notes'
	if (x <= -1):
		return 5 * (x + 2) * (x + 2);
	elif (x <= 1.2):
		return 10 - 5 * x * x;
	else:
		return 5 * (x - 2.4) * (x - 2.4) - 4.4;


def QUARTIC_Poten(x):
	'A quartic polynomial potential'
	x_0 = -0.126000192586256;
	return (x + x_0)**4 - 4 * (x + x_0) * (x + x_0) - (x + x_0) + 2.618555980765;


def DIFF_WIDTH_Poten(x):
	'Potential with differing well widths'
	if (x <= -1):
		return 5 * (x + 2) * (x + 2);
	elif (x <= 0):
		return 10 - 5 * x * x;
	elif (x <= 0.3461):
		return -50 * x * x + 10;
	else:
		return 50 * (x - np.sqrt(0.48)) * (x - np.sqrt(0.48)) - 2;


def equilibrium_sol(poten_fun, kT):
	'Returns the equilibrium solution function of a given potential and temperature'
	normalisation_const = integrate.quad(lambda x: np.exp(-poten_fun(x) / kT), -np.inf, np.inf, epsabs=0.0000001, epsrel=0.0000001)
	def eq_func(x):
		return np.exp(-poten_fun(x) / kT) / normalisation_const[0]
	return eq_func


def exact_energy_diff(poten_fun, kT):
	'Calculates the exact energy difference'
	def integrand_func(x):
		return np.exp(-poten_fun(x) / kT)

	P_left =  integrate.quad(integrand_func, -np.inf, 0, epsabs=0.0000001, epsrel=0.0000001)
	P_right = integrate.quad(integrand_func, 0, np.inf,  epsabs=0.0000001, epsrel=0.0000001)
	return -kT * np.log(P_left[0] / P_right[0])


def output_range(poten_name, poten_fun, kT_range):
	file = open(poten_name + "_exact.csv", "w")
	file.write("Temperature, Energy Diff\n")
	for kT in kT_range:
		file.write(str(kT) + ", " + str(exact_energy_diff(poten_fun, kT)) + "\n")
	file.close()

if __name__ == "__main__":
	eq_fun = equilibrium_sol(QUARTIC_Poten, 0.1)
	x_left = np.linspace(-2, 0, 1000)
	x_right = np.linspace(0, 2, 1000)
	y_left = eq_fun(x_left)
	y_right = eq_fun(x_right)

	plt.plot(x_left, y_left)
	plt.show()
	plt.plot(x_right, y_right)
	plt.show()
