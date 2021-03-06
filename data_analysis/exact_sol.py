import scipy.integrate as integrate
import numpy as np
import matplotlib.pyplot as plt

KT_SHIFT = 4.4
QUARTIC_SHIFT = 2.8256360858458973
DIFF_WIDTH_SHIFT = 2


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
    return (x + x_0) ** 4 - 4 * (x + x_0) * (x + x_0) - (x + x_0) + 2.618555980765;


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
    normalisation_const = integrate.quad(lambda x: np.exp(-poten_fun(x) / kT), -np.inf, np.inf, epsabs=0.0000001,
                                         epsrel=0.0000001)

    def eq_func(x):
        return np.exp(-poten_fun(x) / kT) / normalisation_const[0]

    return eq_func


def exact_energy_diff(poten_name, kT):
    'Calculates the exact energy difference'

    poten_fun = potential_selector(poten_name, kT)[0]

    def integrand_func(x):
        return np.exp(-poten_fun(x) / kT)

    P_left = integrate.quad(integrand_func, -np.inf, 0, epsabs=0.0000001, epsrel=0.0000001)
    P_right = integrate.quad(integrand_func, 0, np.inf, epsabs=0.0000001, epsrel=0.0000001)
    return -kT * np.log(P_left[0] / P_right[0])


def potential_selector(poten_name, kT):
    if poten_name == "KT":
        return [KT_Poten, KT_SHIFT]
    elif poten_name == "QUARTIC":
        return [QUARTIC_Poten, QUARTIC_SHIFT]
    elif poten_name == "DIFF_WIDTH":
        return [DIFF_WIDTH_Poten, DIFF_WIDTH_SHIFT]


def output_range(poten_name, poten_fun, kT_range):
    file = open(poten_name + "_exact.csv", "w")
    file.write("Temperature, Energy Diff\n")
    for kT in kT_range:
        file.write(str(kT) + ", " + str(exact_energy_diff(poten_fun, kT)) + "\n")
    file.close()


if __name__ == "__main__":
    print exact_energy_diff(QUARTIC_Poten, 0.5)
