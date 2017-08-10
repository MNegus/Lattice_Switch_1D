import sys


def combine(poten_name):
    fout = open("combined.csv", "w")

    fout.write("Potential Name, Dynamics Type, No of steps, Timestep, kT, Free energy diff, Std error\n")

    for num in range(1, 166):
        f = open(poten_name + "_" + str(num) + "_data.csv")
        for line in f:
            fout.write(line)
        f.close()

    fout.close()


if __name__ == "__main__":
    combine(sys.argv[1])