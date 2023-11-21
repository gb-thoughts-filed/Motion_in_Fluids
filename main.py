import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit


def time_positions_extractor(file: str):
    times = list(np.loadtxt(file, unpack=True, usecols=0, skiprows=2))
    positions = list(np.loadtxt(file, unpack=True, usecols=1, skiprows=2))
    return times, positions


def posn_extractor(file: str):
    positions = list(np.loadtxt(file, unpack=True, usecols=1, skiprows=2))
    return positions


def velocities_generator(times: list, posns: list):

    velocities = []

    for i in np.arange(len(times) - 1):
        velocity = (posns[i+1] - posns[i])/np.abs(times[i+1] - times[i])
        velocities.append(velocity)

    velocities.insert(0, 0)
    return velocities


def quick_plot(x: list, y: list, title: str):
    plt.plot(x, y, marker="o", ls='')
    plt.title(title)
    plt.show()


if __name__ == "__main__":

    # glycerin ball 1
    g_b1_times, g_b1_positions = time_positions_extractor(
        "GabrielleJunyuGlycerinNov142023_B1.txt")
    g_b1_velocities = velocities_generator(g_b1_times, g_b1_positions)

    # print(g_b1_velocities)
    quick_plot(g_b1_times[20:], g_b1_velocities[20:], "g ball 1")

    # glycerin ball 2
    g_b2_times, g_b2_positions = time_positions_extractor("GabrielleJunyuGlycerinNov142023_B2.txt")
    g_b2_velocities = velocities_generator(g_b2_times, g_b2_positions)

    quick_plot(g_b2_times, g_b2_velocities, "g ball 2")
    # print(g_b2_velocities)

    # glycerin ball 3
    g_b3_times, g_b3_positions = time_positions_extractor("GabrielleJunyuGlycerinNov142023_B3.txt")
    g_b3_velocities = velocities_generator(g_b3_times, g_b3_positions)

    quick_plot(g_b3_times, g_b3_velocities, "g ball 3")
    # print(g_b3_velocities)

    # glycerin ball 4
    g_b4_times, g_b4_positions = time_positions_extractor("GabrielleJunyuGlycerinNov142023_B4.txt")
    g_b4_velocities = velocities_generator(g_b4_times, g_b4_positions)

    quick_plot(g_b4_times, g_b4_velocities, "g ball 4")
    # print(g_b2_velocities)

    # glycerin ball 5
    g_b5_times, g_b5_positions = time_positions_extractor("GabrielleJunyuGlycerinNov142023_B5.txt")
    g_b5_velocities = velocities_generator(g_b5_times, g_b5_positions)

    quick_plot(g_b5_times, g_b5_velocities, "g ball 5")
    # print(g_b5_velocities)

    # water ball 1
    w_b1_times, w_b1_positions = time_positions_extractor("GabrielleJunyuWaterNov22023_B1.txt")
    w_b1_velocities = velocities_generator(w_b1_times, w_b1_positions)

    quick_plot(w_b1_times[1:], w_b1_velocities[1:], "w ball 1")
    # print(w_b1_velocities)

    # water ball 2
    w_b2_times, w_b2_positions = time_positions_extractor("GabrielleJunyuWaterNov22023_B2.txt")

    w_b2_velocities = velocities_generator(w_b2_times, w_b2_positions)

    quick_plot(w_b2_times[1:], w_b2_velocities[1:], "w ball 2")
    # print(w_b2_velocities)

    # water ball 3
    w_b3_times, w_b3_positions = time_positions_extractor("GabrielleJunyuWaterNov22023_B3.txt")
    w_b3_velocities = velocities_generator(w_b3_times, w_b3_positions)

    quick_plot(w_b3_times[1:], w_b3_velocities[1:], "w ball 3")

    # water ball 4
    w_b4_times, w_b4_positions = time_positions_extractor("GabrielleJunyuWaterNov22023_B4.txt")
    w_b4_velocities = velocities_generator(w_b4_times, w_b4_positions)

    quick_plot(w_b4_times[1:], w_b4_velocities[1:], "w ball 4")

    # water ball 5
    w_b5_times, w_b5_positions = time_positions_extractor("GabrielleJunyuWaterNov22023_B5.txt")
    w_b5_velocities = velocities_generator(w_b5_times, w_b5_positions)

    quick_plot(w_b5_times[1:], w_b5_velocities[1:], "w ball 5")
    # print(w_b5_velocities)
