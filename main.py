import os

import numpy as np
from matplotlib import pyplot as plt
from Analysis_Module import *
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
        velocity = (posns[i + 1] - posns[i]) / (
                1000 * (np.abs(times[i + 1] - times[i])))
        velocities.append(velocity)

    velocities.insert(0, 0)
    return velocities


def quick_plot(x: list, y: list, title: str):
    plt.figure(title)
    plt.plot(x, y, marker="o", ls='', label=title)
    plt.title(title)
    plt.legend()


def get_index_geq(lst, value):
    for i in range(len(lst)):
        if lst[i] >= value:
            return i


def process_velocity(file_path, t_range):
    times, positions = time_positions_extractor(file_path)
    i_start = get_index_geq(times, t_range[0])
    i_end = get_index_geq(times, t_range[1])
    velocities = velocities_generator(times,
                                      positions)
    velocities = velocities[i_start: i_end]
    return np.mean(velocities), np.std(velocities) / np.sqrt(len(velocities))


def velocities_from_data_sets(data_folder_path, t_range_dict):
    data_dict = {}
    for path in os.listdir(data_folder_path):
        key = path[-6:-4]
        print(key)
        data_dict[key] = (process_velocity(
            os.path.join(data_folder_path, path), t_range_dict[key]))
    return data_dict


def get_ball_width(file_path):
    width_dict = {}
    ball_num, m1, m2, m3 = np.loadtxt(file_path, delimiter=',', unpack=True,
                                      skiprows=2)
    mean_width = np.mean((m1, m2, m3), axis=0) / 1000
    mean_width_err = np.std((m1, m2, m3), axis=0) / 1000 * np.sqrt(3)
    for i in range(len(ball_num)):
        width_dict["B{}".format(int(ball_num[i]))] = (
            mean_width[i], mean_width_err[i])
    return width_dict


def v_correct(v_dict, w_dict, dim, dim_err):
    v_corr_dict = {}
    for ii in v_dict:
        lst = []
        dd = w_dict[ii][0] / dim
        dd_err = error_prop_multiplication(dd, [w_dict[ii], [dim, dim_err]])
        divisor = (1 - 2.104 * dd + 2.089 * dd ** 2)
        divisor_err = error_prop_addition([2.104 * dd_err, error_prop_exponent(dd, dd_err, 2)])
        v_corr = v_dict[ii][0] / (1 - 2.104 * dd + 2.089 * dd ** 2)
        lst.append(v_corr)
        lst.append(
            error_prop_multiplication(v_corr, [w_dict[ii], [divisor, divisor_err]])
        )
        v_corr_dict[ii] = lst
    return v_corr_dict


def inverse_function(x, a):
    return a * x ** 2


def inverse_sqrt_function(x, a):
    return a * np.sqrt(x)


if __name__ == "__main__":
    dimension_of_container = 9.5 / 100
    dimension_of_container_err = 0.1 / 100

    water_time_dict = {'B1': (1, 3), 'B2': (1, 2), 'B3': (1, 2),
                       'B4': (0.75, 1.75), 'B5': (1.2, 1.8)}
    glycerin_time_dict = {'B1': (20, 80), 'B2': (10, 30), 'B3': (5, 15),
                          'B4': (2, 6), 'B5': (1, 4)}
    water_width_dict = get_ball_width(
        "Ball Widths.xlsx - Ball Widths Glycerin.csv")
    print(water_width_dict)
    v_dict = v_correct(velocities_from_data_sets("Water_data", water_time_dict), water_width_dict,
                       dimension_of_container, dimension_of_container_err)
    print(v_dict)
    v_lst = []
    v_err_lst = []
    w_lst = []
    w_err_lst = []
    for i in water_width_dict:
        v_lst.append(v_dict[i][0])
        v_err_lst.append(v_dict[i][1])
        w_lst.append(water_width_dict[i][0])
        w_err_lst.append(water_width_dict[i][1])
    print(v_lst)
    plot_x_vs_y(w_lst, w_err_lst, v_lst, v_err_lst, "water", inverse_function)
    plot_x_vs_y(w_lst, w_err_lst, v_lst, v_err_lst, "water", inverse_sqrt_function)
    plt.show()
    # # glycerin ball 1
    # g_b1_times, g_b1_positions = time_positions_extractor(
    #     "Glycerin_data/GabrielleJunyuGlycerinNov142023_B1.txt")
    # g_b1_velocities = velocities_generator(g_b1_times, g_b1_positions)
    #
    # # print(g_b1_velocities)
    # quick_plot(g_b1_times[20:], g_b1_velocities[20:], "g ball 1")
    #
    # # glycerin ball 2
    # g_b2_times, g_b2_positions = time_positions_extractor(
    #     "Glycerin_data/GabrielleJunyuGlycerinNov142023_B2.txt")
    # g_b2_velocities = velocities_generator(g_b2_times, g_b2_positions)
    #
    # quick_plot(g_b2_times, g_b2_velocities, "g ball 2")
    # # print(g_b2_velocities)
    #
    # # glycerin ball 3
    # g_b3_times, g_b3_positions = time_positions_extractor(
    #     "Glycerin_data/GabrielleJunyuGlycerinNov142023_B3.txt")
    # g_b3_velocities = velocities_generator(g_b3_times, g_b3_positions)
    #
    # quick_plot(g_b3_times, g_b3_velocities, "g ball 3")
    # # print(g_b3_velocities)
    #
    # # glycerin ball 4
    # g_b4_times, g_b4_positions = time_positions_extractor(
    #     "Glycerin_data/GabrielleJunyuGlycerinNov142023_B4.txt")
    # g_b4_velocities = velocities_generator(g_b4_times, g_b4_positions)
    #
    # quick_plot(g_b4_times, g_b4_velocities, "g ball 4")
    # # print(g_b2_velocities)
    #
    # # glycerin ball 5
    # g_b5_times, g_b5_positions = time_positions_extractor(
    #     "Glycerin_data/GabrielleJunyuGlycerinNov142023_B5.txt")
    # g_b5_velocities = velocities_generator(g_b5_times, g_b5_positions)
    #
    # quick_plot(g_b5_times, g_b5_velocities, "g ball 5")
    # # print(g_b5_velocities)
    #
    # # water ball 1
    # w_b1_times, w_b1_positions = time_positions_extractor(
    #     "Water_data/GabrielleJunyuWaterNov22023_B1.txt")
    # w_b1_velocities = velocities_generator(w_b1_times, w_b1_positions)
    #
    # quick_plot(w_b1_times[1:], w_b1_velocities[1:], "w ball 1")
    # # print(w_b1_velocities)
    #
    # # water ball 2
    # w_b2_times, w_b2_positions = time_positions_extractor(
    #     "Water_data/GabrielleJunyuWaterNov22023_B2.txt")
    #
    # w_b2_velocities = velocities_generator(w_b2_times, w_b2_positions)
    #
    # quick_plot(w_b2_times[1:], w_b2_velocities[1:], "w ball 2")
    # # print(w_b2_velocities)
    #
    # # water ball 3
    # w_b3_times, w_b3_positions = time_positions_extractor(
    #     "Water_data/GabrielleJunyuWaterNov22023_B3.txt")
    # w_b3_velocities = velocities_generator(w_b3_times, w_b3_positions)
    #
    # quick_plot(w_b3_times[1:], w_b3_velocities[1:], "w ball 3")
    #
    # # water ball 4
    # w_b4_times, w_b4_positions = time_positions_extractor(
    #     "Water_data/GabrielleJunyuWaterNov22023_B4.txt")
    # w_b4_velocities = velocities_generator(w_b4_times, w_b4_positions)
    #
    # quick_plot(w_b4_times[1:], w_b4_velocities[1:], "w ball 4")
    #
    # # water ball 5
    # w_b5_times, w_b5_positions = time_positions_extractor(
    #     "Water_data/GabrielleJunyuWaterNov22023_B5.txt")
    # w_b5_velocities = velocities_generator(w_b5_times, w_b5_positions)
    #
    # quick_plot(w_b5_times[1:], w_b5_velocities[1:], "w ball 5")
    # # print(w_b5_velocities)
    # plt.show()
