# -*- coding: utf-8 -*-
"""
Created on Mon Nov  5 17:41:55 2018

@author: henders
"""

import csv
# import matplotlib
import matplotlib.pyplot as plt
import matplotlib.lines as lines
import numpy as np


LASER = True
RADAR = True


def read_nis_data(file_name):
    nis_data = []
    with open('nis_laser.csv', newline='') as csvfile:
        nis_reader = csv.reader(csvfile)
        for row in nis_reader:
            for col in row:
                try:
                    nis_data.append(float(col))
                except ValueError:
                    continue
    return nis_data


def plot_nis_data(nis_data, five_pct_line, file_name):
    # The x axis will represent time.
    k = np.arange(0, len(nis_data))

    # Plot the NIS data.
    fig, ax = plt.subplots(figsize=(20, 10))
    ax.plot(k, nis_data)

    # Draw the 95% line, 95% of the data should be below this line
    l = lines.Line2D([0, len(nis_data)], [five_pct_line, five_pct_line])
    l.set_color('red')
    ax.add_line(l)

    # Set plot params.
    ax.set(xlabel='time (s)', ylabel='NIS',
           title='NIS')
    ax.grid()

    # Save and draw the plot.
    fig.savefig(file_name)
    plt.show()


if LASER:
    nis_laser = read_nis_data('nis_laser.csv')
    plot_nis_data(nis_laser, 5.991, 'nis_laser.png')

if RADAR:
    nis_radar = read_nis_data('nis_radar.csv')
    plot_nis_data(nis_radar, 7.815, 'nis_radar.png')
