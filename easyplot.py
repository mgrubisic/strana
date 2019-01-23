import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd
import seaborn as sns
import os
sns.set()
plt.figure()

def easyplot(arraylike_x, arraylike_y, xlabel = 'X', ylabel = 'Y', title = 'my plot', save_as = False, scatter=False, linestyle = '-', color = 'black', grid = True, legend = None, xlim = None, ylim = None, xticks = None, yticks = None, scaled = False):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    if scaled == True:
        plt.axis('scaled')
    # ax.set_xticks(xticks)
    # ax.set_yticks(yticks)
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    # plt.grid(grid)
    ax.set_title(title)
    if scatter == True:
        ax.scatter(arraylike_x, arraylike_y, color = color, s = 20)
    else:
        ax.plot(arraylike_x, arraylike_y, color = color, ls = linestyle)
    if save_as is not False:
        fig.savefig(save_as + '.pdf')
    plt.clf()
    return


def fileplot(save_as, labels=['x','y'], merge=False, *filenames):
    """plots multiple files as xy,
    can merge 2 files for disp-force plot"""
    plt.xlabel(labels[0])
    plt.ylabel(labels[1])
    plt.savefig(save_as)

def sasid_plot(filename, save_as, skip_rows = 22):
    """ sasid spectra plot
    """
    df = pd.read_table(filename, sep="\t", skiprows = skip_rows, names = ["period", "Q=1"])
    df["Q=2"] = df["Q=1"]*0.5
    ax = df.plot(x="period", y=["Q=2"])
    ax.set_xlabel(r'$T_1$ s')
    ax.set_ylabel(r'$A/g$')
    # ax.legend_.remove()
    plt.savefig(save_as)

# def spec

sp = sasid_plot("spectra1.txt", "spectraq1.png")
