"""ACI-RCDF 1 cubic metre of concrete mix without additives

missing!
! maximum values of coarse aggregate of any size
! granulometry limits for 3/4" coarse agg
"""
import numpy as np
from numpy import genfromtxt


def ordinary_portland(fc, granulometry_array, slump=10, std=50, cement_density=3.15):
    """computes quantity of components for a cubic metre of concrete mix without additives
            according to RCDF class 1 concrete via ACI 318 absolute volumes procedure
            NO AIR INCLUDED
            returns array (weights, volumes, volumes for a 50kg sack of cement.)
    """
    revenimiento = genfromtxt('revenimiento.csv', delimiter=',')
    max_size_agg = granulometry_array[1,1]

    t0 = 3.35
    fcr = max(fc - 35 + t0 * std, fc + t0 * std / np.sqrt(3))
    agua_renglon = np.abs(array-value)).argmin() 
    agua_columna = (revenimiento[2:5,0:1]-slump).argmin()
    contenido_agua = 