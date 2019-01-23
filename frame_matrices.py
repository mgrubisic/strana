"""
builds a set of regular frames from scratch and computes
their stiffness matrix

BUILDING = [height_n column_n trabe_n columna_n
			...			...			...
			....		....		....
			h2			c2 		t2 		c2
			h1 			c1 		t1       c1]

			cn,tn,hn, contain fundamental information about buildings

"""
import numpy as np

storeys = 2
bays = 1
storey_height = 3
bay_length = 5
height_ratio = 2


def Coordinates(storeys, bays, storey_height, bay_length, height_ratio):
    nodes = np.zeros((2, (storeys + 1)*(bays + 1)))
    x = 0
    y = 0
    counter = 0
    for i in range(bays + 1):
        for j in range(storeys + 1):
            nodes[:, counter] = [x, y]
            y += storey_height
            if j == 0:
                y = storey_height*height_ratio
            counter += 1
        y = 0
        x += bay_length
    print(nodes)
    return nodes

Coordinates(storeys, bays, storey_height, bay_length, height_ratio)
# guess element sections and properties and assign to nodes

# compute stiffness matrix

import numpy as np
from numpy.linalg import inv


def stiffness_matrix(angular_moment_beam, angular_moment_column, lateral_moment_column,
                     lateral_shear_column):
    matrix = np.array([[angular_moment_beam + angular_moment_column, angular_moment_beam / 2, lateral_moment_column],
                       [angular_moment_beam / 2, angular_moment_beam +
                           angular_moment_column, lateral_moment_column],
                       [lateral_moment_column, lateral_moment_column, 2 * lateral_shear_column]])
    return matrix


def stiffness_matrix2(angular_moment_beam, angular_moment_column, lateral_moment_column, lateral_shear_column):
    matrix = np.array([[angular_moment_beam + angular_moment_column, angular_moment_beam / 2, lateral_moment_column],
                       [angular_moment_beam / 2, angular_moment_beam, 0],
                       [lateral_moment_column, 0, 0]])
    return matrix


def angular_matrix(ang_stiffness):
    matrix = np.array([[ang_stiffness, ang_stiffness / 2],
                       [ang_stiffness / 2, ang_stiffness]])
    return matrix


def angular_matrix2(ang_stiffness):
    matrix = np.array([[0, 0 / 2],
                       [0 / 2, 0]])
    return matrix


def lateral_stiffness(elastic_modulus, length, inertia):
    ki = 12 * elastic_modulus * inertia / length**3
    kj = 6 * elastic_modulus * inertia / length**2
    return ki, kj


def angular_stiffness(elastic_modulus, length, inertia):
    kii = 4 * elastic_modulus * inertia / length
    return kii


def inertia_rectangle(base, height):
    inertia = base * height**3 / 12
    return inertia


def moment_AB(kij, displ_A, displ_B):  # kij is 2x2 array
    moments = np.dot(kij, np.array([displ_A, displ_B]))
    return moments

end_moment = 20 * 600**2 / 12
moments_state1 = np.array([end_moment, -end_moment])

E = 14000 * np.sqrt(250)
shear, m3 = lateral_stiffness(E, 300, inertia_rectangle(35, 35))
moment_beam = angular_stiffness(E, 600, inertia_rectangle(30, 60))
moment_column = angular_stiffness(E, 300, inertia_rectangle(35, 35))

Forces = np.array([-end_moment, end_moment, 10000])

displacements = np.dot(
    inv([stiffness_matrix(moment_beam, moment_column, m3, shear)]), Forces)

M_bar1 = moment_AB(angular_matrix(moment_column), displacements[
                   0, 2] / 300, displacements[0, 0] + displacements[0, 2] / 300)
M_bar2 = moment_AB(angular_matrix(moment_beam),
                   displacements[0, 0], displacements[0, 1])
M_bar3 = moment_AB(angular_matrix(moment_column), displacements[
                   0, 2] / 300, displacements[0, 1] + displacements[0, 2] / 300)

print(displacements)
print(M_bar1, M_bar2 + moments_state1, M_bar3)
print(M_bar1[1] + moments_state1[0] + M_bar2[0],
      M_bar2[1] + moments_state1[1] + M_bar3[1])


""" 2d frame rigid diaphragm mass, stiffness and Periods computation
"""
import numpy as np

# build basic frame with objects

first floor with 4 columsn + slab + loads + 3 beams 
second floor 4 cols + 4 cols2 + slab + extra loads + 3 beams
def 2Dframe(object):
    def __init__(self, colobj_list, beamobj_list, slabobj_list):
        # col=[[colobj1, colobj2,...]...[ncolobj1,ncolobj2...]], len = DOF
    def massMatrix():

        ..
    def stiffnessMatrix():

    def eigenvalues():

# mass matrix, take mass from frame
# stiffness matrix
# eigenvalues

from scipy.linalg import eig
from numpy import mat, diag, array, pi, sort


def frecuencias_modos(rigideces, masas):
    Masas = diag(masas)
    k1, k2, k3 = rigideces
    Rigideces = array([[k1 + k2, -k2, 0], [-k2, k2 + k3, -k3], [0, -k3, k3]])

    freq2, modos = eig(Rigideces, Masas)

    freq = freq2**.5
    periodos = 2 * pi / freq
    print 'frec. = ', freq
    print 'periodos = ', periodos
    print modos
    return freq, periodos, modos


frecuencias_modos([180e3, 406e3, 406e3], [144., 149., 130.])
