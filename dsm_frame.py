from numpy import array, zeros, delete, diag, dot
from scipy.linalg import inv
# starts numbering from 0
# use numbers for indexing

def Frame_2D(nodes, bars, restraints, loads):
    continuity = zeros((4*len(bars), 2 * len(nodes)))
    lengths = []
    removed_dofs = []
    for bix, bar in enumerate(bars):
        xi = nodes[bar[0]][0]
        xj = nodes[bar[1]][0]
        yi = nodes[bar[0]][1]
        yj = nodes[bar[1]][1]
        length = ((xi - xj)**2 + (yi - yj)**2)**.5
        u = float(xj - xi) / length
        v = float(yj - yi) / length
        lengths.append(length)
        for nix, node in enumerate(bar):
            if node == 0:
                xcoord = node
                ycoord = node + 1
            else:
                xcoord = 2 * node
                ycoord = 2 * node + 1
            if nix == 0:
                continuity[bix, xcoord] = -u
                continuity[bix, ycoord] = -v
            else:
                continuity[bix, xcoord] = u
                continuity[bix, ycoord] = v
    for rest in restraints:
        removed_dofs.append(2 * rest)
        removed_dofs.append(2 * rest + 1)
    for roller in roller_x:
        removed_dofs.append(2 * rest)
    for roller in roller_y:
        removed_dofs.append(2 * rest + 1)
    continuity = delete(continuity, removed_dofs, 1)
    k = diag([elastic_moduli[i] * areas[i] / lengths[i]
              for i in range(len(bars))])
    K = dot(dot(continuity.T, k), continuity)
    displacements = dot(inv(K), array(loads))
    deformations = dot(continuity, displacements)
    internal_forces = dot(k, deformations)
    thermal_loads = dot(continuity.T, [elastic_moduli[
                        i] * areas[i] * thermal_exp_coeff * temperature_gradient for i in range(len(bars))])

    return displacements, deformations, internal_forces, continuity, thermal_loads


nodes = [[0, 0], [4.0, 0], [0, 3.0], [4.0, 3.0]]
bars = [[0, 1], [0, 2], [2, 3], [0, 3], [1, 2]]
restraints = [1, 3]
rollers_x = []
rollers_y = []
areas = [5e-4 for i in range(len(bars))]
elastic_moduli = [2e7 for i in range(len(bars))]
loads = [-10800, -9600, -10800, 9600]

c = truss2d(nodes, bars, restraints, areas, elastic_moduli, loads, 50.)

print c[4]

n2 = [[0, 0], [3.0, 0], [6.0, 0], [9.0, 0], [3.0, 4.0], [6.0, 8.0], [9.0, 4.0]]
b2 = [[0, 1], [1, 2], [2, 3], [0, 4], [4, 5], [
    5, 6], [6, 3], [1, 4], [4, 2], [2, 5], [2, 6]]
elastic_moduli = [2e6 for i in range(len(b2))]
areas = [0.05 for i in range(len(b2))]
l2 = array([-0., -6600., 6600., -6600., -3960.,
            11880., -0., 17160., 7920., 6600.])

r2 = [0, 2]

c2 = truss2d(n2, b2, r2, areas, elastic_moduli, l2)
