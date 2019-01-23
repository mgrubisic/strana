from numpy import array, zeros, delete, diag, dot
from scipy.linalg import pinv2
import tridislab as ts


def truss3d(nodes, bars, restraints, areas, elastic_modulus=2.0e6, loads=None, temperature_gradient=0.0, tcl_name=None, nodetags=[], bartags=[], thermal_exp_coeff=12.0e-6, roller_x=[], roller_y=[], roller_z=[]):
    results_dict = {}
    displacements = []
    deformations = []
    internal_forces = []
    continuity = zeros((len(bars), 3 * len(nodes)))
    lengths = []
    removed_dofs = []
    for bix, bar in enumerate(bars):
        xi = nodes[bar[0]][0]
        xj = nodes[bar[1]][0]
        yi = nodes[bar[0]][1]
        yj = nodes[bar[1]][1]
        zi = nodes[bar[0]][2]
        zj = nodes[bar[1]][2]
        length = ((xi - xj)**2 + (yi - yj)**2 + (zi - zj)**2)**.5
        u = float(xj - xi) / length
        v = float(yj - yi) / length
        w = float(zj - zi) / length
        # print u, v, w
        lengths.append(length)
        for nix, node in enumerate(bar):
            if node == 0:
                xcoord = node
                ycoord = node + 1
                zcoord = node + 2
            else:
                xcoord = 3 * node
                ycoord = 3 * node + 1
                zcoord = 3 * node + 2
            if nix == 0:
                continuity[bix, xcoord] = -u
                continuity[bix, ycoord] = -v
                continuity[bix, zcoord] = -w
            else:
                continuity[bix, xcoord] = u
                continuity[bix, ycoord] = v
                continuity[bix, zcoord] = w
    for rest in restraints:
        removed_dofs.append(3 * rest)
        removed_dofs.append(3 * rest + 1)
        removed_dofs.append(3 * rest + 2)
    for roller in roller_x:
        removed_dofs.append(3 * roller)
    for roller in roller_y:
        removed_dofs.append(3 * roller + 1)
    for roller in roller_z:
        removed_dofs.append(3 * roller + 2)
        results_dict['original continuity'] = continuity
    continuity = delete(continuity, removed_dofs, 1)
    k = elastic_modulus * diag([float(areas[j]) / lengths[j]
                                for j in range(len(bars))])
    # return continuity, k
    K = dot(dot(continuity.T, k), continuity)
    if loads is not None:
        if len(loads) != len(K):
            print 'length of loads must be the same as dof {}'.format(len(K))
        displacements = dot(pinv2(K), array(loads))
        deformations = dot(continuity, displacements)
        internal_forces = dot(k, deformations)
    kt = elastic_modulus * temperature_gradient * \
        thermal_exp_coeff * array([float(i) for i in areas])
    # print kt
    thermal_loads = dot(continuity.T, kt)
    thermal_displ = dot(pinv2(K), thermal_loads)
    thermal_defs = dot(continuity, thermal_displ)
    thermal_forces = dot(k, thermal_defs)
    thermal_resultant = thermal_forces - kt
    results_dict['continuity'] = continuity
    results_dict['temp initial forces'] = -kt
    results_dict['K'] = K
    results_dict['u'] = displacements
    results_dict['deformations'] = deformations
    results_dict['forces'] = internal_forces
    results_dict['temp u'] = thermal_displ
    results_dict['temp loads'] = thermal_loads
    results_dict['k'] = k
    results_dict['lengths'] = lengths
    results_dict['temp forces'] = thermal_forces
    results_dict['temp resultant'] = thermal_resultant

    results_dict['lengths'] = lengths

    if tcl_name is not None:
        with open(tcl_name + '.tcl', 'w') as file:
            file.write('#File automatically generated with truss3d.py \n')
            file.write('model BasicBuilder -ndm 3 -ndf 3 \n')

            for tag, node in zip(nodetags, nodes):
                file.write('node {} {} {} {} \n'.format(
                    tag, node[0], node[1], node[2]))
            for restr in restraints:
                file.write('fix {} 1 1 1 \n'.format(restr))
            file.write(
                'uniaxialMaterial Elastic 1 {} \n'.format(elastic_modulus))
            for tag, bar, area in zip(bartags, bars, areas):
                file.write('element truss {} {} {} {} 1 \n'.format(
                    tag, bar[0], bar[1], area))
            file.write('timeSeries Linear 1 \npattern Plain 1 1 { \n')
            loads_dof = [loads[3 * i:3 * i + 3] for i in range(len(K) / 3)]
            for rest in restraints:
                if rest in nodetags:
                    nodetags.remove(rest)
            # print loads_dof, loads
            for tag, load in zip(nodetags, loads_dof):
                file.write('load {} {} {} {} \n'.format(
                    tag, load[0], load[1], load[2]))

            file.write('} \n')
            file.write(
                'system BandGeneral \nnumberer Plain \nconstraints Plain \n')
            file.write(
                'integrator LoadControl 1.0 \nalgorithm Linear \nanalysis Static \n')
            file.write('analyze 1\n')
            if temperature_gradient is not 0:
                file.write('reset \n')
                file.write('timeSeries Linear 2 \npattern Plain 2 1 { \n')
                thermal_dof = [thermal_loads[3 * i:3 * i + 3]
                               for i in range(len(K) / 3)]
                for rest in restraints:
                    if rest in nodetags:
                        nodetags.remove(rest)
                # print loads_dof, loads
                for tag, load in zip(nodetags, thermal_dof):
                    file.write('load {} {} {} {} \n'.format(
                        tag, load[0], load[1], load[2]))

                file.write('} \n')
                file.write(
                    'system BandGeneral \nnumberer Plain \nconstraints Plain \n')
                file.write(
                    'integrator LoadControl 1.0 \nalgorithm Linear \nanalysis Static \n')
                file.write('analyze 1\n')

    return results_dict


crujias_x = [50 * x for x in range(16)]
crujias_y = [50 * y for y in range(26)]
defectos = [584000, 584000, 0] + [-58400, 0, 0] + [0, 0, 0] * \
    (len(crujias_x) - 1) + [0, -584000., 0] + \
    [0, 0, 0] * (791 - 4 - (len(crujias_x) - 1) - 3)
defectos_constructivos = [0, 0, 0] * 17 + [584e3,
                                           0, 0, -584e3, 0, 0] + [0, 0, 0] * (791 - 19 - 4)
grav_sismo = [4.04, 1.12, -17.98]
asentamientos67 = [243.3e3, 243.3e3, -486.7e3]
asentamientos68 = [-243.3e3, 243.3e3, -486.7e3]
asentamientos83 = [243.3e3, -243.3e3, -486.7e3]
asentamientos84 = [-243.3e3, -243.3e3, -486.7e3]
fzas_asentamientos = [0, 0, 0] * 67 + asentamientos67 + asentamientos68 + \
    [0, 0, 0] * 14 + asentamientos83 + \
    asentamientos84 + [0, 0, 0] * (791 - 5 - 84)

apoyos = [[3, 4], [11, 4], [3, 20], [11, 20]]
nodos, barras = ts.tridislab(crujias_x, crujias_y, 50., apoyos, grav_sismo)

r = truss3d(nodos['xy'], barras['ij'], nodos['restraints'], [7.3 for di in range(len(
    barras['ij']))], 2.0e6, fzas_asentamientos, 3402., 'tridilosa-15x25', nodos['tags'], barras['tags'])
# [(i, barra.index(479)) for i, barra in enumerate(b) if 471 in barra]
# >>> rd1['forces'][:2]
# array([ -7.78044296e-10,  -3.73850680e+04])
# >>> rd1['forces'][18]
# -505345.14
# >>> rd1['forces'][20]
# 47392.041
# >>> rd1['forces'][22]
# 25518.120

# >>> rd1['forces'][18]+584000
# 78654.854
