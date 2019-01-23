from numpy import genfromtxt, sign
from easyplot import easyplot


def return_map(strain_list, yield_stress = 2530.0, Young=2.039e6, Tangent = 4.06e5):
    strains = genfromtxt(strain_list)
    stresses = [0]
    alpha = [0]
    plastic = [0]
    for n, strain in enumerate(strains[1:]):
        trial = Young*(strain - plastic[n])
        surface = abs(trial) - yield_stress - Tangent*alpha[n]
        if surface <= 0:
            plastic.append(plastic[n])
            alpha.append(alpha[n])
            stresses.append(trial)
        else:
            strain_rate = surface/(Young + Tangent)
            stresses.append( (1-strain_rate*Young/abs(trial))*trial )
            plastic.append( plastic[n]+strain_rate*sign(trial) )
            alpha.append( alpha[n] + strain_rate )
    print len(strains), len(stresses)
    easyplot(strains, stresses, 'strains', 'stresses', 'return mapping algorithm', 'plastic.test.1')
    return

return_map('prueba.perez.gavilan.csv')
