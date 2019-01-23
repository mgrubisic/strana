""" Primero simple como DEM
TODO:
esta mal programada la alpha pero sirve
"""
from numpy import pi, array


def pilotes_friccion(diametro, espesores, cohesiones, presiones_efectivas, FR=1.0):
    perimetro = pi*diametro
    cohesiones = array(cohesiones)
    espesores = array(espesores)
    alpha = array([0.5*(presiones_efectivas/cohesiones)**.5])
    print 'coeficiente de aherencia', alpha
    friccion_negativa = perimetro*FR*sum(alpha*cohesiones*espesores)

    return friccion_negativa

fn= pilotes_friccion(0.4, [2.5, 3.0, 3.6], [22.0, 22.0, 22.0], [17.5, 41.29, 55.11])

print fn
print sum(fn)
