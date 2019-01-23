from stresses_soil import Damy
from numpy import array, ones, pi, exp, log10, zeros

def desplante(peso_medio, presion_critica, peso_volumetrico, naf_desde_arriba=0, peso_vol_agua = 9.81):
    nivel_desplante = (peso_medio - 0.5*peso_vol_agua*naf_desde_arriba)/(1.5*peso_volumetrico-0.5*peso_vol_agua)
    return nivel_desplante

# print desplante(32.0, 14.0, 17.0)

def deformaciones(inp, largo_cim, ancho_cim,  rigideces_suelo, espesores, poisson=[0.5]):
    profundidades = array([sum(espesores[:i]) + 0.5*espesores[i] for i in range(len(espesores))])
    if len(poisson) == 1:
        poisson = ones((len(rigideces_suelo)))*0.5
    else:
        poisson = array(poisson)
    sz, sx, sy = Damy(inp, largo_cim, ancho_cim, profundidades, poisson)
    defs = array(espesores)*(sz - poisson*(sx + sy))/rigideces_suelo
    return defs



# exp_inmediata = deformaciones(51.0, 20.0, 30.6, [4955.0, 4905.0, 5005.0], [1.0, 4.0, 5.0])
# comp_inmediata = deformaciones(32.0, 20.0, 30.6, [3980, 4000, 3890.0], [1.0, 4.0, 5.0])
# print comp_inmediata, sum(comp_inmediata)

def Factor_consolidacion(tiempo_years, coeficientes_consolidacion, espesores, fronteras_permeables=2):
    factor_tiempo = tiempo_years*365.25*86400*coeficientes_consolidacion/(espesores/fronteras_permeables)**2
    factor_consolidacion = 1.0
    for n in range(10):
        N = 0.5*pi*(2*n+1)
        factor_consolidacion -= 2*exp(-N**2*factor_tiempo)/N**2
    return factor_consolidacion, factor_tiempo

def consolidacion(inp, largo_cim, ancho_cim, tiempo_anos, espesores, rigideces_consolidacion_primaria, coeficientes_consolidacion,  poisson=[0.5], fronteras_permeables = 2, rigideces_consolidacion_secundaria=1e9, xis = 0.5):
    profundidades = array([sum(espesores[:i]) + 0.5*espesores[i] for i in range(len(espesores))])
    if len(poisson) == 1:
        poisson = ones((len(rigideces_consolidacion_primaria)))*0.5
        xis = ones((len(rigideces_consolidacion_primaria)))*5.0
        fronteras_permeables = ones((len(rigideces_consolidacion_primaria))) * fronteras_permeables
    else:
        poisson = array(poisson)
        xis = array(xis)
        fronteras_permeables = array(fronteras_permeables)
    Ut = zeros((len(espesores)))
    factores_tiempo = zeros((len(espesores)))
    for ix, espesor in enumerate(espesores):
        Ut[ix], factores_tiempo[ix] = Factor_consolidacion(tiempo_anos, coeficientes_consolidacion[ix], espesor, fronteras_permeables[ix])

    sz, sx, sy = Damy(inp, largo_cim, ancho_cim, profundidades, poisson)
    # U_t, factores_tiempo = factor_consolidacion(50.0, coeficientes_consolidacion, espesores)
    def_primaria = Ut*array(espesores)*sz/rigideces_consolidacion_primaria
    def_secundaria = array(espesores)*sz/rigideces_consolidacion_secundaria*log10(1+xis*factores_tiempo)
    defs_totales = def_primaria + def_secundaria
    hundimiento = sum(defs_totales)

    return {'primaria':def_primaria, 'secundaria': def_secundaria, 'totales':defs_totales, 'hundimiento': hundimiento}


cons = consolidacion(19.0, 20.0, 30.6, 50.0, [1.0, 4.0, 5.0] , [6200.0, 6795.0, 7200.0], [2e-7, 1.2e-7, 1e-7],rigideces_consolidacion_secundaria=[11295.0, 12400.0, 12805.0], fronteras_permeables =[1,2,2])

# print(cons['totales'], cons['hundimiento'])
print('expansion inmediata')
exp_inmediata = deformaciones(3.5*17, 24.0, 32.0, [4890.0, 4840.0, 4940.0], [0.5, 4.0, 5.0])
print exp_inmediata, sum(exp_inmediata)
comp_inmediata = deformaciones(44.5, 24.0, 32.0, [3930, 3950, 3840], [0.5, 4.0, 5.0])
print comp_inmediata, sum(comp_inmediata)
cons2 = consolidacion(21.5, 24, 32.0, 50.0, [0.5, 4.0, 5.0] , [1./7.8e-4, 1./5.2e-4, 1./3.4e-4], [2e-7, 1.2e-7, 1e-7], rigideces_consolidacion_secundaria=[11295.0, 12400.0, 12805.0], fronteras_permeables =[1,2,2])
print(cons2['primaria'], sum(cons2['primaria']))