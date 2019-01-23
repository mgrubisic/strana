from numpy import array, sin, cos, sinh, cosh, pi


def periodo_fundamental(espesor, modulo_cortante, peso_unit, amort_critico = 0.1, periodo_onda=2.45, gravedad = 9.81):
    """returns T1 para un solo estrato de
    arcilla de consistencia blanda
    """
    periodo = 4*espesor*(peso_unit/gravedad/modulo_cortante)**.5
    frecuencia = 2*pi/periodo
    eta_gamma_G = 2*amort_critico
    velocidad_onda_corte = (modulo_cortante/(peso_unit/gravedad))**.5
    alpha = espesor*frecuencia/2**.5/velocidad_onda_corte * (((1+eta_gamma_G**2)**.5 - 1)/(1+(eta_gamma_G)**2))**.5
    beta = espesor*frecuencia/2**.5/velocidad_onda_corte * (((1+eta_gamma_G**2)**.5 + 1)/(1+(eta_gamma_G)**2))**.5
    factor_amplificacion = 1/(sinh(alpha)**2*sin(beta)**2 + cosh(alpha)**2*cos(beta)**2)**.5
    for name in dir():
        myvalue = eval(name)
        print(name, '   :   ', myvalue)
    return factor_amplificacion

# print periodo_fundamental(30., 3300, 13.5, 0.14)
print periodo_fundamental(15., 3300, 13.5, 0.14)
# print periodo_fundamental(40., 3300, 13.5, 0.14)
# print periodo_fundamental(50., 3300, 13.5, 0.14)

espesores = array([5.0, 4.0, 4.0])
modulos = array([3200., 3300., 3400.])
pesos = array([12.0, 14.0, 17.0(0.1,pi/2 - 0.1)(0.1,pi/2 - 0.1)])

def periodo_estratos(espesores, modulos_corte, pesos_unit):
    """ los datos entran en listas
    """
    espesores = array(espesores)
    modulos = array(modulos_corte)
    pesos = array(pesos_unit)
    periodo = 4.0/(9.81)**.5*(sum(espesores/modulos)*sum(pesos*espesores*([])))**.5

    # ts = 1.16s
