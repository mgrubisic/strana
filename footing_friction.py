import numpy as np
"""
revisar la seguridad del terreno de cimentacion
por capacidad de carga, de la zapata rectangular,
de concreto reforzado.
"""
ancho_zapata = 1.5
ancho_dado = 0.25
profundidad_desplante = 0.6
profundidad_dado = 0.3
largo_zapata = 1.9
largo_dado = 0.3
carga_axial = 260.  # kN
compacidad_relativa = 0.58
peso_vol_seco = 16.
densidad_solidos = 2.6
ang_fric_interna_campo = 37 * np.pi / 180
factor_resistencia = 0.35  # zona 1
factor_carga = 1.4
factor_carga_suelo = 1.1
peso_vol_concreto = 9.81 * 2.4
peso_vol_agua = 9.81
relacion_vacios = densidad_solidos * peso_vol_agua / peso_vol_seco - 1
peso_saturado = ((densidad_solidos + relacion_vacios) /
                 (1 + relacion_vacios)) * peso_vol_agua
peso_sumergido = peso_saturado - peso_vol_agua
area_zapata = ancho_zapata * largo_zapata
"""
i. NAF = 20 m
ii. NAF = al nivel de la superficie
iii. NAF = 0.5 m bajo el nivel de desplante
"""
naf = 200

if  profundidad_desplante > naf:

    presion_total = (profundidad_desplante - naf) * peso_saturado + naf * peso_vol_seco
    presion_hidraulica = peso_vol_agua * (profundidad_desplante - naf)
    presion_efectiva = presion_total - presion_hidraulica

elif naf > profundidad_desplante and naf > ancho_zapata:

    presion_total = profundidad_desplante * peso_vol_seco
    presion_efectiva = presion_total

elif naf > profundidad_desplante and naf <= ancho_zapata:

    presion_total = profundidad_desplante
    # presion_hidraulica =
    # presion_efectiva =

else:
    print 'error en el naf'
# Vesic , RCDF2004
if compacidad_relativa < 0.67:
    alpha = 0.67 + compacidad_relativa - 0.75 * compacidad_relativa**2
    if alpha > 1:
        alpha = 1
else:
    alpha = 1

phi = np.arctan(alpha * np.tan(ang_fric_interna_campo))
coef_presion_pasiva = np.tan(np.pi / 4 + phi / 2)**2
peso_zapata = peso_vol_concreto * ancho_zapata * \
    largo_zapata * (profundidad_desplante - profundidad_dado)
peso_dado = peso_vol_concreto * ancho_dado * largo_dado * profundidad_dado

"""factores de capacidad de carga y forma"""
Nq = np.exp(np.pi * np.tan(phi)) * coef_presion_pasiva
Ny = 2 * (Nq + 1) * np.tan(phi)
Nc = (Nq - 1) / np.tan(phi)

fc = 1 + 0.25 * ancho_zapata / largo_zapata
fy = 1 - 0.4 * ancho_zapata / largo_zapata
fq = 1 + np.tan(phi) * ancho_zapata / largo_zapata

if naf < profundidad_desplante:
    peso_suelo = peso_sumergido
elif naf > ancho_zapata:
    peso_suelo = peso_vol_seco

peso_relleno = peso_vol_seco


capacidad_carga_ultima = presion_efectiva + factor_resistencia * \
    (presion_total * (Nq * fq - 1) +
     0.5 * peso_suelo * ancho_zapata * Ny * fy)

# sigma Q
carga_actuante_neta = factor_carga*(carga_axial + peso_zapata + peso_dado) + factor_carga_suelo * peso_relleno
presion_ultima = carga_actuante_neta / area_zapata


factor_seguridad_total = capacidad_carga_ultima / presion_ultima

print factor_seguridad_total

# for name in dir():
#     myvalue = eval(name)
#     print(name, '  :  ', myvalue)
