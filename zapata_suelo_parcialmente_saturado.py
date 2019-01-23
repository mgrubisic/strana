""" zapata corrida en suelo parcialmente saturado
"""
import numpy as np
ancho_estrato = 6
peso_vol = 16.3
ancho_zapata = 1.5
largo_zapata = 1.5
peso_suelo = peso_vol * ancho_estrato
cohesion_no_drenada = 11.3
# ancho_dado = 0.25
profundidad_zapata = 0.8
# profundidad_dado = 0.3
# largo_dado = 0.3
# carga_axial = 260  # kN
# compacidad_relativa = 0.58
# peso_vol_seco = 16
# densidad_solidos = 2.6
ang_fric_interna_campo = 29 * np.pi / 180
ang_fric_fredlung = 16.5 * np.pi / 180
factor_resistencia = 0.45  # zona1
phi = ang_fric_interna_campo
succion = 50
# factor_carga = 1.4
# peso_vol_concreto = 9.81 * 2.4
# # peso_vol_agua = 9.81
# relacion_vacios = densidad_solidos * peso_vol_agua / peso_vol_seco - 1
# peso_saturado = ((densidad_solidos + relacion_vacios) /
#                  (1 + relacion_vacios)) * peso_vol_agua
# peso_sumergido = peso_saturado - peso_vol_agua
# area_zapata = ancho_zapata * largo_zapata

# phi = np.arctan(alpha * np.tan(ang_fric_interna_campo))
presion_pasiva = np.tan(np.pi / 4 + phi / 2)**2
cohesion_total = cohesion_no_drenada + succion * np.tan(ang_fric_fredlung)
# peso_zapata = peso_vol_concreto * ancho_zapata * \
#     largo_zapata * (profundidad_zapata - profundidad_dado)
# peso_dado = peso_vol_concreto * ancho_dado * largo_dado * profundidad_dado
presion_desplante = profundidad_zapata * peso_suelo
"""factores de capacidad de carga y forma"""

Nq = np.exp(np.pi * np.tan(phi)) * presion_pasiva
Ny = 2 * (Nq + 1) * np.tan(phi)
Nc = (Nq - 1) / np.tan(phi)

fc = 1 + 0.25 * ancho_zapata / largo_zapata
fy = 1 - 0.4 * ancho_zapata / largo_zapata
fq = 1 + np.tan(phi) * ancho_zapata / largo_zapata


capacidad_carga_ultima = factor_resistencia * \
    (cohesion_total * Nc * fc + presion_desplante * Nq * fq +
     0.5 * peso_suelo * ancho_zapata * Ny * fy)
# presion_ultima = presion_media * factor_carga

# factor_seguridad_total = capacidad_carga_ultima / presion_ultima
for name in dir():
    myvalue = eval(name)
    print(name, '  :  ', myvalue)
