""" calcular las dimensiones del refuerzo para las cintas
de tierra armada para la siguiente configuracion
"""
import numpy as np
altura = 9
phi = 35 * np.pi / 180
peso_vol_suelo = 18  #kN/m3
separacion_vertical = 0.5
separacion_horizontal = 1
delta = 24 * np.pi / 180  # angulo de friccion estatica entre suelo y cinta
esfuerzo_fluencia = 240000  # kPa
tasa_corrosion = 0.025/1000  # m/año
vida_util = 50  # años
factor_seguridad_tension = 3
factor_seguridad_desliz = 3

coeficiente_activo = np.tan(np.pi / 4 - phi / 2)**2
presion_vertical_maxima = altura * peso_vol_suelo
maxima_tension = separacion_horizontal * separacion_vertical * \
    coeficiente_activo * presion_vertical_maxima  # en la base del muro
area = factor_seguridad_tension*coeficiente_activo* presion_vertical_maxima * separacion_vertical * separacion_horizontal/esfuerzo_fluencia
print(maxima_tension)
print(area)
"""proponer dimensiones transversales (rectangular)
"""
ancho = 0.075
peralte = area/ancho
print(peralte)
""" agregamos espesor adicional por la corrosion
"""
peralte += tasa_corrosion * vida_util
print(peralte)

profundidades = np.arange(0.5,altura+0.5,separacion_vertical)
longitud = np.zeros(len(profundidades))
for z in profundidades:
	longitud[z] = (altura - z)/(np.tan(np.pi/4+phi/2))+coeficiente_activo*separacion_horizontal*separacion_vertical*factor_seguridad_desliz/(2*ancho*np.tan(delta))
	print('longitud de la barra a profundidad  '
		+ str(z) + ' m  --  es:  ' + str(longitud[z]) + '  m')

