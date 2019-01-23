from numpy import pi, tan, sin, cos, deg2rad, rad2deg, arange
"""data
"""

def talud_simplificado(altura, base, peso_unit, cohesion_no_drenada, ang_friccion_interna, sobrecarga, nivel_freatico, coeff_sismico_vertical = 0, coeff_sismico_horizontal = 0):
# seco aguas abajo
    phi = deg2rad(ang_friccion_interna)
    for angulo in range(20,70):
        angulo = deg2rad(angulo)

        peso_a = altura**2/tan(angulo)/2 * peso_unit
        sobrecarga_a = sobrecarga * altura /tan(angulo)
        fza_sismica_horizontal_a = coeff_sismico_horizontal * peso_a
        fza_sismica_vertical_a = coeff_sismico_vertical * peso_a
        subpresion_a = nivel_freatico**2 * 9.81/2/sin(angulo)
        cohesion_a = altura/sin(angulo) * cohesion_no_drenada


        peso_b = peso_unit * altura * float(base) /2
        fza_sismica_horizontal_b = coeff_sismico_horizontal * peso_b
        fza_sismica_vertical_b = coeff_sismico_vertical * peso_b
        cohesion_b = base * cohesion_no_drenada
        subpresion_b = base * nivel_freatico * 9.81/2


        Fzas_resistentes = tan(phi)*((peso_a + sobrecarga_a + fza_sismica_vertical_a)*cos(angulo) - fza_sismica_horizontal_a*sin(angulo) - subpresion_a) + cohesion_a + (peso_b + fza_sismica_vertical_b - subpresion_b)*tan(phi) + cohesion_b
        Fzas_actuantes = fza_sismica_horizontal_b + fza_sismica_vertical_b*cos(angulo) + sin(angulo)*(sobrecarga_a + peso_a + fza_sismica_vertical_a)
        factor_seg = Fzas_resistentes / Fzas_actuantes
        print rad2deg(angulo), factor_seg


# talud_simplificado(8., 16., 17., 15., 25., 50., 6., 0, 0.1)
for x in range(1000):
    talud_simplificado(x, 16., 17., 15., 25., 50., 6., 0, 0.1)
