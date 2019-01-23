# TODO: incorporate into submodule of Plate

def placa_tension(ancho, espesor, gramil, separaciones, diametro_tornillo =  2.54, fy_placa = 3520.0, fu_placa=4570.0, fu_tornillos = 8440.0, agujeros_estandar = 0.2, planos_corte = 2, factor_area_U = 0.9):
    """meter la placa critica y separaciones como [5.0, 10.0, 10.0] para tres tornillos
    se consideran 2mm adicional en agujeros
    ancho es el ancho total de placa
    """
    area_tornillo = Pi*diametro_tornillo**2/4
    tornillos = 2 * len(separaciones)
    aplastamiento = tornillos * 0.75 * 2.4 * diametro_tornillo * espesor * fu_placa
    corte_tornillos = tornillos * planos_corte * area_tornillo * 0.6 * fu_tornillos
    flujo_plastico = 0.9 * ancho * espesor * fy_placa
    area_neta = ancho * espesor - 2*(diametro_tornillo + agujeros_estandar) * espesor
    fractura_neta = 0.75 * area_neta * factor_area_U * fu_placa
    print 'ATc ------ Anc ---- ATt ----- Ant ----- 0.6Anc ----- Rrup'
    for index, torn in enumerate(separaciones):
        seps = separaciones[:index+1]
        torns = 2 * len(seps)
        total_corte = 2*sum(seps) * espesor
        total_tension = gramil * espesor
        neta_corte = espesor*(2*sum(seps) - diametro_tornillo*(torns-1))
        neta_tension = (gramil - diametro_tornillo) * espesor
        if neta_tension >= 0.6 * neta_corte:
            Rrup = tornillos/torns* 0.75 * (0.6 * fy_placa * total_corte + fu_placa* neta_tension)
        else:
            Rrup = tornillos/torns*0.75 * (0.6*fu_placa*neta_corte + fy_placa * total_tension)
        areas_ruptura = array([total_corte, neta_corte, total_tension, neta_tension, 0.6 * neta_corte])
        print areas_ruptura, Rrup * 1e-3

def eslabones(cortantes, altura_entrepiso, longitud_crujia, bf = 40, tf = 2.54, depth=40.0, tw=1.0,  tol = 1.0):
    """regresa una lista con las secciones de los eslabones y separaciones de atiesadores, con la fuerza resultante """
    for cortante in cortantes:
        V = cortante*altura_entrepiso/longitud_crujia
        eslab = profile('esl')
        eslab.beamcolumn_W(depth, tw, bf, tf)
        while abs(V - eslab.vp*1e-3 ) >= tol:
            if V>eslab.vp*1e-3:
                tw = tw*1.15
                depth = depth*1.15
                # bf = bf * 1.15
                # tf = tf * 1.15
                eslab.beamcolumn_W(depth, tw, bf, tf)
            else:
                tw = tw*0.95
                depth = depth*0.95
                # bf = bf * 0.95
                # tf = tf *0.95
                eslab.beamcolumn_W(depth, tw, bf, tf)
        eslabon = 1.6*eslab.mp/eslab.vp
        fza_col = 1.25*1.1*sin((longitud_crujia - eslabon*1e-3)/2./altura_entrepiso)*V
        separacion = 30*tw - depth/5
        print '|{:.1f}|{:.2f}|{:.1f}|{:.0f}|{:.1f}|'.format(fza_col, tw, depth, eslabon, separacion)

# print '---'
# print '| $P$  | $t_a$ |  $d$ | $e$ | $s$ |'
# print '| T  | cm |  cm | cm | cm |'
# eslabones([608/4, 557./4, 498.0/4, 410.0/4, 321./4, 175./4], 3.6, 8.0, 40.0, 2.54, 40.0, 0.85)
# eslabones([608/4, 557./4, 498.0/4, 410.0/4, 321./4, 175./4], 3.6, 12, 40.0, 2.54, 40.0, 0.85)

# import print_options
# print_options.set_float_precision(4)
