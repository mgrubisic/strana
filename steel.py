#-*- coding: utf-8 -*-
from __future__ import division
"""TODO [x1, x2, x3] for dims with indicial representation of quantities
notation and general equations with are simple and work in all dimensions.. what does this mean?
TODO importar perfiles del IMCA, ASCI con propiedades.
TODO que pasa si varias placas tienen diferentes E, fy? esto es un ejercicio academico interesante
"""
from numpy import array, diff, sign, sin, interp, set_printoptions, dot
from numpy import sum as sum
from numpy import pi as Pi
set_printoptions(precision=4)
from pandas import read_csv, DataFrame
import pprint
pp = pprint.PrettyPrinter(indent=4)

def correccion_inelasticidad(esf, fy):
    return fy*(1-fy/4/esf) if esf > fy/2 else esf

def Tinertia(x, xt, y, yt):
    """ para dos bloques uno sobre otro en forma de t invertida
        yt
        |
     y  |
    ____|____ xt
    -> x
    """
    area = x*xt+yt*y
    centroide = (x*xt**2/2+(xt+0.5*y)*yt*y)/area
    Iy = yt*y**3/12
    Ix = x*xt**3/12
    I1 = Ix + x*xt*abs(centroide - xt/2)**2 + Iy + yt*y*abs(xt+0.5*y - centroide)**2
    I2 = yt**3*y/12 + x**3*xt/12
    return centroide, area, I1, I2


def Ninertia(*objects):
    """ returns centroid and inertia of n objects
    [C1, A1, I1], [C2, A2, I2] [...]
    """
    area = sum([objects[i][1] for i in range(len(objects))])
    for obj in objects:
        if len(obj) != 3:
            print 'please enter objects in this form [C1, A1, I1]..'
            return
    centroid = sum([objects[i][1]*objects[i][0] for i in range(len(objects))])/area
    inertia = sum([objects[i][2]+objects[i][1]*(objects[i][0] - centroid)**2 for i in range(len(objects))])
    S1 = inertia/centroid
    return centroid, inertia, S1

# print Tinertia(15*0.266, 0.266, 2.03, 0.266)

class Placa(object):
    def __init__(self, fy=2530.0, espesor=1.0, x2=1.0, x1=1e12, bordes_apoyados=2, bordes_empotrados=0, Young=2.039e6, Poisson=0.30):
        """ en placas largas  x1=inf
        """
        factores_placa = {}
        relaciones_aspecto = [0]
        coef_flexion_pura = [0]
        self.fy = fy
        self.Young = Young
        self.Poisson = Poisson
        self.nu = Poisson
        self.E = Young
        self.espesor = float(espesor)
        self.x1 = x1
        self.x2 = x2
        self.relacion_aspecto = self.x1/self.x2
        self.esbeltez = self.espesor/self.x2
        self.Ix = x2*espesor**3/12
        self.Sx = self.Ix/(0.5*espesor)
        self.Iy = espesor*x2**3/12
        self.Sy = self.Iy/(0.5*x2)
        self.area = espesor*x2
        self.area_efectiva = self.area
        self.factor_esbeltez = (self.Young/self.fy)**.5
        if self.relacion_aspecto >= 1.0:
            self.k_cortante = 5.34 + 4.0/self.relacion_aspecto**2
        else:
            self.k_cortante = 4.0 + 5.34/self.relacion_aspecto**2
        # self.k_flexion = interp(self.relacion_aspecto, relaciones_aspecto, coef_flexion_pura)
        self.k_compresion = 1.0
        self.k_atiesador_trabe_armada = 0.75
        if (bordes_apoyados+bordes_empotrados) > 4:
            print('revise las condiciones de frontera')
            return
        elif bordes_apoyados == 2:
            self.k_compresion = 4.0
            relaciones_aspecto = [0.4, 0.5, 0.6, 0.666, 0.75, 0.8, 1.0, 1.5]
            coef_flexion_pura = [29.1, 25.6, 24.1, 23.9,24.1, 24.4, 25.6, 24.1]
        elif bordes_apoyados == 1 and bordes_empotrados == 1:
            relaciones_aspecto = [0.4, 0.5, 0.6, 0.65, 0.66, 0.7, 0.8, 0.9, 1.0]
            coef_flexion_pura = [29.5, 26.0, 24.65, 24.48, 24.48, 24.6, 25.3, 26.6, 28.3]
            # self.k_flexion = interp(self.relacion_aspecto, relaciones_aspecto, coef_flexion_pura)
            self.k_compresion = 5.42
        elif bordes_apoyados == 0 and bordes_empotrados == 2:
            self.k_compresion = 6.97
            relaciones_aspecto = [0.3, 0.35, 0.4, 0.45, 0.47, 0.48, 0.50, 0.60, 0.70]
            coef_flexion_pura = [47.3, 43.0, 40.7, 39.7, 39.6, 39.6, 39.7, 41.8, 45.8]
        elif bordes_apoyados == 1 and bordes_empotrados == 0:
            self.k_compresion = 0.425
        elif bordes_apoyados == 0 and bordes_empotrados == 1:
            self.k_compresion = 1.277
        elif bordes_apoyados == 0 and bordes_empotrados == 4:
            relaciones_aspecto2 = [0.4, 0.5, 0.6, 0.75, 0.8, 1.0, 1.5]
            coef_compresion = [8.4, 6.3, 5.2, 4.3, 4.2, 4.0, 4.3]
            self.k_compresion = interp(self.relacion_aspecto, relaciones_aspecto2, coef_compresion)

        try:
            self.k_flexion = interp(self.relacion_aspecto, relaciones_aspecto, coef_flexion_pura)
            factores_placa['esbeltez'] = self.esbeltez
            factores_placa['relacion_aspecto'] = self.relacion_aspecto
            factores_placa['k_compresion'] = self.k_compresion
            factores_placa['k_cortante'] = self.k_cortante
            factores_placa['k_flexion'] = self.k_flexion
            factores_placa['k_atiesador_trabe_armada'] = self.k_atiesador_trabe_armada
        except UnboundLocalError as e:
                print('algun factor no esta asignado')
                print(e)
        self.factores_placa = factores_placa
        return

    def Resistencia(self, atiesadores = 0, espesor_atiesador=0.63, relacion_areas_atiesadores = 0.1, posicion='1/2', carga_axial_critica=0, atiesado = False):
        """ 1/4" espesor, delta=0.1, atiesado = True; para la resistencia pospandeo
        TODO limitaciones de las relaciones de aspecto para atiesadores que no estan a la mitad
        TODO resistencia a flexion con atiesadores
        """
        resistencias = {}
        atiesador = {}
        self.esf_unitario = self.esbeltez**2 * Pi**2 * self.Young /(12*(1 - self.Poisson**2))
        self.cortante_fluencia = self.fy / 3**0.5
        self.esf_critico_compresion = self.k_compresion * self.esf_unitario
        self.cortante_critico = self.k_cortante * self.esf_unitario
        self.esf_critico_flexion = self.k_flexion * self.esf_unitario
        self.esf_ultimo_compresion = correccion_inelasticidad(self.esf_critico_compresion, self.fy)
        self.esf_ultimo_flexion = correccion_inelasticidad(self.esf_critico_flexion, self.fy)
        self.cortante_ultimo = (0.8*self.cortante_fluencia*self.cortante_critico)**.5 if self.cortante_critico > self.cortante_fluencia/2 else self.cortante_critico
        self.compresion_trabe_armada = correccion_inelasticidad(self.esf_unitario*self.k_atiesador_trabe_armada, self.fy)

        k_pospandeo = 0.43*4 if atiesado == False else 4.0
        # 4* algebra to make eq for flanges work
        factor_resistencia_pospandeo = 1.052/self.esbeltez/k_pospandeo**.5*(carga_axial_critica/self.Young)**.5
        self.ancho_efectivo = self.x2*(1- 0.22/factor_resistencia_pospandeo)/factor_resistencia_pospandeo if factor_resistencia_pospandeo > 0.673 else self.x2
        print 'la = {:.3f}, be = {:.2f}'.format(factor_resistencia_pospandeo, self.ancho_efectivo)
        self.area_efectiva = self.ancho_efectivo * self.espesor
        self.resistencia_pospandeo = self.area_efectiva*self.fy

        self.Ix = self.espesor**3*self.ancho_efectivo/12
        self.Iy = self.ancho_efectivo**3*self.espesor/12

        if atiesadores == 1:
            largo_atiesador = relacion_areas_atiesadores*self.espesor*self.x2/espesor_atiesador
            if posicion == '1/2':
                relacion_rigideces = 24.4 + 112*relacion_areas_atiesadores*(1+relacion_areas_atiesadores)
                inercia_min = 0.092*self.espesor**3*self.x2*relacion_rigideces
                inercia_atiesador = Tinertia(30*self.espesor, self.espesor, largo_atiesador, espesor_atiesador)[2]
                esf_atiesador_elastico = Pi**2*self.Young/(12*(1-self.Poisson**2))*0.425*(espesor_atiesador/largo_atiesador)**2
                esf_atiesador = correccion_inelasticidad(esf_atiesador_elastico, self.fy)
                self.esf_critico_compresion = 4*self.esf_critico_compresion
                self.esf_ultimo_compresion = correccion_inelasticidad(self.esf_critico_compresion, self.fy)
            elif posicion == '1/5':
                self.esf_critico_flexion = self.esf_unitario*129
                self.esf_ultimo_flexion = correccion_inelasticidad(self.esf_critico_flexion, self.fy)
                self.factores_placa['k_flexion'] = 129
            else:
                relaciones_rigidez ={'1/3':(16+200*relacion_areas_atiesadores), '1/4':(16+200*relacion_areas_atiesadores), '1/5':(50+200*relacion_areas_atiesadores)}
                res_posterior_pandeo = {'1/3': 4, '1/4': 6,'1/5': 7}

        elif atiesadores == 2:
            self.esf_critico_compresion = 9*self.esf_critico_compresion
            self.esf_ultimo_compresion = correccion_inelasticidad(self.esf_critico_compresion, self.fy)
            largo_atiesador = relacion_areas_atiesadores*self.espesor*self.x2/espesor_atiesador
            relacion_rigideces = 96.0 + 610.0*relacion_areas_atiesadores + 975.0*relacion_areas_atiesadores**2
            inercia_min = 0.092*self.espesor**3*self.x2*relacion_rigideces
            inercia_atiesador = Tinertia(30*self.espesor, self.espesor, largo_atiesador, espesor_atiesador)[2]
            inercia_min = 0.092*self.espesor**3*self.x2*relacion_rigideces
            inercia_atiesador = Tinertia(30*self.espesor, self.espesor, largo_atiesador, espesor_atiesador)[2]
            while inercia_min > inercia_atiesador:
                print 'el atiesador no tiene inercia suficiente {:.1f} > {:.1f}'.format(inercia_min, inercia_atiesador)
                relacion_areas_atiesadores += 0.015
                largo_atiesador = relacion_areas_atiesadores*self.espesor*self.x2/espesor_atiesador
                relacion_rigideces = 96.0 + 610.0*relacion_areas_atiesadores + 975.0*relacion_areas_atiesadores**2
                inercia_min = 0.092*self.espesor**3*self.x2*relacion_rigideces
                inercia_atiesador = Tinertia(30*self.espesor, self.espesor, largo_atiesador, espesor_atiesador)[2]
            esf_atiesador_elastico = Pi**2*self.Young/(12*(1-self.Poisson**2))*0.425*(espesor_atiesador/largo_atiesador)**2
            esf_atiesador = correccion_inelasticidad(esf_atiesador_elastico, self.fy)

        elif atiesadores > 2:
            print('no se tiene programado mas de 2 atiesadores')
        try:
            atiesador['largo_atiesador'] = largo_atiesador
            atiesador['esfuerzo_ultimo_atiesador'] = esf_atiesador
            atiesador['esfuerzo_elastico_atiesador'] = esf_atiesador_elastico
            atiesador['espesor'] = espesor_atiesador
            atiesador['rel_rigideces'] = relacion_rigideces
            atiesador['inercia']= inercia_atiesador
            atiesador['inercia_min'] = inercia_min
            atiesador['rel_areas'] = relacion_areas_atiesadores
            atiesador['area_atiesador'] = largo_atiesador * espesor_atiesador

        except UnboundLocalError as e:
            print('no hay atiesadores asignados')
            print(e)

        resistencias['esf_unitario'] = self.esf_unitario
        resistencias['compresion_critico'] = self.esf_critico_compresion
        resistencias['compresion_ultimo'] = self.esf_ultimo_compresion
        resistencias['compresion_atiesador_trabe_armada'] = self.compresion_trabe_armada
        resistencias['cortante_critico'] = self.cortante_critico
        resistencias['cortante_ultimo'] = self.cortante_ultimo
        resistencias['cortante_fluencia'] = self.cortante_fluencia
        resistencias['flexion_critico'] =  self.esf_critico_flexion
        resistencias['flexion_ultimo'] =  self.esf_ultimo_flexion
        resistencias['carga_ultima_compresion'] = self.esf_ultimo_compresion * self.area
        resistencias['factor_pospandeo'] = factor_resistencia_pospandeo
        resistencias['ancho_efectivo'] = self.ancho_efectivo
        resistencias['area_efectiva'] = self.area_efectiva
        resistencias['resistencia_pospandeo'] = self.resistencia_pospandeo
        self.resistencias = resistencias
        self.atiesador = atiesador


class Profile(object):
    """TODO: general geometric shape (fem?)
    """
    def __init__(self, name='default profile', fy=3520.0, area=1e-16, Zx=0., Zy=0., Young=2.036e6, shear_modulus=784e3):
        self.name = name
        self.area = area
        self.A = area
        self.area_efectiva = area
        self.fy = fy
        self.E = Young
        self.Young = Young
        self.G = shear_modulus
        self.shear_modulus = shear_modulus
        self.nu = Young / shear_modulus / 2 - 1
        self.Poisson = self.nu
        self.Zx = Zx
        self.Zy = Zy
        self.Py = area * fy
        self.Mp = self.fy*array([Zx, Zy])

    def Geometry(self, inertia_x=0, inertia_y=0, polar_inertia=0, warping_constant=0, ymax=0, xmax=0, hor_dist_shear_center=0, vertical_dist_shear_center=0):
        """ x,y max are measured from centroid
        """
        self.Ix = inertia_x
        self.Iy = inertia_y
        self.flexural_stiffness1 = inertia_x * self.Young
        self.polar_inertia = polar_inertia
        self.J = polar_inertia
        self.Cw = warping_constant
        self.X0 = hor_dist_shear_center
        self.Y0 = vertical_dist_shear_center
        self.gyration_radius_x = (inertia_x / self.area)**.5
        self.rx = (inertia_x / self.area)**.5
        self.gyration_radius_y = (inertia_y / self.area)**.5
        self.ry = (inertia_y / self.area)**.5
        self.r0 = (self.X0**2 + self.Y0**2 + self.rx**2 + self.ry**2)**0.5
        self.H = 1 - (self.X0**2 + self.Y0**2) / self.r0**2
        try:
            self.Sx = inertia_x/ymax
            self.Sy = inertia_y/xmax
            self.My = self.fy*array([self.Sx, self.Sy])
        except ZeroDivisionError as z:
            print z
            print 'ymax xmax no especificado'

    def Plate_girder(self, web, top_flange, bottom_flange, atiesador, separacion_atiesadores=1e12, length=100.0, uniform_load=1.0, shear=0, moment=0, filete=0.6, longitud_apoyo_o_placa=0, ranura_detalle = 7.0, altura_filete=0, FR=0.9):
        """
        TODO limit state summary
        TODO unequal flanges section modulus
        TODO warping constant correct value
        TODO unequal flanges plastic modulus
        """
        self.fy = web.fy
        self.plate_girder_results = {}
        self.web = web
        self.top_flange = top_flange
        self.bottom_flange = bottom_flange
        self.ta = web.espesor
        self.peralte_alma = web.x2
        self.tf1 = top_flange.espesor
        self.tf2 = bottom_flange.espesor
        self.bf1 = top_flange.x2
        self.bf2 = bottom_flange.x2
        self.areas = array([bottom_flange.area, web.area, top_flange.area])
        self.area = sum(self.areas)
        self.altura = self.tf1+self.tf2+self.peralte_alma
        self.centroids = array([self.tf2/2, self.tf2+self.peralte_alma/2,self.tf2+self.peralte_alma+self.tf1/2])
        self.centroid = dot(self.centroids,self.areas)/self.area
        self.inertias1 = array([self.tf2**3*self.bf2/12, self.peralte_alma**3*self.ta/12, self.tf1**3*self.bf1/12])
        self.inertias2 = array([self.bf2**3*self.tf2/12, self.peralte_alma*self.ta**3/12, self.bf1**3*self.tf1/12])
        self.to_centroid = self.centroid - self.centroids
        self.Ix = sum(self.inertias1 + self.areas*self.to_centroid**2)
        self.Iy = sum(self.inertias2)
        self.rx = (self.Ix/self.area)**.5
        self.ry = (self.Iy/self.area)**.5
        # self.J = self.Ix + self.Iy
        self.J =((self.bf1+self.bf2)*self.tf1**3 + self.peralte_alma*self.ta**3)/3
        self.Cw = (self.tf1*self.peralte_alma**2)/12*(self.bf1**3*self.bf2**3/(self.bf1**3+self.bf2**3))
        self.Sx = self.Ix/(self.peralte_alma/2 + self.tf1)
        self.Sy = self.Iy/(self.bf1/2)
        self.Zx = self.bf1*self.tf1*(self.peralte_alma-self.tf1) + 0.25*self.ta*(self.peralte_alma-2*self.tf1)**2
        self.Zy = 0.5*self.bf1**2*self.tf1 + 0.25*self.ta**2*(self.peralte_alma-2*self.tf1)
        self.My = self.fy*array([self.Sx, self.Sy])
        self.Mp = self.fy*array([self.Zx, self.Zy])
        self.relacion_aspecto = separacion_atiesadores/self.peralte_alma
        self.k_trabe_armada = 5.00 + 5.00/(self.relacion_aspecto)**2
        self.factor_kv = (self.k_trabe_armada*web.Young/web.fy)**.5
        if 1./web.esbeltez <= 1.1*self.factor_kv:
            Cv = 1.0
            print '7.2.3'
        elif 1./web.esbeltez > 1.37*self.factor_kv:
            Cv = 1.51*(self.k_trabe_armada*web.Young/web.fy)*web.esbeltez**2
            print '7.2.5'
        else:
            Cv = 1.1*(self.k_trabe_armada*web.Young/web.fy)**.5*web.esbeltez
            print '7.2.4'
        self.VR = array([0,0])
        if 1.0/web.esbeltez<=1.1*self.factor_kv:
            self.VR[1] = 0.6*web.area*web.fy
        else:
            self.VR[1] = 0.6*web.area*web.fy*(Cv+(1-Cv)/(1.15*(1+(self.relacion_aspecto)**2)**.5))
        self.VR[0] = 0.6*web.area*web.fy*Cv
        self.Cv = Cv
        print 'kv, Cv', self.k_trabe_armada, self.Cv
        self.peso = 7850e-4*self.area
        self.flujo_plastico_alma = (5*(top_flange.espesor+filete)+longitud_apoyo_o_placa)*web.fy*web.espesor
        self.abollamiento_alma = 0.8*web.espesor**2*(web.Young*web.fy*top_flange.espesor/web.espesor)**.5*(1+3.0*longitud_apoyo_o_placa*(web.espesor/top_flange.espesor)**1.5/web.x2)
        self.desplazamiento_patin = 67.5e6*web.espesor**3*top_flange.espesor/web.x2**2*(1+0.4*((web.x2/web.espesor)/(length/top_flange.x2))**3)

        area_accion_columna = 2*atiesador.area + 12*web.espesor**2
        I_accion_columna = 12*web.espesor**4/12 + 2*(atiesador.Iy + atiesador.area*(atiesador.x2/2 + web.espesor/2)**2)
        print 'resultados accion columna Iy, A, Ix'
        print atiesador.Iy, area_accion_columna, I_accion_columna

        self.col_atiesador = Profile('col_atiesador', atiesador.fy, area_accion_columna)
        self.col_atiesador.Geometry(I_accion_columna, I_accion_columna, 2*I_accion_columna)
        self.col_atiesador.Buckling(self.peralte_alma*3/4)
        self.col_atiesador.aplastamiento = 1.8*0.75*atiesador.fy*(atiesador.x2 - ranura_detalle)*atiesador.espesor

        Iat = atiesador.espesor*(web.espesor+2*atiesador.x2)**3/12
        j = 2.5/(separacion_atiesadores/web.x2)**2 - 2 if 2.5/(separacion_atiesadores/web.x2)**2 - 2 >= 0.5 else 0.5
        Iat1 = min(separacion_atiesadores, web.x2)*web.espesor**3*j
        Iat2 = web.x2**4*(atiesador.fy/self.Young)**1.5/40
        Imin = 0
        try:
            Imin = Iat1 + (Iat2 - Iat1)*(shear-FR*self.VR[0])/(FR*self.VR[1]-FR*self.VR[0])
            Imin = 0 if Imin < 0 else Imin
            if shear==0 and moment==0:
                shear = length*uniform_load/2
                moment = uniform_load*length**2/8
        except ZeroDivisionError as e:
            print e, 'no hay atiesadores asignados'

        self.plate_girder_results['Actions'] = {'Mu':moment*1e-5, 'Vu':shear*1e-3}
        self.plate_girder_results['Cv'] = self.Cv
        self.plate_girder_results['Mp'] = self.Mp*1e-5
        self.plate_girder_results['My'] = self.My*1e-5
        self.plate_girder_results['VR'] = self.VR*1e-3
        self.plate_girder_results['pandeo'] = 'rango inelastico' if 1./web.esbeltez > 1.00*self.factor_kv else 'rango elastico'
        self.plate_girder_results['kv'] = self.factor_kv
        self.plate_girder_results['relacion_aspecto'] = self.relacion_aspecto
        self.plate_girder_results['k_plate_girder'] = self.k_trabe_armada
        self.plate_girder_results['atiesadores'] = {'Inercia_atiesador':Iat, 'Inercia minima requerida':Imin, 'Adecuado rigideces':Iat>=Imin, 'Iat1':Iat1, 'Iat2':Iat2, 'Pandeo atiesador':1./atiesador.esbeltez>=0.56*(atiesador.Young/atiesador.fy)**.5, 'Resistencia del atiesador':((2*atiesador.area*atiesador.compresion_trabe_armada)*1e-3,'ton'), 'Accion columna del atiesador (t)':self.col_atiesador.buckling_results, 'Aplastamiento del atiesador (t)':self.col_atiesador.aplastamiento}
        self.properties = {'rx, ry':array([self.rx, self.ry]), 'Ix, Iy, J, Cw':array([self.Ix, self.Iy, self.J, self.Cw]), 'Sx, Sy (Myx, Myy)':(array([self.Sx, self.Sy,]), self.My), 'Zx, Zy (Mpx Mpy)':(array([self.Zx, self.Zy]), self.Mp), "peso kg/m" :self.peso}
        self.plate_girder_results['estados_limite'] = {'flujo_plastico_alma':self.flujo_plastico_alma, 'abollamiento_alma':self.abollamiento_alma,'desplazamiento_patin': self.desplazamiento_patin, 'deflexion':(5*uniform_load*length**4/384/self.Young/self.Ix,length*1.0/240)}

    def HSS(self, alma_izq, alma_der, patin_arriba, patin_abajo, cortante=0, momento=0, separacion_atiesadores=1e9):
        """ TODO: unequal flanges or webs correct geometric properties
        TODO: FLEXURE for HSS, correct alma, patines etc
        """
        self.placas = [alma_der, alma_izq, patin_arriba, patin_abajo]
        for placa in self.placas:
            self.area_efectiva += placa.area
        self.area = self.area_efectiva
        d = alma_izq.x2
        b = patin_arriba.x2
        if d != alma_der.x2:
            print 'la altura de las almas no coincide, CUIDADO'
        if b != patin_abajo.x2:
            print 'los anchos de los patines no corresponden, CUIDADO'
        self.d = d
        self.b = b
        self.Ix = sum([alma_izq.Iy, alma_der.Iy, patin_arriba.area_efectiva*(0.5*(d+patin_arriba.espesor))**2 + patin_arriba.Ix, patin_abajo.area_efectiva*(0.5*(d+patin_abajo.espesor))**2 + patin_abajo.Ix])
        self.Iy = sum([patin_arriba.Iy, patin_abajo.Iy, alma_izq.area_efectiva*(0.5*(b+alma_izq.espesor))**2+alma_izq.Ix,alma_der.area_efectiva*(0.5*(b+alma_der.espesor))**2+alma_der.Ix ])
        self.rx = (self.Ix/self.area)**.5
        self.ry = (self.Iy/self.area)**.5
        self.J = self.Ix + self.Iy
        self.Cw = 1e-12
        self.Sx = self.Ix/(d+patin_arriba.espesor)
        self.Sy = self.Iy/(b+alma_der.espesor)
        self.Zx = patin_abajo.area_efectiva*(d+patin_abajo.espesor) + alma_der.area_efectiva*d/2
        self.Zy = alma_der.area_efectiva*(b+alma_der.espesor) + patin_abajo.area*b/2
        self.My = self.fy*array([self.Sx, self.Sy])
        self.Mp = self.fy*array([self.Zx, self.Zy])
        self.properties = {'rx, ry':array([self.rx, self.ry]), 'Ix, Iy, J, Cw':array([self.Ix, self.Iy, self.J, self.Cw]), 'Sx, Sy (Myx, Myy)':(array([self.Sx, self.Sy,]), self.My), 'Zx, Zy (Mpx Mpy)':(array([self.Zx, self.Zy]), self.Mp)}
        self.top_flange = patin_arriba
        self.bottom_flange = patin_abajo
        self.web = alma_izq


    def Buckling(self, effective_length, symmetries=('x', 'y'), n=1.4, FR=0.9, kx=1.0, ky=1.0, kz=1.0):
        """TODO more explicit forms of failure in dictlike form
        """
        try:
            factor_compresion = 4*self.web.esbeltez**.5 if 4*self.web.esbeltez**.5 <=0.76 else 16*self.web.esbeltez
            if factor_compresion < 0.39:
                factor_compresion = 0.39
            print 'kc = {}'.format(factor_compresion)

            if 0.5/self.top_flange.esbeltez <= 0.64*factor_compresion**.5*self.top_flange.factor_esbeltez:
                self.tipo_patin_superior_compresion = 'tipo 1,2 o 3'
            else:
                self.tipo_patin_superior_compresion = 'tipo 4'
            if 0.5/self.bottom_flange.esbeltez <= 0.64*factor_compresion**.5*self.bottom_flange.factor_esbeltez:
                self.tipo_patin_inferior_compresion = 'tipo 1,2 o 3'
            else:
                self.tipo_patin_inferior_compresion = 'tipo 4'

            if 1./self.web.esbeltez <= 1.49*self.web.factor_esbeltez:
                self.tipo_alma_compresion = 'tipo 1,2 o 3'
            else:
                self.tipo_alma_compresion = 'tipo 4'

            self.relaciones_esbeltez_compresion = {'patin_superior':self.tipo_patin_superior_compresion, 'patin_inferior':self.tipo_patin_inferior_compresion, 'esbeltez_alma':self.tipo_alma_compresion}
        except AttributeError as error:
            print error

        try:
            for placa in self.placas:
                if 1.0/placa.esbeltez < 1.40*placa.factor_esbeltez:
                    self.tipo_seccion = 'no-esbelta'
                else:
                    self.tipo_seccion = 'tipo 4'
        except AttributeError as error:
            print error

        Feft = 1e9
        Fex = Pi**2 * self.Young / (kx * effective_length / self.rx)**2
        Fey = Pi**2 * self.Young / (ky * effective_length / self.ry)**2
        if 'x' in symmetries and 'y' in symmetries:
            Fez = 1.0 / (self.Ix + self.Iy) * (Pi**2 * self.Young * self.Cw / (kz * effective_length)**2 + self.G * self.J)
        elif 'x' in symmetries:
            Fez = (self.G * self.J + Pi**2 * self.E * self.Cw /
                (kz * effective_length)**2) / (self.area * self.r0**2)
            Feft = (Fex + Fez) / (2 * self.H) * \
                (1 - (1 - (4 * Fex * Fez * self.H) / (Fex + Fez)**2)**.5)
        elif 'y' in symmetries:
            Fez = (self.G * self.J + Pi**2 * self.E * self.Cw /
                (kz * effective_length)**2) / (self.area * self.r0**2)
            Feft = (Fey + Fey) / (2 * self.H) * \
                (1 - (1 - (4 * Fey * Fez * self.H) / (Fey + Fez)**2)**.5)
        else:
            print 'non-symmetrical profiles are not yet implemented (root solving)'
        critical_stresses = {'sex': Fex, 'sey': Fey, 'sez': Fez, 'seft': Feft}
        critical_loads = {'pex': Fex * self.area, 'pey': Fey * self.area, 'pez': Fez * self.area}
        failure_stress, stress_value =min(critical_stresses, key=critical_stresses.get),min(critical_stresses.values())
        equivalent_lambda = (self.fy / stress_value)**.5
        reduccion_esbeltez = 1.0 / (1 + equivalent_lambda**(2 * n))**(1.0 / n)
        """ return areas efectivas
        """
        try:
            self.web.Resistencia(carga_axial_critica = reduccion_esbeltez*self.fy, atiesado=True)
            self.top_flange.Resistencia(carga_axial_critica = reduccion_esbeltez*self.fy, atiesado=False)
            self.bottom_flange.Resistencia(carga_axial_critica = reduccion_esbeltez*self.fy, atiesado=False)
            self.area_efectiva = self.web.area + self.top_flange.area_efectiva + self.bottom_flange.area_efectiva
        except AttributeError as error:
            print error
        try:
            if self.tipo_seccion == 'tipo 4':
                self.area_efectiva = 0
                for placa in self.placas:
                    placa.Resistencia(carga_axial_critica=reduccion_esbeltez*self.fy, atiesado=True)
                    self.area_efectiva += placa.area_efectiva
        except AttributeError as error:
            print error

        print 'area efectiva = {}'.format(self.area_efectiva)
        self.Rc = FR * reduccion_esbeltez * self.fy * self.area_efectiva
        self.critical_stresses = critical_stresses
        self.critical_loads = critical_loads

        self.Py = self.area * self.fy
        self.buckling_results = {'stresses':self.critical_stresses, 'loads':self.critical_loads, 'resistance':('plano_falla (carga ultima Rc):', failure_stress, self.Rc/1e3 ,'-Ton'), 'misc':{'Py':self.Py, 'lambda':equivalent_lambda, 'xi':reduccion_esbeltez, 'peso kg':self.area*0.00785*effective_length}}

    def Flexure(self, effective_length=100.0, C_coefficient=1.0, tipo_seccion='tipo 1', moments=[], FR = 0.9):
        """ TODO: moment distribution and autocompute C factor
        moments = [Mmax, MA, MB, MC]
        simply supported beams C=0.6.
        """
        self.tipo_patin_inferior_flexion = tipo_seccion
        self.tipo_alma_flexion = 'tipo 1'
        try:
            if 0.5/self.top_flange.esbeltez <= 0.3*self.top_flange.factor_esbeltez:
                self.tipo_patin_superior_flexion = 'tipo 1'
            elif 0.5/self.top_flange.esbeltez <= 0.38*self.top_flange.factor_esbeltez:
                self.tipo_patin_superior_flexion = 'tipo 2'
            elif 0.5/self.top_flange.esbeltez <= 1.00*self.top_flange.factor_esbeltez:
                self.tipo_patin_superior_flexion = 'tipo 3'
            else:
                self.tipo_patin_superior_flexion = 'tipo 4'
            if 0.5/self.bottom_flange.esbeltez <= 0.3*self.bottom_flange.factor_esbeltez:
                self.tipo_patin_inferior_flexion = 'tipo 1'
            elif 0.5/self.bottom_flange.esbeltez <= 0.38*self.bottom_flange.factor_esbeltez:
                self.tipo_patin_inferior_flexion = 'tipo 2'
            elif 0.5/self.bottom_flange.esbeltez <= 1.00*self.bottom_flange.factor_esbeltez:
                self.tipo_patin_inferior_flexion = 'tipo 3'
            else:
                self.tipo_patin_inferior_flexion = 'tipo 4'
            if 1./self.web.esbeltez <= 2.45*self.web.factor_esbeltez:
                self.tipo_alma_flexion = 'tipo 1'
            elif 1./self.web.esbeltez <= 3.76*self.web.factor_esbeltez:
                self.tipo_alma_flexion = 'tipo 2'
            elif 1./self.web.esbeltez <= 5.70*self.web.factor_esbeltez:
                self.tipo_alma_flexion = 'tipo 3'
            else:
                self.tipo_alma_flexion = 'tipo 4'

            self.relaciones_esbeltez_flexion = {'patin_superior':self.tipo_patin_superior_flexion, 'patin_inferior':self.tipo_patin_inferior_flexion, 'esbeltez_alma':self.tipo_alma_flexion}
        except AttributeError as e:
            print e

        self.inertias = array([self.Ix, self.Iy])

        self.Me = Pi/C_coefficient/effective_length * (self.Young*self.G*self.J * self.inertias + self.inertias*self.Cw*(Pi*self.Young/effective_length)**2)**0.5

        try:
            if self.tipo_alma_flexion or self.tipo_patin_inferior_flexion or self.tipo_patin_superior_flexion == 'tipo 4' or 'tipo 3':
                self.MR = 1.15*FR*self.My[0]*(1 - 0.28*self.My[0]/self.Me[1]) if self.Me[1]>2.0/3*self.My[0] else self.Me[1]
                if self.tipo_alma_flexion == 'tipo 4':
                    self.Me[1] = self.Me[1] * (1 - self.areas[1]/self.areas[2]/(1200+300*self.areas[1]/self.areas[2])*(self.peralte_alma/self.ta-5.6*(self.Young/self.fy)**.5)) if (1 - self.areas[1]/self.areas[2]/(1200+300*self.areas[1]/self.areas[2])*(self.peralte_alma/self.ta-5.6*(self.Young/self.fy)**.5)) <= 1.0 else self.Me[1]
            else:
                self.MR = 1.15*FR*self.Mp[0]*(1 - 0.28*self.Mp[0]/self.Me[1]) if self.Me[1]>2.0/3*self.Mp[0] else self.Me[1]

            self.flexure_results = {'Mp':self.Mp*1e-5, 'buckling moments':self.Me*1e-5, 'MR':self.MR*1e-5, 'My':self.My*1e-5, 'reduccion_esbeltez': self.tipo_alma_flexion == 'tipo 4', 'ancho/grueso':self.relaciones_esbeltez_flexion}
        except AttributeError as e:
            print e

def flexion_biaxial(Muox, MRx, Muoy=0, MRy=1.0):
    cumple = False
    try:
        print '{:.0f}/{:.0f} + {:.0f}/{:.0f} = {:.2f} + {:.2f} = {:.2f}'.format(Muox*1e-5, MRx*1e-5, Muoy*1e-5, MRy*1e-5, Muox/MRx, Muoy/MRy, Muox/MRx+ Muoy/MRy)
        cumple = Muox/MRx+ Muoy/MRy < 1.0
    except ZeroDivisionError as e:
        print e
    return 'Cumple?', cumple

def flexocompresion_biaxial(Pu, Pr, Muox, MRx, Muoy=0, MRy=1.0):
    cumple = False
    try:
        print '{:.0f}/{:.0f} + {:.0f}/{:.0f} + {:.0f}/{:.0f} = {:.2f} + {:.2f} + {:.2f} = {:.2f}'.format(Pu*1e-3, Pr*1e-3, 0.85*Muox*1e-5, MRx*1e-5, 0.6*Muoy*1e-5, MRy*1e-5, Pu/Pr, 0.85*Muox/MRx, 0.6*Muoy/MRy, Pu/Pr + 0.85*Muox/MRx + Muoy/MRy)
        print '{:.0f}/{:.0f} + {:.0f}/{:.0f} =  {:.2f} + {:.2f} = {:.2f}'.format(Muox*1e-5, MRx*1e-5, Muoy*1e-5, MRy*1e-5, Muox/MRx, Muoy/MRy, Muox/MRx + Muoy/MRy)
        cumple1 = 1.*Pu/Pr + 0.85*Muox/MRx+ 0.6*Muoy/MRy < 1.0
        cumple2 = Muox/MRx+Muoy/MRy < 1.0
    except ZeroDivisionError as e:
        print e
    print 'Cumple con1?', cumple1
    print 'Cumple sin Pu', cumple2
    return



def Soldadura(altura_filete, numero_filetes, cortantes, momentos_area, inercias, resistencia_electrodo=4900.0, verbose=True):
    """ how in the world was this computed
    resistencia is the strength per unit length (cm)
    """
    resistencia = 0.707*0.75*resistencia_electrodo*altura_filete
    flujo = array([cortantes])*array([momentos_area])/array([inercias])
    if verbose == True:
        print 'Resistencia kg/cm = {:.0f} vs.'.format(numero_filetes*resistencia), flujo
        print '{:.0f}/{:.0f} = {:.2f}, Cumple? {}'.format(sum(flujo), numero_filetes*resistencia, sum(flujo)/numero_filetes/resistencia, sum(flujo)<= numero_filetes*resistencia)
        print '-------'
    comparacion = numero_filetes*resistencia > sum(flujo)
    return sum(flujo)

# Soldadura(0.6, 2, [140e3, 38e3], [166*50,138*38], [963e3,706e3])

class Trabe_carril(object):
    def __init__(self, claro, primaria, secundaria, separacion, espesor_placa_union = 1.00, porcentaje_anclaje = 3.0, FR=0.9):
        """ el anclaje de la placa a la trabe primaria
        es con respecto al epesor, son minimo 3 cm, pero se recomienda aumentar
        """
        self.claro = claro
        self.separacion = separacion
        self.primaria = primaria
        self.secundaria = secundaria
        self.espesor_placa_union = espesor_placa_union
        anclaje = espesor_placa_union*porcentaje_anclaje
        self.ancho_efectivo_placa_union = espesor_placa_union*1.0*(primaria.web.Young/primaria.web.fy)**.5
        primaria.altura_efectiva, secundaria.altura_efectiva = 1.03*(primaria.web.Young/primaria.web.fy)**.5*array([primaria.ta, secundaria.ta])
        print '1aef {} 2 aef{}'.format(primaria.altura_efectiva, secundaria.altura_efectiva)
        longitud_placa = separacion - primaria.bf1/2 + anclaje
        area_placa = espesor_placa_union*longitud_placa
        self.area_efectiva_placa = espesor_placa_union*(anclaje + self.ancho_efectivo_placa_union)
        Ix_efectivo_placa = espesor_placa_union**3*(anclaje + self.ancho_efectivo_placa_union)/12
        Iy_placa = espesor_placa_union*longitud_placa**3/12
        primaria.centroid_t, primaria.area_efectiva, primaria.Ixt, primaria.Iyt = Tinertia(primaria.bf1, primaria.tf1, primaria.altura_efectiva, primaria.ta)
        secundaria.centroid_t, secundaria.area_efectiva, secundaria.Ixt, secundaria.Iyt = Tinertia(secundaria.bf1, secundaria.tf1, secundaria.altura_efectiva, secundaria.ta)
        self.Cx, self.Ix, self.Sx_inf = Ninertia([primaria.centroid, primaria.area, primaria.Ix], [primaria.altura + espesor_placa_union/2, self.area_efectiva_placa , Ix_efectivo_placa])

        self.Sx_sup = self.Ix/(primaria.altura + espesor_placa_union - self.Cx)
        self.Cy, self.Iy, self.Sy_izq = Ninertia([primaria.bf1/2, primaria.area_efectiva, primaria.Iyt], [primaria.bf1/2 + (separacion - longitud_placa/2), area_placa, Iy_placa], [primaria.bf1/2 + separacion, secundaria.area_efectiva, secundaria.Iyt])
        print 'centroide y2', primaria.bf1/2+separacion+secundaria.bf1/2 - self.Cy
        print 'area efectiva', secundaria.area_efectiva + espesor_placa_union*secundaria.bf2/2
        print 'cy =',primaria.bf1/2+separacion+secundaria.bf1/2 - self.Cy - secundaria.bf1/2
        self.Sy_der = self.Iy/(primaria.bf1/2+separacion+secundaria.bf1/2 - self.Cy)
        self.Mxsup, self.Mxinf, self.Myizq, self.Myder = FR*primaria.fy*array([self.Sx_sup, self.Sx_inf, self.Sy_izq, self.Sy_der])
        self.Mypatin = primaria.bottom_flange.Sy*FR*primaria.bottom_flange.fy
        # print 'Peso--', self.peso*1e-3
        print "Syizq Syder Sxinf Sxsup Sy patin (cm3)"
        print self.Sy_izq, self.Sy_der, primaria.bottom_flange.Sy, self.Sx_inf, self.Sx_sup
        print "Myizq, Myder, Mxinf, Mxsup, Mxpatin (tm)"
        print 1e-5*array([self.Myizq, self.Myder, self.Mypatin, self.Mxinf, self.Mxsup])
        print "Ix/Ix_trabe, Iy/Iy_trabe"
        print primaria.Ix/self.Ix, primaria.Iy/self.Iy

    def Fatiga(self, momento=1.0, cortante=1.0, impacto = 54e3, frenado = 13e3, peso_carro=15e3, peso_puente_grua=125e3, carga_gancho=30e3, ruedas=2, carga_rueda=50e3, ruedas_traccion = 1, altura_riel=8.9, filete = 0.6, ciclos_equivalentes_horizontales = 1e6, ciclos_laterales = 0.5e6, FC = 1.4):
        """ momento y cortante incluye CV CM sin factorizar
        """
        self.carga_impacto = array([0.25*momento, 0.25*cortante])
        self.carga_lateral = max(array([carga_gancho, 0.2*(carga_gancho + peso_carro), 0.1*(carga_gancho + peso_puente_grua)]))/carga_rueda/2/ruedas * array([momento, cortante])
        self.carga_frenado = 0.2*ruedas_traccion*carga_rueda*array([(altura_riel+self.primaria.altura)/2, (altura_riel+self.primaria.altura)/self.claro])
        # el momento maximo por frenado esta mal por que supone que ocurre al centro, lo cual no necesariamente
        self.excentricidad_cargas = 1 + altura_riel/self.primaria.altura
        self.carga_lateral = self.excentricidad_cargas * self.carga_lateral
        self.Muox = FC*(self.carga_frenado[0] + self.carga_impacto[0] + momento)
        self.Muoy = FC*self.carga_lateral[0]

        amplitud_equivalente = (ciclos_laterales/ciclos_equivalentes_horizontales)**(1./3)
        print 'excentricidad de cargas', self.excentricidad_cargas
        print 'relacion de carga (amplitud equivalente)', amplitud_equivalente
        print '--------'
        print 'REVISION RESISTENCIAS'
        print 'esquina superior lado riel sin impacto'
        print flexion_biaxial(FC*(self.carga_frenado[0]+momento), self.Mxsup, self.Muoy, self.Myizq
        )
        print 'patin inferior sin carga lateral'
        print flexion_biaxial(self.Muox, self.Mxinf, 0.0, self.Myizq)
        print 'patin inferior, cargas vivas sin impacto'
        print flexion_biaxial(FC*(self.carga_frenado[0]+momento), self.Mxinf, FC*self.carga_frenado[0], self.Mypatin)
        print 'flexion esquina superior lado de la trabe secundaria'
        print flexion_biaxial(0.0, self.Mxsup, self.Muoy, self.Myder)
        print '--------'
        print 'CORTANTE'
        print '{:.0f}/{:.0f} = {:.2f}, Cumple? {}'.format(FC*cortante*1e-3, 0.9*self.primaria.VR[0]*1e-3, FC*cortante/0.9/self.primaria.VR[0], FC*cortante/0.9/self.primaria.VR[0] < 1.0)
        print '--------'
        print 'FLUJO PLASTICO ALMA'
        print '{:.0f} < {:.0f}, Cumple? {}'.format(1.5*carga_rueda, self.primaria.flujo_plastico_alma, self.primaria.flujo_plastico_alma>1.5*carga_rueda)
        print '--------'
        print 'ABOLLAMIENTO ALMA'
        print '{:.0f} < {:.0f}, Cumple? {}'.format(1.5*carga_rueda,self.primaria.abollamiento_alma*0.75,  self.primaria.abollamiento_alma*0.75>1.5*carga_rueda)
        print '--------'
        print 'ATIESADOR ACCION COLUMNA'
        print '{:.0f} < {:.0f}, Cumple? {}'.format(FC*cortante, self.primaria.col_atiesador.Rc, self.primaria.col_atiesador.Rc>FC*cortante)
        print '--------'
        print 'ABOLLAMIENTO ATIESADOR'
        print '{:.0f} < {:.0f}, Cumple? {}'.format(FC*cortante/2, self.primaria.col_atiesador.aplastamiento, self.primaria.col_atiesador.aplastamiento>FC*cortante/2)
        print '--------'
        cortante_bidireccional = [FC*cortante, FC*self.carga_lateral[1]]
        print 'RESISTENCIA SOLDADURAS propuestas (kg/cm2)'
        print 'a) patin inferior con alma'
        q_patin_inferior = Soldadura(filete, 2, FC*cortante, self.primaria.bf2*self.primaria.tf2*(self.Cx-self.primaria.tf2/2), self.Ix)
        print 'b) patin superior con alma'
        q_patin_superior = Soldadura(filete, 2, cortante_bidireccional, [(self.primaria.bf1*self.primaria.tf1+self.area_efectiva_placa)*(self.primaria.altura - self.Cx),],[self.Ix, self.Iy])
        print 'c-d) placa respaldo con patin primaria'
        q_placa_primaria = Soldadura(filete, 2, cortante_bidireccional, [self.area_efectiva_placa*(self.primaria.altura - self.Cx + self.espesor_placa_union), self.primaria.area_efectiva*(self.Cy - self.primaria.bf1/2)], [self.Ix, self.Iy])
        print 'e) primera soldadura placa respaldon con patin secundaria'
        q_placa_secundaria = Soldadura(filete, 1, FC*self.carga_lateral[1], (self.secundaria.area_efectiva + self.espesor_placa_union*self.secundaria.bf2/2)*(self.primaria.bf1/2+self.separacion+self.secundaria.bf1/2 - self.Cy - self.secundaria.bf1/2) , self.Iy)
        print 'f) segunda soldadura placa con secundaria'
        q_segunda_placa_secundaria = Soldadura(filete, 1, FC*self.carga_lateral[1], (self.primaria.bf1/2+self.separacion+self.secundaria.bf1/2 - self.Cy - self.secundaria.bf1/2)*(self.secundaria.altura_efectiva*self.secundaria.ta + self.secundaria.bf1*self.secundaria.tf1/2), self.Iy)
        print 'soldadura atiesador'
        q_atiesador = Soldadura(filete, 2, FC*cortante/2, 1, self.primaria.web.x2)
        # carga impacto se obtiene ya de lineas de influencia
        print '--------'
        print 'CARGAS DE FATIGA'
        print 'carga lateral superior'
        momento_y = 0.5*amplitud_equivalente*FC*self.carga_lateral[0]
        print momento_y
        print 'carga lateral inferior'
        # print 0.5*amplitud_equivalente
        print 'cortante lateral'
        cortante_y = 0.5*amplitud_equivalente*FC*self.carga_lateral[1]
        print cortante_y
        print 'soldadura atiesador'
        cortante_atiesador = FC*(cortante/2+self.carga_frenado[1]+self.carga_impacto[1])
        print cortante_atiesador
        print '-------'
        print 'Rangos de esfuerzo metal base (kg/cm2)'
        print 'metal base, ultima fibra inferior'
        print momento/self.Sx_inf
        print 'fibra superior patin inferior'
        print momento/self.Ix*(self.Cx - self.primaria.tf2)
        print 'patin superior (bidireccional sin tension)'
        print '- {:.1f} +- {:.1f} = {:.1f}'.format(momento/self.Ix*(self.primaria.altura + self.espesor_placa_union - self.Cx), momento_y/self.Iy*(self.Cy-self.primaria.bf1/2), momento/self.Ix*(self.primaria.altura + self.espesor_placa_union - self.Cx) + momento_y/self.Iy*(self.Cy-self.primaria.bf1/2))
        print 'soldadura placa con primaria (sin tension)'
        print '- {:.1f} +- {:.1f} = {:.1f}'.format(momento/self.Ix*(self.primaria.altura + self.espesor_placa_union - self.Cx), momento_y/self.Iy*(self.Cy-self.primaria.bf1), momento/self.Ix*(self.primaria.altura + self.espesor_placa_union - self.Cx) +momento_y/self.Iy*(self.Cy-self.primaria.bf1))
        print 'primera soldadura placa con secundaria (reversible)'
        print '{:.1f}'.format(momento_y/self.Iy*self.primaria.bf1/2+self.separacion+self.secundaria.bf1 - self.Cy - self.secundaria.bf1)
        print 'segunda soldadura placa con secundaria (reversible)'
        print '{:.1f}'.format(momento_y/self.Iy*self.primaria.bf1/2+self.separacion+self.secundaria.bf1 - self.Cy - self.secundaria.bf1/2)
        print '--------'
        print 'CARGAS FATIGA SOLDADURA'
        # se hace otra vez flujos de cortante pero ahora con cortante_y en vez de carga lateral y con cortante sin FC sin impacto sin carga lateral?
        print 'a) patin inferior con alma'
        print Soldadura(filete, 2, cortante, self.primaria.bf2*self.primaria.tf2*(self.Cx-self.primaria.tf2/2), self.Ix, verbose=False)
        print 'b) patin superior con alma'
        print Soldadura(filete, 2, [cortante, cortante_y], [(self.primaria.bf1*self.primaria.tf1+self.area_efectiva_placa)*(self.primaria.altura - self.Cx),],[self.Ix, self.Iy], verbose=False)
        print 'c-d) soldadura placa con primaria'
        print Soldadura(filete, 2, [cortante, cortante_y], [self.area_efectiva_placa*(self.primaria.altura - self.Cx + self.espesor_placa_union), self.primaria.area_efectiva*(self.Cy - self.primaria.bf1/2)], [self.Ix, self.Iy], verbose = False)
        print 'e) primera soldadura placa con secundaria'
        print Soldadura(filete, 1, cortante_y, (self.secundaria.area_efectiva + self.espesor_placa_union*self.secundaria.bf2/2)*(self.primaria.bf1/2+self.separacion+self.secundaria.bf1/2 - self.Cy - self.secundaria.bf1/2) , self.Iy, verbose=False)
        print 'f) segunda soldadura placa con secundaria'
        print Soldadura(filete, 1, cortante_y, (self.primaria.bf1/2+self.separacion+self.secundaria.bf1/2 - self.Cy - self.secundaria.bf1/2)*(self.secundaria.altura_efectiva*self.secundaria.ta + self.secundaria.bf1*self.secundaria.tf1/2), self.Iy, verbose=False)
        print  'soldadura atiesador con alma'
        # el 140 de donde sale?
        print cortante_atiesador/2/self.primaria.altura/filete/0.707
        # carga impacto se o

# TRABE CARRIL
# -----------
fy = 3520.0
atie = Placa(fy, 1.588, 22.5*3/4)
atie.Resistencia()
patines = Placa(fy, 1.588, 55.0)
alma = Placa(fy, 1.50, 250.0-2*1.6)
pp.pprint(atie.resistencias)
tarmada = Profile()
tarmada.Plate_girder(alma, patines, patines, atie, longitud_apoyo_o_placa=29.2)
tarmada.Flexure(1800.0)
# pp.pprint(tarmada.properties)
# pp.pprint(tarmada.plate_girder_results)
# pp.pprint(tarmada.flexure_results)
# pp.pprint(tarmada.properties)
# pp.pprint(tarmada.flexure_results)
print(0.9*tarmada.Sx*3520e-5)
spatines = Placa(fy, 1.14, 16.5)
salma = Placa(fy, 0.89, 52.5-2*1.14)
tsecundaria = Profile()
tsecundaria.Plate_girder(salma, spatines, spatines, atie)
tsecundaria.Flexure(1800.0)
pp.pprint(tsecundaria.flexure_results)
pp.pprint(tsecundaria.properties)

print 'cp = {}'.format(tsecundaria.area)
tc = Trabe_carril(1800, tarmada, tsecundaria, 150.0, 1.0)
tc.Ix
tc.Iy
tc.Fatiga(436*1e5, 112e3, carga_gancho=40e3, altura_riel=8.9, filete=0.8, FC=1.0)
tc.Ix

fy = 3520.0
atie = Placa(fy, 3, 30.0)
atie.Resistencia()
tp = 1.73
patines = Placa(fy, tp, 28)
alma = Placa(fy, 1.08, 46.3-2*tp)
# col = Profile()
# col.Plate_girder(alma, patines, patines, atie)
# col.Buckling(450.0)


trabe_ppal = Profile()
trabe_ppal.Plate_girder(alma, patines, patines, atie)
trabe_ppal.Ix, trabe_ppal.Zx
trabe_ppal.VR
2*1.54*trabe_ppal.Mp/trabe_ppal.altura
# print col.properties
# print col.buckling_results
# print col.relaciones_esbeltez_compresion
# col.Flexure(450.,0.89)
# print col.flexure_results

# flexocompresion_biaxial(69e3, 1520e3, 12e5/0.85, 416e5, 33e5/0.6, 250e5)
# flexocompresion_biaxial(69e3, 1520e3, 42e5/0.85, 416e5, 95e5/0.6, 250e5)

# pats = Placa(fy, 1.1, 14)
# alma_trabe = Placa(fy, 0.7, 40)
# trabe = Profile()
# trabe.Plate_girder(alma_trabe, pats, pats, atie)
# print trabe.plate_girder_results
# trabe.Flexure(250.0, 1.0)
# print trabe.flexure_results
placas = Placa(fy, 2.54, 74.5)
cajon1 = Profile()
effective_len = 1.18
cajon1.HSS(placas, placas, placas, placas)
cajon1.Buckling(100)
cajon1.Flexure(100, 0.6)
pp.pprint(cajon1.properties)

pp.pprint(cajon1.buckling_results)
pp.pprint(cajon1.flexure_results)

FR = 0.9
alpha = 0.1
Pu = 20e3
Mu = alpha*400*Pu

PE = cajon1.buckling_results['loads']['pex']

B1 = 1./(1-Pu/FR/PE)
B2 = 1./(1-Pu*3*1.2*4/PE/2)
print B1, B2
print 'Muox', B1*B2*Mu/1e5

flexocompresion_biaxial(Pu, cajon1.Rc, B1*B2*Mu, cajon1.MR)
# TODO ayudas de diseno basado en peso por metro lineal, > 200 es muy pesado el perfil..etc
# """TODO: fem or variational methods to find the interaction coefficients for plate_girders? is this really important?
# """
# def Coef_restriccion(self, k = None):
#     if k is None:
#         if self.seccion == 'HSS':
#             if (self.ta/self.tp * self.x1/self.d) <= 1:
#                 self.coef_restriccion = (self.ta/self.tp)**3 * 0.38/(1 - (self.ta/self.tp * self.d/self.x1)**2)
#                 self.k_compresion = (2 + 2.0/(10*self.coef_restriccion + 3))**2
#                 self.esbeltez = self.ta/(self.d - 2*self.tp)
#             else:
#                 return 'las relaciones de esbeltes sobrepasan t*c/tc*d >1, el factor placa con interaccion no puede ser calculado'

#         elif self.seccion == 'I' or self.seccion == 'W':
#             if (9.4*(self.ta/self.tp * self.x1/2/self.d)**2) <=1:
#                 print(9.4*(self.ta/self.tp * self.x1/2/self.d)**2)
#                 self.placa_critica = 'alma'
#                 print('la placa del alma es critica')
#                 self.coef_restriccion = (self.ta/self.tp)**3 *(0.16+0.0056*(2*self.d/self.x1)**2)*1.0/(1-9.4*(self.ta/self.tp * self.x1/2/self.d)**2)
#                 # print((self.ta/self.tp)**3)
#                 self.k_compresion = (2 + 2.0/(10*self.coef_restriccion + 3))**2
#                 self.esbeltez = self.ta/(self.d - 2*self.tp)
#             else:
#                 self.placa_critica = 'patin'
#                 print('la placa del patin es critica')
#                 self.coef_restriccion = 4 * (self.tp/self.ta)**3 * self.d/self.x1 * 1.0/(1 - 0.106 * (2*self.tp/self.ta * self.d / self.x1)**2)
#                 self.k_compresion = (0.65 + 2.0/(3*self.coef_restriccion + 4))**2
#                 self.esbeltez = self.tp/(self.x1/2 - self.ta/2)
#         else:
#             print('esa seccion no se encuentra programada')
#             # return
#     else:
#         self.k_compresion = k

# placahss = Placa(3515.0, 0.63, 15.0)
# ptr = Profile()
# ptr.HSS(placahss, placahss, placahss, placahss)
# ptr.Buckling(150.0)
# print ptr.properties
# print ptr.buckling_results

ej5 = Placa(2530.0, 0.2, 100.0)
ej5.Resistencia()
pprint(ej5.resistencias)
