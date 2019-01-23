from numpy import deg2rad, rad2deg
class Footing(object):
    """ x1<x2
    """
    def __init__(self, x1, x2, depth = 0.):
        self.x1 = x1
        self.x2 = x2
        self.x3 = depth
        self.depth = depth
        self.area = x1 * x2

    def loading(self, axial_load, M2, M1, load_factor = 1.1):
        self.axial_load = axial_load
        self.e1 = M2 / axial_load
        self.e2 = M1 / axial_load
        self.x1_eff = self.x1 - 2 * self.e1
        self.x2_eff = self.x2 - 2 * self.e2
        self.eff_area = self.x1_eff * self.x2_eff
        self.ult_pressure = axial_load/self.eff_area*load_factor

    def bearing_capacity(self, cohesion=0, friction=0, unit_weight = 17.0, resistance_factor = 0.7):
        self.cohesion = cohesion
        self.form_factor_fc = 1 + (self.x1_eff/self.x2_eff + self.x3/self.x1_eff)/4
        self.ult_resistance = (2 + pi)*cohesion*self.form_factor_fc*resistance_factor + self.x3*unit_weight

    def tilting(self, moment, elastic_shear_modulus = 3300.0, elastoplastic_shear=2700.0, poisson=0.5):
        self.equivalent_radius = ((self.x2*self.x1**3)/(3*pi))**(1./4)
        self.tilt = 3*(1-poisson)*moment*(elastic_shear_modulus/elastoplastic_shear - 1)/(8*elastic_shear_modulus*self.equivalent_radius**3)


def momento_volteo(alturas, pesos, coeficiente_sismico, ductilidad=2):
    # niveles = len(alturas)
    alturas = array(alturas)
    pesos = array(pesos)
    peso_total = sum(pesos)
    fuerzas_inercia = coeficiente_sismico*1./ductilidad*alturas*pesos*peso_total/sum(pesos*alturas)
    momentos_volteo = fuerzas_inercia*alturas
    return sum(momentos_volteo)

# losa = Footing(20.0, 30.6, 3.0)
# peso = 51e3
# # Mx = momento_volteo([3, 6, 9, 12, 15, 18, 21], [peso/7 for i in range(7)], 0.32)
# # print(losa.area)
# losa.loading(peso, 85337.0, 0.3*85337.0)
# losa.bearing_capacity(21.4)
# print(losa.ult_pressure)
# # print(losa.form_factor_fc)
# print(losa.ult_resistance)
altura = 24+3.5
losa = Footing(24.0, 36.0, 3.5)
presion = 81.0
peso = presion*losa.area
momento = 2./3*altura*peso*0.16
print(momento)
losa.loading(peso, momento, 0.3*momento, 1.1)
losa.bearing_capacity(21.4)
print(losa.form_factor_fc)
print(losa.x1_eff, losa.x2_eff)
print(losa.ult_pressure)
print(losa.ult_resistance)


losa.tilting(momento)
print(losa.equivalent_radius, losa.tilt)

class Footing(object):

    def __init__(self, unit_weight, shape='rectangular'):
        self.unit_weight = unit_weight
        self.shape = shape

    def geometry(self, height, minor_len, mayor_len, slope=0):

        self.height = height
        self.minor_len = minor_len
        self.mayor_len = mayor_len

        self.volume = height * minor_len * mayor_len
        self.footingarea = minor_len * mayor_len

    def static_plasticity_factors(self):

        self.fc = 1 + \
            (self.minor_len / self.mayor_len + self.depth / self.minor_len) / 4
        self.fq = 1 + self.minor_len / self.mayor_len * self.fric_slope
        self.fy = 1 - 0.4 * self.minor_len / self.mayor_len


        # try self.moment2:
            # except
        if not hasattr(self, moment2):
            if not hasattr(self, moment3):
                return

        self.fcd = 1 + (self.eff_minor_len / self.eff_mayor_len +
                        self.depth / self.eff_minor_len) / 4

        self.fqd = 1 + self.eff_minor_len / \
            self.eff_mayor_len * self.fric_slope

        self.fyd = 1 - 0.4 * self.eff_minor_len / self.eff_mayor_len

    def acting_loads(self, axial, moment2, moment3):
        """ axial load is acting on plastic centroid
            moment 2 around x2

            self.
            effective lengths
            real form factors
        """
        self.axial_load = axial
        self.moment2 = moment2
        self.moment3 = moment3
        excentricity_len2 = moment3 / axial
        excentricity_len3 = moment2 / axial
        self.eff_minor_len = self.minor_len - 2 * excentricity_len2
        self.eff_mayor_len = self.mayor_len - 2 * excentricity_len3

    def foundation_builder(self, Footing, Column, Soil_layers, footing_depth, groundwater_depth):
        """ arranges the foundation's parts on z-axis
        """


class SoilStrata():

class FootingArray():

class Foundation(FootingArray, Columns, SoilStrata, foundation_depth, groundwater_level):
    """ arranges foundation's parts on z-axis
    """


class Column(object):

    def __init__(self, unit_weight):
        self.unit_weight = unit_weight

    def geometry(self, height, minor_len, mayor_len, slope=0):
        self.height = height
        self.minor_len = minor_len
        self.mayor_len = mayor_len
        self.volume = height * minor_len * mayor_len
        self.weigth = self.unit_weight * gravity * self.volume


class Soil(object):
    """ general model for homogenous isotropic
    """

    def __init__(self, thickness):
        self.thickness = thickness

    def shear_strength(self, cohesion=0, internal_friction_angle=0):
        self.cohesion = cohesion
        self.phi = internal_friction_angle

    def in_situ(self, in_situ_friction_angle, compaction=1):

        if compaction > 0.67:
            alpha = 1
        elif compaction < 0.67 and compaction > 0:
            alpha = compaction + 0.67 - 0.75*compaction**2
            if alpha > 1:
                alpha = 1
        else:
            raise Exception('compaction can not be negative')

        self.in_situ_phi = in_situ_friction_angle
        self.compaction = compaction
        self.phi = np.arctan(alpha * np.tan(in_situ_friction_angle))
        self.passive_coeff = np.tan(np.pi/4-self.phi/2)**2

    def physical_prop(self, unit_weight, dry_unit_weight=1, saturation_rate=0, void_ratio=0):
        self.unit_weight = unit_weight
        self.unit_weigth = unit_weight * gravity
        self.saturated_weigth =
        self.bulk_weigth =

        """ saturated weigth
        dry weigth
        buoyant weigth
        """
    def plasticity_factors(self):
        self.Nq = np.exp(np.pi*np.tan(self.phi))*self.passive_coeff
        self.Ny = 2 * (self.Nq + 1)* np.tan(self.phi)
        self.Nc = (self.Nq - 1)/np.tan(self.phi)


    # def mechanical_properties(self):
    #     self.pressure_bottom = self.eff_unit_weigth * self.thickness

# class Soil_strata(object):
    # """ add soils in order from top to bottom
    # """
#   def __init__(self, [list_of_Soil_strata_objects_in_order]):

inertialfor = Inertial_for(9.81)

partially_saturated_sand = Soil(0.45)
partially_saturated_sand.in_situ(deg2rad(36), 0.65)
partially_saturated_sand.physical_prop(2.62, 16/gravity, 0.7)

saturated_sand = Soil(3)
