import numpy as np
from numpy.polynomial import Polynomial as P
from numpy import genfromtxt


def rect_geometry(moment):
    """estimates geometry for rectangular beam
    under PURE BENDING

=    Returns:
    height, width (cm)
    """
    height = moment / 1e5
    width = height / 2

    return height, width


def design_rect_simply(height, width, moment, fc=250,
                       fy=4200, cover=3, resistance_factor=0.9):
    """ RCDF2004
        Computes rebars for given moment(load_factor*moment) in a rectangular section

        Returns
        -------
        [(number rebars, #(/8"")), steel area (cm2), error (%),
        underly(minimally) or overly reinforced, balanced steel area (cm2),
        balanced moment]

        Notes:

    """
    # rebar_areas = np.genfromtxt('rebar_chart.csv', delimiter=',',
    # dtype=float)
    depth = height - cover
    q_min = 0.7 * np.sqrt(250) / (0.68 * fc)
    min_rebar = depth * width * q_min * 0.68 * fc / fy
    print('Min steel area: ' + str(min_rebar) + ' cm2')
    fc = 0.8 * fc
    beta = 1.05 - 0.8 * fc / 1400

    if beta < 0.65:
        beta = 0.65
    elif beta > 0.85:
        beta = 0.85
    fc = 0.85 * fc

    min_moment = resistance_factor * fc * width * \
        depth**2 * q_min * (1 - 0.5 * q_min)
    print('Min moment: ' + str(min_moment / 1e5) + ' t-m')
    balanced_area = fc * depth * width * 6000 * beta / (fy * (6000 + fy))
    balanced_ratio = balanced_area / (depth * width)
    q = balanced_ratio * fy / fc
    bal_moment = width * depth**2 * fc * resistance_factor * q * (1 - 0.5 * q)
    q_equation = P([moment, -resistance_factor * width * depth**2 * fc, 0.5])
    print('Balanced area: ' + str(balanced_area) + ' cm2')
    print('Balanced moment: ' + str(bal_moment / 1e5) + ' t-m')
    print('Max area: ' + str(balanced_area * 0.9) + ' cm2 (gravity)')
    print('Max area: ' + str(balanced_area * 0.75) + ' cm2 (seism)')
    print('Max moment: ' + str(width * depth**2 * fc * resistance_factor *
                               0.9 * q * (1 - 0.5 * 0.9 * q) / 1e5) + ' t-m (gravity)')
    print('Max moment: ' + str(width * depth**2 * fc * resistance_factor *
                               0.75 * q * (1 - 0.5 * 0.75 * q) / 1e5) + ' t-m (seism)')
    1 / fy * width * depth * fc * (1 - np.sqrt(1 - 2 * moment / (resistance_factor * width * depth**2 * fc)))
    steel_area = q * width * depth * fc / fy
    print('Steel area= ' + str(steel_area))
    """ np.nonzero(steel_area-rebar_chart) get indices
    np.nanmin, np.allclose investigate and investigate classes.
    """

    if steel_area < balanced_area:
        print('Underly reinforced, steel area : ' + str(steel_area) + ' cm2')
    else:
        print('Overly reinforced, steel area : ' + str(steel_area) + ' cm2')

    return steel_area

def design_rect_doubly(height, width, moment, fc = 250,
                       fy = 4200, cover = 3, resistance_factor = 0.9):
    """ RCDF2004
        Computes rebars for given moment(load_factor*moment) in a doubly
        reinforced rectangular section

        Returns:
        [(number rebars, #(/8"")), steel area (cm2), error (%),
        underly(minimally) or overly reinforced, balanced steel area (cm2), balanced moment(t-m)
        max moment (gravity), max moment (seism)]

        Notes:
    """
        # rebar_areas = np.genfromtxt('rebar_chart.csv', delimiter=',',
        # dtype=float)
    depth=height - cover
    q_min=0.7 * np.sqrt(250) / (0.68 * fc)
    min_rebar=depth * width * q_min * 0.68 * fc / fy
    print('Min steel area: ' + str(min_rebar) + ' cm2')
    fc=0.8 * fc
    beta=1.05 - 0.8 * fc / 1400

    if beta < 0.65:
        beta=0.65
    elif beta > 0.85:
        beta=0.85
    fc=0.85 * fc

    min_moment=resistance_factor * fc * width * \
        depth**2 * q_min * (1 - 0.5 * q_min)
    print('Min moment: ' + str(min_moment / 1e5) + ' t-m')
    balanced_area=fc * depth * width * 6000 * beta / (fy * (6000 + fy))
    balanced_ratio=balanced_area / (depth * width)
    q=balanced_ratio * fy / fc
    bal_moment=width * depth**2 * fc * resistance_factor * q * (1 - 0.5 * q)
    q_equation=P([moment, -resistance_factor * width * depth**2 * fc, 0.5])
    print('Balanced area: ' + str(balanced_area) + ' cm2')
    print('Balanced moment: ' + str(bal_moment / 1e5) + ' t-m')
    print('Max area: ' + str(balanced_area * 0.9) + ' cm2 (gravity)')
    print('Max area: ' + str(balanced_area * 0.75) + ' cm2 (seism)')
    print('Max moment: ' + str(width * depth**2 * fc * resistance_factor *
                               0.9 * q * (1 - 0.5 * 0.9 * q) / 1e5) + ' t-m (gravity)')
    print('Max moment: ' + str(width * depth**2 * fc * resistance_factor *
                               0.75 * q * (1 - 0.5 * 0.75 * q) / 1e5) + ' t-m (seism)')
    str(1 / fy * width * depth * fc * (1 - np.sqrt(1 - 2 * moment / (resistance_factor * width * depth**2 * fc))))
    steel_area=q * width * depth * fc / fy
    print('Steel area= ' + str(steel_area))
    """ np.nonzero(steel_area-rebar_chart) get indices
    np.nanmin, np.allclose investigate and investigate classes.
    """

    if steel_area < balanced_area:
        print('Underly reinforced, steel area : ' + str(steel_area) + ' cm2')
    else:
        print('Overly reinforced, steel area : ' + str(steel_area) + ' cm2')

    return steel_area


def analysis_rect_simply_underly(height, width, steel_area, fc=250,
                               fy=4200, cover=3, resistance_factor=0.9):
    """moment (kg-cm) given steel area (cm2) (RCDF)
    """
    depth=height - cover
    ratio=steel_area / (width * depth)
    fc=0.85 * fc
    q=ratio * fy / fc
    return resistance_factor * width * depth**2 * fc * q * (1 - 0.5 * q)/1e5


def analysis_rect_simply_overly(height, width, steel_area, fc=250,
                              fy=4200, cover=3, resistance_factor=0.9):
    """computes moment (kg-cm) given steel area (cm2) (RCDF)
    """
    beta=1.05 - 0.8 * fc / 1400

    if beta < 0.65:
        beta=0.65
    elif beta > 0.85:
        beta=0.85

    depth=height - cover
    fc=0.68 * fc
    f_eq=P([-6000 * depth * width * beta * fc / steel_area, 6000, 1])
    print(f_eq.roots())
    fs=np.select([f_eq.roots() > 0], [f_eq.roots()], 0)
    a=steel_area * fs / (width * fc)
    return a * width * fc * (depth - 0.5 * a) * resistance_factor


class simple_concrete_column(object):
    def __init__(self, x1, x2, fc):
        """ x1 x2 fc is maximum compressive resistance
        """
        self.area = x1*x2
        self.A = x1*x2
        self.fpc = fc
        self.fc = 0.68*fc
        self.x1 = x1
        self.x2 = x2
        self.I1 = float(x2)*x1**3/12
        self.I2 = float(x1)*x2**3/12
        self.Pr = self.area*self.fc
        self.perimeter = 2*(x1 + x2)
        self.Ec = 14000*(fc)**.5
        self.E = 14000*(fc)**.5
        self.flexural_stiffness1 = self.Ec * self.I1
        self.flexural_stiffness2 = self.Ec * self.I2
        self.plastic_axial_load = self.area*self.fc
        self.Pu = self.area*self.fc
        self.Zx = x1*x2**2/4
        self.Zy = x2*x1**2/4

    def rebar(self, designation=3, number_x1=2, number_x2 = 2, cover=3.0, fy = 4200.0, Es = 2.036e6):
        """accomodates rebar in perimeter
        """
        self.fyr = fy
        self.Es = Es
        rebars = 2*(number_x1 + number_x2 - 2)
        diameter = 2.54/8*designation
        rebar_area = pi*diameter**2/4
        self.rebar_area = rebar_area
        self.rebars = rebars
        self.rebar_plastic_axial_load = rebar_area*rebars*fy
        self.rebar_Pu = rebar_area*rebars*fy
        self.rebar_area = rebar_area
        x1eff = self.x1 - 2*cover
        x2eff = self.x2 - 2*cover
        x1sep = x1eff/(number_x1-1)
        x2sep = x2eff/(number_x2-1)
        bar_inertia = pi*(diameter/2)**4/4
        inertia_rebars = bar_inertia *rebars
        self.x1bars2center = array([0.5*self.x1 - cover - n*x1sep for n in range((number_x1)/2)])
        self.x2bars2center = array([0.5*self.x2 - cover - n*x2sep for n in range((number_x2)/2)])
        self.x1resistance = (number_x1-2)*rebar_area*fy
        self.x2resistance = (number_x2-2)*rebar_area*fy
        self.rebar_MRx = 2*rebar_area*fy*(2*sum(self.x2bars2center) + (number_x1 - 2)*self.x2bars2center[0])
        self.rebar_MRy = 2*rebar_area*fy*(2*sum(self.x1bars2center) + (number_x2 - 2)*self.x1bars2center[0])
        self.rebar_I1 = 2*number_x1*rebar_area*x1eff*0.5 + inertia_rebars + 4*sum(self.x2bars2center[1:]**2)
        self.rebar_I2 = 2*number_x2*rebar_area*x2eff*0.5 + inertia_rebars + 4*sum(self.x1bars2center[1:]**2)
        self.rebar_flexural_stiffness1 = self.rebar_I1*Es
        self.rebar_flexural_stiffness2 = self.rebar_I2*Es



from __future__ import division
class rcBeam(object):
    """ properties at i
    requires nominal strength f'c
    """

    def __init__(self, nominal_compresive_strength_fc, unitWeight):
        self.fc = nominal_compresive_strength_fc
        self.unitWeight = unitWeight

    def rectangular(self, height, width, length, cover=3, centroid_comp=3):
        """centroid_comp is distance from edge to compression rebar centroid
        """
        self.shape = 'Rectangular'
        self.height = height
        self.width = width
        self.length = length
        self.cover = cover
        self.depth = height - cover
        self.area = height * width
        self.volume = self.area * length
        self.weight = self.volume * self.unitWeight
        self.centroid_comp = centroid_comp
        self.inertia2 = height ** 3 * width / 12
        self.inertia3 = width**3 * height / 12
        self.effective_area = self.depth * self.width


class rcSlab(object):

    def __init__(self, area, thickness, unitWeight):
        self.area = area
        self.thickness = thickness
        self.unitWeight = unitWeight
        self.volume = area*thickness
        self.weight = unitWeight* self.volume



beam1 = rcBeam(200, 2.4)

beam1.rectangular()


colobjects = []
class rc_beam(object):
    """ properties at i
    requires nominal strength f'c
    """

    def __init__(self, nominal_compresive_strength_fc):
        self.area = []
        self.inertia2 = []
        self.inertia3 = []
        self.fc = nominal_compresive_strength_fc

    def rectangular(self, height, width, cover=3, centroid_comp=3):
        """centroid_comp is distance from edge to compression rebar centroid
        """
        self.shape = 'Rectangular'
        self.height = height
        self.width = width
        self.cover = cover
        self.depth = height - cover
        self.area = height * width
        self.centroid_comp = centroid_comp
        self.inertia2 = height ** 3 * width / 12
        self.inertia3 = width**3 * height / 12
        self.effective_area = self.depth * self.width

    def T_shape(self, height, web_width, flange_width, flange_thickness,
                cover=3, centroid_comp=3):
        """height includes web + flange

        PENDING: inertia2,3
        """
        self.shape = 'T'
        self.height = height
        self.cover = cover
        self.depth = height - cover
        self.area = flange_width * flange_thickness + \
            (height - flange_thickness) * web_width
        self.web_width = web_width
        self.flange_width = flange_width
        self.flange_thickness = flange_thickness
        self.width = web_width
        self.centroid_comp = centroid_comp
        self.effective_area = self.flange_width * self.depth

    def analysis_RCDF(self, tension_area=0, compression_area=0,
                      fyt=4200, fyc=4200, FR=0.9, steel_elastic_modulus=2e6, tolerance=0.01):
        """ Returns moment2 according to NTC2004 Mexican BC

        self:
        ----
        fcc (effective strength)
        As, Asc
        p, pc
        beta

        bal_ratio
        bal_steel, bal_moment
        min_steel , min_moment,
        max_steel_gravity, max_moment_gravity,
        max_steel_seism, max_moment_seism

        state (overly, underly)
        """
        self.As = tension_area
        self.Asc = compression_area
        self.p = tension_area / self.area
        self.pc = compression_area / self.area

        if 0.8 * self.fc <= 280:
            self.beta = 0.85
        else:
            self.beta = 1.05 - 0.8 * self.fc / 1400
            if self.beta < 0.65:
                self.beta = 0.65
        # if rcdf then
        self.fcc = 0.8 * 0.85 * self.fc
        # if aci318 then
        # self.fcc = 0.8 * self.fc
        self.bal_ratio = self.fcc * 6000 \
            * self.beta / (fyc * (6000 + fyc))
        self.bal_steel = self.bal_ratio * self.area
        q = self.bal_ratio * fyc / self.fcc

        if self.p > self.bal_ratio:
            print('Overly-reinforced')
            self.state = 'overly'
        elif self.p < self.bal_ratio:
            print('Underly-reinforced')
            self.state = 'underly'

        self.min_steel = 0.7 * self.fc**.5 * self.effective_area / fyc
        self.min_moment = self.Moment(self.min_steel)
        self.bal_moment = self.width * (self.effective_area / self.width)**2 * \
            self.fcc * FR * q * (1 - 0.5 * q)
        self.max_steel_gravity = 0.9 * self.bal_steel
        self.max_moment_gravity = self.Moment(
            self.max_steel_gravity, 0, fyt, fyc, FR, steel_elastic_modulus, tolerance)
        self.max_steel_seism = 0.75 * self.bal_steel
        self.max_moment_seism = self.Moment(
            self.max_steel_seism, 0, fyt, fyc, FR, steel_elastic_modulus, tolerance)
        self.moment2 = self.Moment(
            tension_area, compression_area, fyt, fyc, FR, steel_elastic_modulus, tolerance)

        return self.moment2

    def Moment(self, tension_area, compression_area=0, fyt=4200, fyc=4200,
               FR=0.9, steel_elastic_modulus=2e6, tolerance=0.01):
        """ returns resistant moment
        applies general equilibrium and continua

        self:
        --------
        tension_strain
        compression_strain
        """
        c = self.height * 0.2
        Tension = (self.depth - c) * 0.003 / c * steel_elastic_modulus
        Compression = self.beta * c * self.width * self.fcc

        while abs(Tension - Compression) / Tension > tolerance:
            tension_strain = (self.depth - c) * 0.003 / c
            compression_strain = (c - self.centroid_comp) * 0.003 / c
            self.tension_strain = tension_strain
            self.compression_strain = compression_strain
            # Tension = (self.depth - c) * 0.003 / c * steel_elastic_modulus
            if tension_strain < fyt / steel_elastic_modulus:
                fyt = tension_strain * steel_elastic_modulus
            if compression_strain < fyc / steel_elastic_modulus:
                fyc = compression_strain * steel_elastic_modulus
                Tension = self.p * self.area * fyt
# something is strange here, perhaps could be clearer about upper bound on strenghts
            if self.shape == 'Rectangular':
                Compression = self.beta * c * self.width * \
                    self.fcc + fyc * self.pc * self.area
                moment2 = FR * (self.beta * c * self.width * self.fcc *
                                (c - self.beta * c / 2) +
                                self.Asc * fyc * (c - self.centroid_comp) +
                                self.As * fyt * (self.depth - c))
            elif self.shape == 'T':
                if c * self.beta <= self.flange_thickness:
                    Compression = self.beta * c * self.flange_width * self.fcc
                    moment2 = FR * Compression * \
                        (c - self.flange_thickness / 2)
                else:
                    comp_flange = self.flange_width * \
                        self.flange_thickness * self.fcc
                    comp_rect = (self.beta * c - self.flange_thickness) * \
                        self.web_width * self.fcc
                    Compression = comp_flange + comp_rect
                    moment2 = FR * (comp_flange * (c - self.flange_thickness / 2) +
                                    comp_rect * (c - (c * self.beta - self.flange_thickness) / 2) +
                                    self.As * fyt * (self.depth - c))
            if Tension > Compression:
                c *= (1 + (Tension - Compression) / Tension)
            else:
                c *= (1 + (Tension - Compression) / Compression)
            print('block depth:  ' + str(c))
            print('T :  ' + str(Tension))
            print('C : ' + str(Compression))

        return moment2

    def design_RCDF(self, moment):
        """ returns As, Asc for given section, properties and moment
        self:
        ----
        As (steel area in tension)
        Asc (steel area in compression)
        """
        self.As = 1
        moment2 = self.Moment(self.As)

        while moment2 < moment:
            self.As += 0.5
            moment2 = self.Moment(self.As)
            if self.As > self.bal_steel:
                self.As = self.bal_steel
                self.Asc = 0.5
                moment2 = self.Moment(self.As, self.Asc)
                if moment2 < moment:
                    self.Asc += 0.5
                    moment2 = self.Moment(self.As, self.Asc)
        return self.As, self.Asc


beam1 = rc_beam(200)
beam1.rectangular(50,30)
print(beam1.area, beam1.shape)


print(analysis_rect_simply_underly(65, 30, 15, 200, 4200, 5,))
print(design_rect_doubly(60, 30, 44000))
