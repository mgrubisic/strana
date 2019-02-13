from dataclasses import dataclass, field
# from typing import NewType, Any, Generic
from numpy import array, zeros, isclose
sys.path.append('C:/Users/vyraj/Dropbox/str/scripts/')

# all loading is truly dynamic in nature, static is when there is no or only 1 time step, IDAs are a general case of nonlinear dynamic analyses
@dataclass
class Concrete:
    """regular resistance concrete"""
    resistance: float = 30.0
    stiffness: float = 14000*resistance**0.5
    unit_weight: float = 24.0
    material: str = 'Concrete01'
    strain_crack: float = - 0.002
    strain_ult: float = - 0.003

@dataclass
class Conc_Beam(Concrete):
    """bernoulli beam"""
    length: float = 1.0
    section: list = field(default_factory=[1.0, 1.0])
    cracked_inertia: float = 0.5
    steel: list field(default_factory=[0.01, 0.01])
    residual_strength: float = 0.1 # this is in reality a steel_residual_str
    def __post_init__(self):
        self.inertia = self.cracked_inertia*self.section[0]**3*self.section[1]/12
        self.area = self.section[0]*self.section[1]
        self.volume = self.length * self.area
        self.weight = self.volume*self.unit_weight
        self.mass = self.weight/9.81
    def build(self, mattag):
        section_tag = 10*mattag
        uniaxialMaterial(self.material, mattag, self.resistance, self.strain_crack, self.residual_strength**self.resistance, self.strain_ult)
        uniaxialMaterial('Steel01', section_tag, 420e3, 200e9, 0.05)
        section('RCSection2d', section_tag, mattag, mattag, section_tag, 0.95*self.section[1], self.section[0], 0.05*self.section[1], self.area*self.steel[0], self.area*self.steel[1], 0, 5, 5, 1)
        return section_tag


@dataclass
class Conc_Column(Concrete):
    """simple uniaxial column"""
    height: float = 1.0
    section: list = field(default_factory=list)
    cracked_inertia: float = 0.7
    steel: float = 1.0
    residual_strength: float = 0.1 # this is in reality a steel_residual_str
    def __post_init__(self):
        self.inertia = self.cracked_inertia*self.section[0]**3*self.section[1]/12
        self.area = self.section[0]*self.section[1]
        self.volume = self.length * self.area
        self.weight = self.volume*self.unit_weight
        self.mass = self.weight/9.81
    def build(self, mattag):
        section_tag = 10*mattag
        uniaxialMaterial(self.material, mattag, self.resistance, self.strain_crack, self.residual_strength**self.resistance, self.strain_ult)
        uniaxialMaterial('Steel01', section_tag, 420e3, 200e9, 0.05)
        section('RCSection2d', section_tag, mattag, mattag, section_tag, 0.95*self.section[1], self.section[0], 0.05*self.section[1], self.area*self.steel, self.area*self.steel, self.area*self.steel, 5, 5, 1)
        return section_tag


@dataclass
class Frame2D:
    """1 storey 1 bay equal column frame (col, beam, uniform_load)"""
    column: object
    beam: object
    uniform_load: float = 1.00
    nodes = {}
    def build(self, self_weight = False):
        if self_weight == False:
            self.mass = self.uniform_load*self.beam*length/9.81
        else:
            self.mass = self.uniform_load*self.beam*length/9.81 + self.beam.mass + self.column.mass
        self.width = self.beam.length
        self.height = self.column.height
        return

@dataclass
class Building2D:
    """identical frame2D: number_of_bays, framelist
    example 3 bays 4 storeys Building2D(3, 4, [fr1, fr2, fr3, fr4]) """
    model('basic', '-ndm', 2)
    bays: int
    frames: list
    # nodes = {}
    def build(self):
        # TODO: check at the end whether columns x's line up
        # column_check = []
        # if False in column_check:
        #     print "columns do no line up, check frame.beam.lengths"
        geomTransf('Linear', 1)
        geomTransf('PDelta', 2)
        nodetag = 0
        eletag = 0
        height = 0
        masters = []
        for fr, frame in enumerate(frames):
            beam_section_tag = frame.beam.build(mattag = id(frame.beam))
            # use unique object id as a material (internally) and returns a section tag to be used for element()
            column_section_tag = frame.column.build(mattag = id(frame.column))
            # diagonal_section_tag = frame.diag.build(id(frame.diag))

            for bay in range(bays):
                node(nodetag, bay*frame.width, height)
                if fr == 0:
                    fix(nodetag, [1,1,1])
                    nodetag += 1
                    pass
                if fr == len(frames):
                    continue
                else:
                    col = frame.column
                    element('beamWithHinges', eletag, nodetag - bays - 1, nodetag, column_section_tag, 0.1, col.stiffness, col.area, col.inertia, 2)
                    eletag += 1

                if bay == 0:
                    mass(nodetag, bays*frame.mass)
                    masters.append(nodetag)
                    continue
                else:
                    beam = frame.beam
                    element('beamWithHinges', eletag, nodetag-1, nodetag, beam_section_tag, 0.1, beam.stiffness, beam.area, beam.inertia, 1)
                    eletag += 1
                    equalDOF(masters[fr], nodetag, 1)

                if bay == bays:
                # if fr == 0 loop gets skipped
                    for b2 in range(bays):
                        timeSeries('Linear', 1)
                        pattern('Plain', fr, 1)
                        eleLoad('-ele', eletag - bays - b2, '-type', 'beamUniform', frame.beam.load)

                height += frame.height
                nodetag += 1

        return masters

    # def pushover(self, system="BandGeneral", numberer="RCM", constraints="Plain" control = "LoadControl", algorithm="Newton", analysis="static", integrator="DisplacementControl"):
        # select master nodes ()
        # select eigenvector first mode shape ()
        # push and monitor, output drifts phlogs?
    # def ida(...)
    """transient nonlinear analysis is an unscaled ida."""
    # def adaptive_pushover...


col1 = Conc_Column(10,[15,5])
beam1 = Conc_Beam(10, [15,5])
frame1 = Frame2D(col1, col1, beam1)
frame1.build()
# print(frame1.nodes)
b1 = Building2D(3, [frame1, frame1])
b1.build()


kaushik = [4.2, 3.2, 3.2, 3.2],[3,3,3]
diaz = [3*fr for fr in range(9)] [3,3,3]
