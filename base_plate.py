from numpy import array
""" TODO: import plate object, import profile object, solid concrete object import forces object? [P3, M1, M2]
Coalescence of 4 objects and their interaction!
is this an Operator? is it linear? what are its Eigenvalues

"""

def moment_plate(x1, x2, anchor_distance, d, bp, axial_load = 1.0, moment = 0, fc = 300.0, fy = 4200.0):
    fpu = 0.85*fc
    critial_len = max(0.5*(x1-0.95*d), 0.5*(x2 - 0.8*bp), 0.25*(d*bp)**.5)
    excentricity = 1.0*moment/axial_load
    ecrit = x1/2 - axial_load/2/x2/fpu
    if (anchor_distance + x1/2)**2 < 2*axial_load*(excentricity + anchor_distance)/x2/fpu:
        print 'plate size insufficient'
        return
    bearing_len = x1 - 2*excentricity if excentricity <= ecrit else (anchor_distance + x1/2) - ((anchor_distance + x1/2)**2 - 2*axial_load*(excentricity + anchor_distance)/x2/fpu)**.5
    rod_tension = x2*fpu*bearing_len - axial_load if excentricity > ecrit else 0
    thickness_check = (4*rod_tension*anchor_distance/x2/0.9/fy)**.5
    print 'thickness should be', thickness_check
    return rod_tension

def  rod_check(t1,t2=0, depth=1.0, distance_to_border=1.0,  rod_number=4, fc=300.0,  fu=4900.0):
    """ TODO: rod &  conc resistance factors fn
    """
    tension = t1+t2
    # f1 =  1.0
    # f2 =1.0
    # f3 = 1.0
    area = 3.1416*(rod_number*1.0/8*2.54)**2
    rod_tension = 0.75**2*area*fu
    print 'rod resistance to tension', rod_tension
    extraction = 0.7*1.4*8*fc*area
    print 'extraction resistance', extraction
    cone_resistance = 7*(fc*depth**3)**.5
    print 'cone resistance', cone_resistance
    # concrete_cracking =
    if distance_to_border < 0.4*depth:
        border_resistance = 0.7*42*distance_to_border*(area*fc)**.5

    else:
        border_resistance = 10e8

    print 'border_resistance', border_resistance
    resistance = min(rod_tension, extraction, cone_resistance, border_resistance)

    return resistance > tension


print moment_plate(80.0, 80., 30.0, 50., 50., 125000., 60e5)
print moment_plate(80.0, 80., 30.0, 50., 50., 125000., 20e5)
print rod_check(22e3, 9e3, 50.0, 10, 6)



# def flexural_baseplate(x1, x2, thickness=1.0, anchor_distance=30.0, d=50.0, bp=50.0, axial_load = 0, M1 = 0, M2 = 0, fc = 300.0, fy = 4200.0):
#     """x1 is the smaller dimension (really? yes),
#     loads are given in kg, kg-cm
#     anchor_distance is measured to column centerline in cm
#     """
#     fpu = 0.85*fc
#     dims = array([x1, x2])
#     qmax = fpu*dims
#     cantilever_len = max([0.5*(x1-0.95*d), 0.5*(x2-0.8*bp), 0.25*(d*bp)**.5])
#     excentricities = array([M1/axial_load, M2/axial_load])
#     critical_excentricity = 0.5*(dims - axial_load/qmax)
#     # effective_compression = [dim - 2*ex if ex < ecrit else (anchor_distance + dim/2) - ((anchor_distance+dim/2)**2 - 2*axial_load*(anchor_distance + ex)/dim/fpu)**.5 for dim,ex,ecrit in dims, excentricities, critical_excentricity]
#     # bearing_len = [x1 - 2*ex if ex <= critical_excentricity else (anchor_distance + x1/2) - ((anchor_distance + x1/2 )**2 - 2*axial_load*(ex+anchor_distance)/qmax)**.5 for ex in excentricities]
#     # rod_tension = [2*axial_load*(ex+anchor_distance)/qmax*Y for ex,Y in excentricities, bearing_len]
#     # # check whether thickness is governed by bearing or tension inferface
#     # m = (x1 - 0.95*d)/2
#     # n = (x2 - 0.8*bp)/2
#     # bearing_interface = 1.5*m*(fpu/fy)**.5 if max(bearing_len) >= m else 2.11(fpu*max(bearing_len)*(m - max(bearing_len)/2)/fy)**.5, 'bearing'
#     # tension_interface = 2.11*n*(max(rod_tension)*(0.5*(x1 - d) - anchor_distance)/x2/fy)**.5, 'tension'
#     print effective_compression
#     # print bearing_interface, tension_interface
#     print 1.5*m*(fpu/fy)**.5, 2.11*n*(max(rod_tension)*(anchor_distance - d/2)/x2/fy)**.5

#     # plate_thickness = 


    # print bearing_len, rod_tension, plate_thickness
    # M_u = [axial_load/x2/Y*(cantilever_len**2/2) if Y >= cantilever_len else axial_load/x2/Y*Y*(cantilever_len-Y/2) for Y in effective_compression]
    # MR = 0.9*thickness**2*fy/4
    # print M_u, MR
    # thickness_check = (4*max(M_u)/fy/0.9)**.5
    # print thickness_check


# plate = flexural_baseplate(100, 100, 2.86 , 30, 50, 50, 125e3, 20e5, 60e5)

# def rod_check(diameter=0...)
