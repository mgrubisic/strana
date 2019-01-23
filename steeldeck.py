import stl_beamcolumn

def steeldeck(profile, effective_width=100.0, length=1000.0,  fc = 300.0, slab_thickness = 10.0, separation_between_slab_beam=0, tolerance=500.0):
    """returns resisting moment for a composite steel beam-concrete slab configuration
    assumes Ec = 14000*sqrt(fc)
    """
    compression = 1e5
    tension = 0.
    c = slab_thickness
    while abs(compression - tension) > tolerance:
        print 'c = {}'.format(c)
        a = 0.85*c
        compression = 0.68*fc*a*effective_width
        centroids = array([slab_thickness + separation_between_slab_beam + profile.tf/2,
                           slab_thickness + separation_between_slab_beam + profile.depth/2,
                           slab_thickness + separation_between_slab_beam + profile.depth - profile.tf/2 ])
        strains = 0.003*(centroids/c - 1)
        areas = array([profile.area_flange, profile.area_web, profile.area_flange])
        tensions = array([0,0,0])
        for index, strain in enumerate(strains):
            force = profile.E*strain*areas[index]
            yield_force = profile.fy*areas[index]
            if force < yield_force:
                tensions[index] = force
            else:
                tensions[index] = yield_force

        tension = sum(tensions)
        if compression > tension:
            c = c - 0.01
        else:
            c = c + 0.01

        centroid_tension = sum(tensions*areas)/sum(tensions)
        lever_arm = centroid_tension - 0.85*c/2
        Mn = compression*lever_arm
        # static deflections
        Ec = 14000*(fc)**.5
        n = Ec/profile.E
        btr = effective_width*n
        area_conc = slab_thickness*btr
        new_centroid = (profile.area*profile.depth/2 + area_conc*(slab_thickness/2+separation_between_slab_beam + profile.depth))/(area_conc+profile.area)
        inertia_concrete = btr * slab_thickness**3/12
        Ic = inertia_concrete + area_conc*(slab_thickness/2 + separation_between_slab_beam + profile.depth - new_centroid)**2
        Is = profile.Ix + profile.area*(new_centroid - profile.depth/2)**2
        Itr = Ic + Is

    print 'Ic = {:.1f} cm4'.format(Ic)
    print 'Is = {:.1f} cm4'.format(Is)
    print 'c = {:.2f} cm, a = {:.2f},  lever_arm = {:.1f} cm'.format(c, a,  lever_arm)
    print 'COMP = {:.1f} T, TENSION = {:.1f} T'.format(compression*1e-3, tension*1e-3)
    print 'Mn = {:.1f} Tm'.format(Mn*1e-5)
    print 'new centroid: {:.1f}'.format(new_centroid)
    print 'Itr: {:.0f} cm4'.format(Itr)
    print 'construction -- deflection (cm) per unit distributed load in T/m {:.4f} = '.format(10*5.0*length**4/384.0/profile.E/profile.Ix)
    print 'service -- deflection (cm) per unit distributed load in T/m {:.4f} = '.format(10*5.0*length**4/384./profile.E/Itr)