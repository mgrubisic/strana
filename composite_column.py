import stl_beamcolumn
 
def interaction(profile, axial_load, moments=[[0, 0], [0, 0]], FR=0.9):
    """ design moments are [[Mxi, Myi], [Mxj, Myj]] and already factored by PDelta B2 coeff.
    """
    Mi = axial_load / profile.Pc + 0.80 * \
        moments[0][0] / profile.Mp[0] + 0.8 * moments[0][1] / profile.Mp[1]
    upper = Mi if Mi < FR else Mi, 'upper moments exceed capacity'
    Mj = axial_load / profile.Pc + 0.80 * \
        moments[1][0] / profile.Mp[0] + 0.8 * moments[1][1] / profile.Mp[1]
    lower = Mj if Mj < FR else Mj, 'lower moments exceed capacity'
    Mtot = axial_load / profile.Pc + \
        max(moments[0][0], moments[1][0]) / profile.Mr[0] + \
        max(moments[0][1], moments[1][1]) / profile.Mr[1]
    complete = Mtot if Mtot < FR else Mtot, 'complete column does not meet capacity'
    profile.interaction_results = {
        'lower': lower, 'upper': upper, 'complete': complete}

def composite_column_interaction(profile, concrete_column, effective_len = 100.0, k_x = 1.0,  k_y = 1.0, axial_load = 0,  etabs_csv_filename = None, output_filename = None, input_units = 'Tm'):
    """returns resistances given a profile and concrete column objects
    optionally check against etabs csv resultsfile
    """
    conc = concrete_column
    steel = profile
    Py = (conc.A - steel.A)*conc.fc + steel.A*steel.fy + conc.rebar_Pu
    fmy = steel.fy + 0.7*conc.fyr*conc.rebar_area/steel.A + 0.6*0.8*conc.fpc*conc.A/steel.A
    Em = steel.E + 0.2*conc.E*conc.A/steel.A
    slenderness = k_y*effective_len/steel.ry*(fmy/Pi**2/Em)**.5
    Rc = 0.85*steel.A*fmy/(1+slenderness**(2.8)**(1./1.4))*1e-3
    Mpx = ( (conc.Zx - steel.Zx)*conc.fc/2 + steel.Zx*steel.fy + conc.rebar_MRx )*1e-5
    Mpy = ((conc.Zy - steel.Zy)*conc.fc/2 + steel.Zy*steel.fy + conc.rebar_MRy)*1e-5
    Pemx = Pi**2*Em*steel.A/(k_x*effective_len/steel.rx)**2*1e-3
    Pemy = Pi**2*Em*steel.A/(k_y*effective_len/steel.ry)**2*1e-3
    B1x = 1./(1 - axial_load/0.9/Pemx)
    B1y = 1./(1 - axial_load/0.9/Pemy)

    if etabs_csv_filename is not None:
        output = []
        file = read_csv(etabs_csv_filename + '.csv', skiprows=[1] )
        combinations = file['Load Case/Combo'].unique()
        extremos = ['A', 'B', 'Completa']
        for comb in combinations:
            df = file[file['Load Case/Combo'] == comb]
            leng = (df['Station'].max() - df['Station'].min())*1e2
            # detect a change in sign for curvature
            signs = (diff(sign(df['M3'])))*1
            m3 = array([df['M3'].iloc[0], df['M3'].iloc[-1]])
            m2 = array([df['M2'].iloc[0], df['M2'].iloc[-1]])
            Pu = abs(df['P'].iloc[0])
            for extremo in extremos:
                if extremo == 'A':
                    Muox = abs(m3[0])
                    Muoy = abs(m2[0])
                    Muoxp = '-'
                    Muoyp = '-'
                    B1x = '-'
                    B1y = '-'
                    resistance = Pu/0.9/Py + 0.85*Muox/0.9/Mpx + 0.6*Muoy/0.9/Mpy
                    if resistance <= 1.0:
                        ELU = 'Cumple'
                    else:
                        ELU = 'No cumple'
                elif extremo == 'B':
                    Muox = abs(m3[1])
                    Muoy = abs(m2[1])
                    Muoxp = '-'
                    Muoyp = '-'
                    B1x = '-'
                    B1y = '-'
                    resistance = Pu/0.9/Py + Muox/0.9/Mpx + Muoy/0.9/Mpy
                    if resistance <= 1.0:
                        ELU = 'Cumple'
                    else:
                        ELU = 'No cumple'
                else:
                    B1x = 1./(1 - Pu/0.9/Pemx)
                    B1y = 1./(1 - Pu/0.9/Pemy)
                    Muox = '-'
                    Muoy = '-'
                    Muoxp = B1x*max(abs(m3))
                    Muoyp = B1y*max(abs(m2))
                    resistance = Pu/Rc + Muoxp/0.9/Mpx + Muoyp/0.9/Mpy
                    if resistance <= 1.0:
                        ELU = 'Cumple'
                    else:
                        ELU = 'No cumple'
                output.append([comb, extremo, Rc, Mpx, Mpy, Pu, Muox, Muoy, B1x, B1y, Muoxp, Muoyp, resistance, ELU])
    interaction_dataframe = DataFrame(output, columns = ['CC', 'Extremo', 'Rc', 'Mpx', 'Mpy', 'Pu', 'Muox', 'Muoy', 'B1x', 'B1y', 'Muoxp', 'Muoyp', 'R', 'ELU'])
    interaction_dataframe.to_csv(output_filename + '.csv', sep=',', encoding='utf-8')
    return interaction_dataframe
