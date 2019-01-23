import stl_beamcolumn

def etabs_flexure(profile_object, csv_filename, output_filename, input_units='Tm'):
    """insert a etabs-generated csv file for the beam
    to be analyzed, units can be specified.
    """
    if input_units == 'Tm':
        factor = 1e-5
    else:
        print 'those units are not implemented yet'
    output = []
    file = read_csv(csv_filename + '.csv', skiprows=[1] )
    elements = file['Element'].unique()
    combinations = file['Load Case/Combo'].unique()
    for comb in combinations:
        df = file[file['Load Case/Combo'] == comb]
        for ele in elements:
            df2 = df[df['Element'] == ele]
            leng = (df2['Station'].max() - df2['Station'].min())*1e2
            # detect a change in sign for curvature
            signs = (diff(sign(df2['M3'])))*1
            moments = array([df2['M3'].iloc[0], df2['M3'].iloc[-1]])
            if sum(signs) == 0: # no cambio de signo, curvatura simple
                C_factor = 0.6 + 0.4*min(abs(moments))/max(abs(moments))
            else: # curvatura doble
                C_factor = 0.6 - 0.4*min(abs(moments))/max(abs(moments))
                if C_factor < 0.4:
                    C_factor = 0.4
            Mmax = max(abs(df['M3']))
            Mu, MR = flexure(profile_object, leng, C_factor)
            if MR*factor > abs(Mmax):
                ELU = 'Cumple'
            else:
                ELU = 'No cumple'
            output.append([comb, ele, leng, min(moments), max(moments), C_factor, Mmax, Mu*factor, MR*factor, ELU])
        
    column_dataframe = DataFrame(output, columns = ['CC', 'Tramo', 'L', 'M1', 'M2', 'C', 'Mmax', 'Mu', 'MR', 'ELU'])
    column_dataframe.to_csv(output_filename + '.csv', sep=',', encoding='utf-8')
    return column_dataframe

    

def etabs_interaction(profile, csv_filename, output_filename, input_units = 'Tm'):
    """insert a etabs-generated csv file for the column
    to be analyzed, units can be specified.
    """
    if input_units == 'Tm':
        factor = 1e-5
    else:
        print 'those units are not implemented yet'
    output = []
    file = read_csv(csv_filename + '.csv', skiprows=[1] )
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
        Rc = profile.Rc*1e-3
        Py = profile.Py*1e-3
        Pex = profile.critical_loads['pex']*1e-3
        Pey = profile.critical_loads['pey']*1e-3
        Mpx = profile.Mp[0]*factor
        Mpy = profile.Mp[1]*factor
        for extremo in extremos:
            if sum(signs) == 0: # no cambio de signo, curvatura simple
                C_factor = 0.6 + 0.4*min(abs(m3))/max(abs(m3))
            else: # curvatura doble
                C_factor = 0.6 - 0.4*min(abs(m3))/max(abs(m3))
                if C_factor < 0.4:
                    C_factor = 0.4
            MRx = flexure(profile, leng, C_factor)[1]*factor
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
                resistance = Pu/0.9/Py + 0.85*Muox/0.9/Mpx + 0.6*Muoy/0.9/Mpy
                if resistance <= 1.0:
                    ELU = 'Cumple'
                else:
                    ELU = 'No cumple'
            else:
                C_factor = 1.0
                B1x = 1./(1 - Pu/0.9/Pex)
                B1y = 1./(1 - Pu/0.9/Pey)
                Muox = '-'
                Muoy = '-'
                Muoxp = B1x*max(abs(m3))
                Muoyp = B1y*max(abs(m2))
                resistance = Pu/Rc + Muoxp/MRx + Muoyp/0.9/Mpy
                if resistance <= 1.0:
                    ELU = 'Cumple'
                else:
                    ELU = 'No cumple'
            output.append([comb, extremo, Rc, Mpx, Mpy, Pu, Muox, Muoy, B1x, B1y, C_factor, MRx,  Muoxp, Muoyp, resistance, ELU])
    column_dataframe = DataFrame(output, columns = ['CC', 'Extremo', 'Rc', 'Mpx', 'Mpy', 'Pu', 'Muox', 'Muoy', 'B1x', 'B1y', 'C', 'MRx', 'Muoxp', 'Muoyp', 'R', 'ELU'])
    column_dataframe.to_csv(output_filename + '.csv', sep=',', encoding='utf-8')
    return column_dataframe