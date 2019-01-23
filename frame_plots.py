"""Opensees EDP plotter"""
from pandas import read_csv, merge, read_table
from scipy.interpolate import interp1d, UnivariateSpline
from numpy import genfromtxt, deg2rad
import os
import csv
import subprocess
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
from easyplot import easyplot
home = os.path.expanduser("~")
desktop = os.path.join(home, 'Desktop')
dropbox = os.path.join(home, 'Dropbox')
cwd = os.getcwd()
print(cwd)

class frame_object(object):
    def __init__(self, absolute_path_to_framefile, storey_list, bay_list, name = None, units_length = 'm', units_mass= 'kg', units_force = 'N'):
        self.path = absolute_path_to_framefile
        self.storeys = storey_list
        self.total_height = storey_list[-1]
        self.bays = bay_list
        self.name = name
        self.units_length = units_length
        self.units_mass = units_mass
        self.units_force = units_force
        self.dofs = len(self.storeys)

    def elwood_drift(self, axial_load, transverse_steel_ratio, fy=490, spacing=0.30, depth_core = 0.60, drift=None):
        """SIe units
        kN MPa m
        """
        if drift is not None:
            self.col_ult_drift = drift
        else:
            self.col_ult_drift = 0.03 + 4*transverse_steel_ratio + ..
        return

    def eigen_analysis(self, model=1, order=2, plotmodes=2, mode_output=1, outputdir=os.getcwd()):
        # returns tuple dataframe with frequencies and dict with modes
        os.chdir(self.path + '/_output/des_output/designEigenA/')
        modelfile = 'model({})_order({}).out'.format(model, order)
        frequencies_factors = pd.read_table(modelfile, delim_whitespace=True, skiprows=5, skip_blank_lines=True, usecols=range(4), header=None, names = ['Mode', 'Period','Participation Factor', 'Gamma'])[:len(self.storeys)]

        modal_results = pd.read_table(modelfile, delim_whitespace=True, skiprows=7 + len(self.storeys), skip_blank_lines =True, usecols = range(6), header=None, names = ['Storey', 'Height','Mass', 'Displ', 'nDist', 'Force'])


        eigendict = {}
        eigendict['periods'] = frequencies_factors
        modeshapes = pd.DataFrame()
        modal_masses = []
        for mode in range(self.dofs):
            eigendict[mode] = modal_results.ix[mode*self.dofs:(mode + 1)*self.dofs - 1]
            phi = eigendict[mode]['Displ']
            mass = eigendict[mode]['Mass']
            modal_masses.append(sum(mass*phi)**2/sum(phi*mass*phi))

        eigendict['periods']['modal mass'] = modal_masses

        self.eigendict = eigendict
        self.fundamental_period = eigendict['periods']['Period'][0]
        os.chdir(outputdir)
        eigendict['periods'].to_csv('eigendata-periods-{}.csv'.format(self.name))
        for mode in range(mode_output):
            eigendict[mode].to_csv('eigendata-mode{}-{}.csv'.format(mode+1,self.name))
        fig, axes = plt.subplots(nrows=1, ncols = plotmodes, sharex=True)
        # fig.yticks = [range(self.dofs)]
        fig.suptitle('{} Modal shapes of {}'.format(plotmodes, self.name))
        for mode, axis in enumerate(axes):
            modeshape = eigendict[mode]['Displ']
            # normalized_mode = modeshape/modeshape[0]
            axis.barh(range(self.dofs), modeshape)
        fig.savefig('{}-{}mode-shapes.pdf'.format(self.name, plotmodes))

        return eigendict


    def capacity_curves(self, outputdir = os.getcwd(), ultimate_displ = 1, bilin_method = 1, elastic_slope_sensitivity = 100, postyield_slope_sensitivity=100, palette='Greys_d'):
        """ sensitivity == 2nd point for a slope in same units as model
        bilin can be
        1. equal-energy with elastic slope and initial slope
        2. equal-energy and postyield slope with origin points
        3. ?
        ultimate displ:
        1. 'rcdf' 0.03 collapse prevention
        2. 'elwood' predefined collapse displacement
        ---
        TODO: swtich to index positioning instead of displacement in mm (error in fu[du])
        """
        if ultimate_displ == 1:
            ultimate_displ = 'rcdf'
        elif ultimate_displ == 2:
            ultimate_displ = 'elwood'
        else:
            print( 'no such ultimate disp criterion is implemented')
            return
        if bilin_method == 1:
            bilin_method = 'energy-elastic slope'
        elif bilin_method == 2:
            bilin_method = 'energy-postyield slope'
        elif bilin_method == 3:
            bilin_method = 'least squares'
        else:
            print( 'no such bilin method implemented')
            return
        # sns.set_palette(sns.color_palette("BrBG", 10))
        # sns.set_palette(sns.color_palette("BrBG", 10))
        sns.set_palette(sns.color_palette(palette))
        os.chdir(self.path + '/_output/SPO_output/')

        shear_cols_pd = genfromtxt('Vbase1.out', delimiter=" ", dtype=None)
        shear_cols = genfromtxt('Vbase1_npd.out', delimiter=" " ,dtype=None)
        shearPD = np.sum(shear_cols_pd, 1)
        shear = np.sum(shear_cols, 1)
        displ = range(1, len(shear)+1)
        storey_drifts = pd.DataFrame()
        try:
            os.chdir('IS_drifts&shears1')
            for storey in range(len(self.storeys)):
                storey_drifts['storey {}'.format(storey+1)] = pd.read_table('Drift-Storey1_{}.out'.format(storey+1), names = ['storey {}'.format(storey+1)], skiprows=[-1])
        except Exception as e:
            print(e)
            print('Storey drift files do not exist or are incomplete, run SPO again with correct increments')

        self.storey_drifts = storey_drifts

        if ultimate_displ == 'rcdf':
            ult_drift = 0.03
        elif ultimate_displ == 'elwood':
            ult_drift = self.col_ult_drift
        else:
            ult_drift = ultimate_displ
        du = int(storey_drifts[storey_drifts > ult_drift].idxmax().min())
        fu = shear[du]
        fuPD = shearPD[du]
        shear = shear[:du]
        shearPD = shearPD[:du]
        displ = displ[:du]


        Ke = shear[du/elastic_slope_sensitivity]/displ[du/elastic_slope_sensitivity]
        KePD= shearPD[du/elastic_slope_sensitivity]/displ[du/elastic_slope_sensitivity]

        if bilin_method == 'energy-elastic slope':
            dy = np.argmin([abs(np.trapz(shear) - np.trapz([0, Ke*d, fu],[0,d,du])) for d in range(du)])
            dyPD = np.argmin([abs(np.trapz(shearPD) - np.trapz([0, KePD*d, fuPD],[0, d, du])) for d in range(du)])
            fy = Ke*dy
            fyPD = KePD*dyPD
            Kt = (fu-fy)/(du-dy)
            KtPD = (fuPD - fyPD)/(du - dyPD)

        elif bilin_method == 'energy-postyield slope':
            Kt = (fu - shear[-postyield_slope_sensitivity])/(du - displ[-postyield_slope_sensitivity])
            KtPD = (fuPD - shearPD[-postyield_slope_sensitivity])/(du - displ[-postyield_slope_sensitivity])
            dy = np.argmin([abs(np.trapz(shear) - np.trapz([0, fu - Kt*(du - d) ,fu], [0, d, du])) for d in range(du)])
            dyPD = np.argmin([abs(np.trapz(shearPD) - np.trapz([0, fuPD - KtPD*(du - d) ,fuPD], [0, d, du])) for d in range(du)])
            fy = fu - Kt*(du - dy)
            fyPD = fuPD - KtPD*(du - dy)
            # dy = np.argmin([np.abs(np.trapz(shear, dx = du/len(displ)) - np.trapz([0, Ke*d, fu],[0, d, du], du/len(displ))) for d in np.arange(0,du,du/len(displ))])
            # dyPD = np.argmin([np.abs(np.trapz(shearPD, dx= du/len(displ)) - np.trapz([0, KePD*d, fuPD],[0, d, du], du/len(displ))) for d in np.arange(0,du,du/len(displ))])

        elif bilin_method == 'least squares':
            """where elastic slope stops touching, where shear stops being linear"""
            di = elastic_slope_sensitivity
            dy = np.argmin([abs(np.trapz(shear) - np.trapz([0, Ke*d, fu],[0,d,du])) for d in range(du)])
            fy = Ke*dy
            Kt = (fu-fy)/(du-dy)

            KtPD, b = np.polyfit(displ[di:], shearPD[di:], 1)
            dyPD = int( b/(Ke - KtPD) )
            fyPD = shear[dyPD]

        else:
            print( 'specified bilinearization method is not yet implemented')

        ductility = float(du)/dy
        ductilityPD = float(du)/dyPD
        theta = (Kt - KtPD)/Ke
        alpha = Kt/Ke
        thetaE = (Ke - KePD)/Ke
        thetaI = (Kt - KtPD)/Ke

        Taux = self.fundamental_period*((1-alpha)/(1-alpha-thetaE + thetaI))**0.5
        thetaaux = (thetaI - alpha*thetaE)/(1 - alpha - thetaE + thetaI)


        yield_drift = float(dy)/self.total_height
        yield_drift_PD = float(dyPD)/self.total_height
        ult_drift = float(du)/self.total_height

        plt.figure()
        fig, ax = plt.subplots()
        rows = ['Normal', r'$P\Delta$']
        cols = [r'$\Delta_y$',r'$\gamma_y$' ,r'$V_y$', r'$\Delta_u$', r'$\gamma_u$',r'$V_u$', r'$\mu$']
        data = [[dy, '{:.3f}'.format(yield_drift), "{:.1f}".format(fy), du,'{:.2f}'.format(ult_drift), "{:.2e}".format(fu), "{:.1f}".format(ductility)],[dyPD, '{:.3f}'.format(yield_drift_PD), "{:.2e}".format(fyPD), du , "{:.3f}".format(ult_drift), "{:.2e}".format(fuPD), "{:.1f}".format(ductilityPD)]]
        ax.table(cellText = data, rowLabels = rows, colLabels = cols, loc= 'top')
        plt.plot(shear, linestyle= '-' )
        plt.plot(shearPD,  linestyle= '-')
        # plt.plot([0,0.5,2],[0,0.5*Ke,shear[-1]])
        bilin = plt.plot([0, dy, du], [0, fy, fu], linestyle = '--', label="$T_1$ = {:.3f} s \n".format(self.fundamental_period) + r"$\alpha$ = {:.3f}".format(alpha) + '\n' +  r"$\theta_A - \alpha$ = {:.3f}".format(thetaaux-alpha) +  "\n $K_E$ = {:.1f} \n $K_I$ = {:.1f}".format(Ke ,Kt) )
        bilinPD = plt.plot([0, dyPD, du], [0, fyPD, fuPD],linestyle = '--',label = "$K_I'$ = {:.1f} \n".format(KtPD) + r"$\theta_E$ = {:.3f}".format(thetaE) + "\n" + r"$\theta_I$ = {:.3f} ".format(thetaI) + "\n" + r"$ T_A$ = {:.3f} s".format(Taux) +"\n" + r"$\theta_A$ = {:.3f}".format(thetaaux))
        plt.xlabel('Roof displacement ({})'.format(self.units_length))
        plt.ylabel('Shear ({})'.format(self.units_force))
        plt.title('Bilinear cap. curves of a {} frame bilin-{}, du : {}'.format(self.name, bilin_method, ultimate_displ), y = 0.92)
        ax.get_xaxis().set_minor_locator(ticker.AutoMinorLocator())
        ax.get_yaxis().set_minor_locator(ticker.AutoMinorLocator())
        ax.grid(b=True, which='major', color='w', linewidth=1.0)
        ax.grid(b=True, which='minor', color='w', linewidth=0.5)
        legend = ax.legend(loc='lower right', shadow=True)
        frame = legend.get_frame()

        os.chdir(outputdir)
        plt.savefig('cap-curves-{}-bilin-{}-du-{}.png'.format(self.name, bilin_method, ultimate_displ))
        plt.clf

        ccdict = {}
        ccdict['bilin method'] = bilin_method
        ccdict['ult displ'] = ultimate_displ
        ccdict['yield displ'] = dyPD
        ccdict['ductility'] = ductilityPD
        ccdict['yield shear'] = fyPD
        ccdict['ult shear'] = fuPD
        ccdict['elastic slope'] = KePD
        ccdict['postyield slope'] = KtPD
        ccdict['stability coeff'] = theta
        ccdict['elastic stab coeff'] = thetaE
        ccdict['inelastic stab coeff'] = thetaI
        ccdict['auxiliar stab coeff'] = thetaaux
        ccdict['auxiliar period'] = Taux
        ccdict['alpha'] = alpha
        ccdict['axial load'] = theta - alpha
        ccdict['auxiliar axial load'] = thetaaux - alpha

        self.ccdict = ccdict
        writer = csv.writer(open("capacity-data-{}-{}-{}.csv".format(self.name, ultimate_displ, bilin_method), "w"))
        for key, value in self.ccdict.items():
            writer.writerow([key, value])

        return ccdict


    def idaplot(self, outputdir = os.getcwd(), records = range(1), collapse_im = 3.0, collapse_drift=0.08,  EDP = 'Interstorey drift (rad)', IM = r'$A(T_1)/g$', color='blue'):
        # can be run in IDA_output directory
        """ requires a set of idas to be run
        returns idadict['records' = 1, 2, ..] = DataFrame([edp, true edp]), index = IM
        """
        if color is 'blue':
            sns.set_palette(sns.color_palette("Blues_d"))
        elif color is 'red':
            sns.set_palette(sns.color_palette('Reds_d'))
        else:
            print( 'that color palette is not yet implemented')

        idadict = {}
        intensities = []
        fig = plt.figure()
        # axes = fig.add_subplot(1,1,1)
        # fig, axes = plt.subplots(111)
        plt.xlabel(EDP)
        plt.ylabel(IM)
        plt.xlim(0, collapse_drift)
        plt.ylim(0, collapse_im)
        for record in records:
            os.chdir(self.path + '/_output/IDA_output/')
            # namelist=['intensity','displ env','Vbenv','Dmin','Vb(Dmin)','Dmax','Vb(Dmax)','D(Vbmin)','Vbmin','D(Vbmax)','Vbmax','Sa','Sa/g']
            namelist = ['intensity']+['storey {}'.format(st+1) for st in range(len(self.storeys))]

            envelopes_df = pd.read_table('record({})/IDA_EDRIFT_rec({}).out'.format(record+1 ,record+1), delim_whitespace=True, skiprows=0,  header=None, names=namelist)
            # envelopes_df = pd.read_table('record({})/IDA_ENV_rec({}).out'.format(record+1,record+1), delim_whitespace = True, skiprows = 1, index_col=0, usecols=range(1), header=None, names=['Denv'])
            for ix, st in enumerate(self.storeys):
                envelopes_df['storey {}'.format(ix+1)] = envelopes_df['storey {}'.format(ix+1)]/st

            envelopes_df['max drift'] = envelopes_df.iloc[:, 1:].max(axis=1)
            # envelopes_df['stable drift'] = envelopes_df['drift'].apply(lambda drift: 'NaN' if drift > collapse_drift else drift)
            # plt.plot(envelopes_df['stable drift'], envelopes_df['Sa/g'])
            # envelopes_df['smooth IM'] = envelopes_df['Int(PGA)'].interpolate('spline', order=3)
            # , 'base shear envelope', 'disp min', 'base shear (disp min)'])
            #x = c12['stable drift'].tolist(0,0)
            x = envelopes_df['max drift'].tolist()
            y = envelopes_df['intensity'].tolist()
            x.insert(0,0)
            y.insert(0,0)
            plt.plot(x,y, linewidth=1, alpha=0.7 )
            # xnew = np.linspace(0, x[-1])
            # print( len(x))
            # f1 = interp1d(x,y, kind='linear', bounds_error=False, fill_value = 'extrapolate')
            # f2 = interp1d(x,y, kind='linear', bounds_error=False, fill_value='extrapolate')
            # f3 = interp1d(x,y, kind='slinear')
            # f4 = interp1d(x,y, kind='quadratic')
            # try:
            # f2 = interp1d(x, y, kind='cubic')
            # plt.plot(xnew, f2(xnew), '-')
            # except Exception as e:
            #     print( e)
            #     print( 'analysis has too few points (0 or 1) or plot failed')
            # else:
            # plt.plot(x, y, '--')
            # plt.plot(x, y,'o', xnew, f2(xnew), '-', xnew, f3(xnew), 'b-', xnew, f4(xnew), '--')
            # plt.plot(xnew, f2(xnew))
            # IM = envelopes_df['Intensity']
            # EDP = envelopes_df['Displ Env']/self.storeys[0]
            # envelopes_df['edp'] = EDP
            intensities.append(max(y))
            idadict['record {}'.format(record+1)] = pd.DataFrame()
            idadict['spline record {}'.format(record+1)] = pd.DataFrame()
            idadict['record {}'.format(record+1)]['max drift'] = envelopes_df['max drift']
            idadict['record {}'.format(record+1)]['Sa/g'] = envelopes_df['intensity']
            # idadict['spline record {}'.format(record+1)]['max drift'] = xnew
            # idadict['spline record {}'.format(record+1)]['Sa/g'] = f2(xnew)
            # idadict[record]['max drift'] = envelopes_df['max drift']
            # idadict[record]['spline max drift'] =
            # idadict[record]['im'] = envelopes_df['Sa/g']
            # idadict[record]['im'] = IM
            # idadict[record]['true edp'] = idadict[record]['edp'].apply(lambda disp: 0 if disp > collapse_max drift else disp)
            # print( idadict[record])
            # trueEDP = idadict[record]['true edp']
            # envelopes_df.iloc[:,1:3]. plot(label = 'record {}'.format(record))
            # envelopes_df['stable Drift'].T.plot(label = 'record {}'.format(record))
            # plt.plot(envelopes_df['Drift'], envelopes_df['Int'], 'b--', envelopes_df['stable Drift'], envelopes_df['Int'], 'g-')
            # plt.plot(envelopes_df['stable Drift'], envelopes_df['smooth IM'])
        idadict['intensties'] = intensities
        idadict['median cc'] = np.median(intensities)
        os.chdir(outputdir)
        writer = csv.writer(open("idadict-{}.csv".format(self.name), "w"))
        writer.writerow(['median cc', idadict['median cc']])
        writer.writerow(['intensities Sa/g', intensities])
        # for key, value in idadict.items():
        #     writer.writerow([key, value])
        self.idadict = idadict
        plt.title(self.name + ' IDA curves for {} record(s), Median CC = {:.2f} g'.format(len(records), np.median(intensities)))
        plt.savefig(self.name + '-IDA-records(1-{})'.format(len(records)) + '.png')
        plt.clf()
        return idadict


    def dmg_evol(self, outputdir, ida_records, sensitivity = 0.05, lineweight = 2, dmg_size = 10, fontsize=5, ncols=10):
        dmg_evolution(self.bays, self.storeys, self.name, self.path, outputdir, ida_records, sensitivity, lineweight, dmg_size, fontsize, ncols)

    def dmg_distr(self, alpha, outputdir):
        dmgplot(self.bays, self.storeys, self.path + '/_output/', outputdir = outputdir, analysis = 'SPO', alpha = alpha)

    def collapse_assessment(self, outputdir=os.getcwd()):
        """outputs csv with collapse data computation"""
        try:
            modal_mass_mode1 = self.eigendict['periods']['modal mass'][0]
        except:
            print( 'no eigendict found, run eigen_analysis() on the frame first')
        

        self.ccdict['modal mass 1'] = modal_mass_mode1
        self.ccdict['Say'] = self.ccdict['yield shear']/modal_mass_mode1

        os.chdir(outputdir)
        writer = csv.writer(open("ccdict-{}.csv".format(self.name), "w"))
        for key, value in self.ccdict.items():
            writer.writerow([key, value])

        return 

    def cc_spectrum(self, spectradir):
        """computes median collapse capacity given a eigen_analysis, capacity_curve_model, collapse_assessment
        requires yield str spectra to be in spectra dir
        returns self.ccdict
        """
        ductilities = [1.00*mu for mu in range(1,11)]
        axial_loads = [0, 0.025, 0.05, 0.075, 0.1, 0.125, 0.150, 0.175, 0.2]
        try:
            duct = min(range(len(ductilities)), key=lambda i: abs(ductilities[a] - self.ccdict['ductility']))
        except Exception as e:
            print( e)
            print( 'please run a capacity model first')
            print( 'make sure a ductility is defined')

        try:
            axial_load = min(range(len(axial_loads)), key=lambda i: abs(axial_loads[a] - self.ccdict['axial load']))
        except Exception as e:
            print( e)
            print( 'please run a capacity model first')
            print( 'make sure an axial load is defined')

        self.ccdict['model ductility'] = duct
        self.ccdict['model axial load'] = axial_load
        os.chdir(spectradir)
        try:
            r_spectra = pd.read_table('R_th({})_mu({:.0f}).txt'.format(axial_load, duct*1000), delim_whitespace=True, index_col=0, names=['T', 'R'])
        except Exception as e:
            print( e)
            print( 'no such spectra found, ensure axial loads and ductility are correct')
            print( 'mu should be in thousands, th in float with specific )significant figures')

            self.R_median_spectra = r_spectra['R'].ix[round(self.fundamental_period, 1)]
            self.Sae_collapse = self.R_median_spectra*self.ccdict['Say']

        return


    def collapse_eval():
        """
        requires: eigen, cap_curves, idaplot, ccspectrum
        returns: percentile values and comparisons of cc
        """
        return




def dmg_evolution(bays, storeys, framename, input_path, outputdir, IDA_record=1, sensitivity = 0.05, grid_lineweight = 2, dmg_size = 10, fontsize=5, ncols = 10):

    namelist=['intensity','displ env','Vbenv','Dmin','Vb(Dmin)','Dmax','Vb(Dmax)','D(Vbmin)','Vbmin','D(Vbmax)','Vbmax','Sa','Sa/g']
    for record in range(IDA_record):
        os.chdir(input_path + '/_output/IDA_output/')
        envelopes_df = pd.read_table('record({})/IDA_ENV_rec({}).out'.format(record+1,record+1), delim_whitespace=True, skiprows=1,  header=None, names=namelist)
        directory = os.chdir(input_path + '/_output/IDA_output/' + 'record({})/TH/'.format(record + 1))
        intensities = envelopes_df['intensity']

        nrows = (len(intensities) - 1)//ncols + 1
        fig, axes = plt.subplots(nrows, ncols)
        # , sharex=True, sharey=True)
        fig.suptitle('{} Damage evolution for record {}'.format(framename, record + 1))
        # for indexing purposes
        axes = axes.flatten()

        bays = [0] + bays
        storeys = [0] + storeys
        framegrid = [(x,y) for y in storeys for x in bays]


        for intensity, value in enumerate(intensities):
            # plt.subplot(nrows, ncols, intensity)
            axis = axes[intensity]
            logfile = 'PH_log_rec({})_SF({}).out'.format(record + 1, intensity + 1)
            try:
                phlogs = pd.read_table(logfile, delim_whitespace = True, names = ['step','step','nodeID','type','end'], skiprows = 3)[['nodeID','type','end']]
            except IOError:
                print( logfile + 'does not exist.')
                continue
            axis.axis('scaled')
            X = axis.set_xticks(bays)
            Y = axis.set_yticks(storeys[1:])
            xoffset = bays[-1]*sensitivity
            yoffset = bays[-1]*sensitivity*0.8
            xlim = axis.set_xlim([-xoffset, bays[-1] + xoffset])
            ylim = axis.set_ylim([-yoffset, storeys[-1] + yoffset])
            xticks = axis.set_yticklabels(range(1,len(storeys)), fontsize=fontsize)
            yticks = axis.set_xticklabels(range(len(bays)), fontsize=fontsize)
            axis.grid(True, lw = grid_lineweight, color = 'black', alpha=0.5)


            for index, node in enumerate(phlogs['nodeID']):
                elem = phlogs['type'][index]
                end = phlogs['end'][index]
                if elem == 'column':
                    x = framegrid[node - 1][0]
                    if end == 'j':
                        y = framegrid[node + len(bays) - 1][1] - yoffset
                    else:
                        y = framegrid[node - 1][1] + yoffset
                elif elem == 'beam':
                    # this transf is necessary since we switch from 4 to 3 elements
                    # and want to start at the first storey not the floor
                    # warning! : integer arithmetic is needed so no __future__ division imports
                    beam_coord_transform = node -1 + 2*len(bays) - len(framegrid)
                    beamskips = (node -1 + len(bays) - len(framegrid))/3
                    # print( type(beam_coord_transform), type(beamskips))
                    # print( 'actual nodeID {}'.format(beam_coord_transform + 1))
                    # print( 'floor skips {}'.format(beamskips))
                    y = framegrid[beam_coord_transform + beamskips][1]
                    if end =='j':
                        x = framegrid[beam_coord_transform + 1 + beamskips][0] - xoffset
                    else:
                        x = framegrid[beam_coord_transform + beamskips][0] + xoffset
                else:
                    print( 'no such element, check input file for mispelled name')
                # else:
                #     print( 'no such element found, check input plastic hinge )log file for mispelled characters'

                    # print( elem, end, x, y)
                    # print( index, node)
                axis.scatter(x, y, color = 'r' if elem == 'column' else 'b', alpha=0.8, s = dmg_size if elem=='column' else dmg_size/2)
                axis.set_title(r'$S_a = {:.3f}$ g'.format(value), fontsize=fontsize)
        os.chdir(outputdir)
        fig.savefig('dmg-evol-{}-record({}).pdf'.format(framename, record + 1))
        fig.clf()
    return

# axes.xaxis.set_major_locator(ticker.FixedLocator(X))
# axes.xaxis.set_major_formatter(ticker.FixedFormatter(( [x for x in range(4)] )))
# axes.yaxis.set_major_locator(ticker.FixedLocator(Y))
# axes.yaxis.set_major_formatter(ticker.FixedFormatter(( [x for x in range(20)] )))

def dmgplot(bays, storeys, input_path = '../_output/', outputdir = os.getcwd(), analysis = 'SPO', IDA_record = 1, IDA_SF = 1, alpha = None, sensitivity = 0.05):
    # bays = [x1, x2 ... xn]
    # storeys = [y1, y2 ... ym]
    # filepath = '/Users/username/.../_output/'
    # alpha =0.00, 0.03
    if analysis == 'SPO':
        logfile = 'PH_log1.out'
        if alpha is not None:
            directory = os.chdir(input_path + '{}_output_alpha={}/'.format(analysis, alpha))
        else:
            directory = os.chdir(input_path + '{}_output/'.format(analysis))


    elif analysis == 'IDA':
        logfile = 'PH_log_rec({})_SF({}).out'.format(IDA_record, IDA_SF)
        directory = os.chdir(input_path + '{}_output/'.format(analysis) + 'record({})/TH/'.format(IDA_record))
    else:
        print( 'no such analysis, check input analysis value')
    print( directory)
    print( logfile)
    phlogs = pd.read_table(logfile, delim_whitespace = True, names = ['step','step','nodeID','type','end'], skiprows = 3)[['nodeID','type','end']]
    print( phlogs)
    fig = plt.figure()
    axes = fig.add_subplot(1,1,1)
    plt.axis('scaled')
    X = axes.set_xticks([0]+bays)
    Y = axes.set_yticks(storeys)
    xoffset = bays[-1]*sensitivity
    yoffset = bays[-1]*sensitivity*0.8
    axes.set_xlim([-xoffset, bays[-1] + xoffset])
    # axes.set_xlim(-1000,22000)
    axes.set_ylim([-yoffset, storeys[-1] + yoffset])
    axes.set_yticklabels([x+1 for x in range(len(storeys))])
    axes.set_xticklabels([x for x in range(len(bays) + 1)])
    plt.grid(True, lw = 6, color = 'black', alpha=0.5)

    bays = [0] + bays
    storeys = [0] + storeys
    framegrid = [(x,y) for y in storeys for x in bays]

    for index, node in enumerate( phlogs['nodeID'] ):
        elem = phlogs['type'][index]
        end = phlogs['end'][index]
        print( elem)
        if elem == 'column':
            x = framegrid[node - 1][0]
            if end == 'j':
                y = framegrid[node + len(bays) - 1][1] - yoffset
            else:
                y = framegrid[node - 1][1] + yoffset
        elif elem == 'beam':
            # this transf is necessary since we switch from 4 to 3 elements
            # and want to start at the first storey not the floor
            # warning! : integer arithmetic is needed so no __future__ division imports
            beam_coord_transform = node -1 + 2*len(bays) - len(framegrid)
            beamskips = (node -1 + len(bays) - len(framegrid))/3
            print( 'actual nodeID {}'.format(beam_coord_transform + 1))
            print( 'floor skips {}'.format(beamskips))
            y = framegrid[beam_coord_transform + beamskips][1]
            if end =='j':
                x = framegrid[beam_coord_transform + 1 + beamskips][0] - xoffset
            else:
                x = framegrid[beam_coord_transform + beamskips][0] + xoffset
        else:
            print( 'no such element found, check input ph log file for )mispelled characters')
        plt.scatter(x,y, color = 'r' if elem == 'column' else 'b', alpha=0.6, s = 20 if elem=='column' else 10)
    plt.title('Damage distribution of a {} bay {} storey frame in {}'.format(len(bays) - 1 ,len(storeys) - 1, analysis))
    os.chdir(outputdir)
    plt.savefig('dmg-distr-{}st{}b-{}-SF{}.png'.format(len(storeys) - 1, len(bays) - 1, analysis, IDA_SF))
    plt.clf()
    return framegrid


def easyplot(arraylike_x, arraylike_y, xlabel = 'X', ylabel = 'Y', title = 'my plot', save_as = False, scatter=False, linestyle = '-', color = 'black', grid = True, scaled=False):
# , legend = None, xlim = None, ylim = None, xticks = None, yticks = None, scaled = False, scatter = False):
    fig = plt.figure()
    axes = fig.add_subplot(1,1,1)
    if scaled == True:
        plt.axis('scaled')
    # axes.set_xticks(xticks)
    # axes.set_yticks(yticks)
    # axes.set_xlim(xlim)
    # axes.set_ylim(ylim)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    # plt.grid(grid)
    plt.title(title)
    if scatter == True:
        plt.scatter(arraylike_x, arraylike_y, color = color, s = 20)
    else:
        plt.plot(arraylike_x, arraylike_y, color = color, ls = linestyle)
    if save_as is not False:
        plt.savefig(save_as + '.png')
    return fig
