"""Opensees EDP plotter"""
from pandas import read_csv, merge, read_table, DataFrame
from scipy.interpolate import interp1d, interp2d, UnivariateSpline
from numpy import genfromtxt, deg2rad, trapz, argmin, linspace, array, zeros
import os
import csv
import subprocess
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
from easyplot import easyplot
sns.set()
home = os.path.expanduser("~")
desktop = os.path.join(home, 'Desktop')
dropbox = os.path.join(home, 'Dropbox')
cwd = os.getcwd()
print('Current working dir is '+cwd)


class Response(object):
    """generalized force-deformation F-u plot
    """
    def __init__(self, filename, input_folder = os.getcwd(),  output_folder=os.getcwd()):
        self.filename = filename
        self.input_folder = input_folder
        self.output_folder = output_folder
    # def FDEplot(self, force_filename, disp_filename, output_path=os.getcwd()):
    #     os.chdir(output_path)
    #     force = f
    def run_opensees(self, tcl_folder_location = os.getcwd()):
        os.chdir(tcl_folder_location)
        subprocess.call("opensees " + self.filename + ".tcl" , shell=True)
        os.chdir(cwd)

    def time_history(self, file_name, column_names = ['time', 'force'], save_folder = os.getcwd(), save_name = 'generic_th_plot.pdf'):
        """ plots a generic | time-force relationship, TIME should always be on the first column
        """
        os.chdir(self.input_folder)
        fig = plt.figure()
        df1 = read_table(file_name, delim_whitespace = True, names = column_names)
        # print(df1)
        df1.plot(x='time', y=column_names[1:])
        os.chdir(self.input_folder + save_folder)
        plt.savefig(save_name)
        plt.close(fig)
        os.chdir(cwd)
        return df1

    def moment_curvature_plot(self, curvature_filename, moments_filename, savename = 'moment-rotation.pdf', sign=-1):
        fig = plt.figure()
        deformation = read_csv(curvature_filename, names = ['strain', 'curvature'], sep = ' ')

        momentframe = read_csv(moments_filename, names = ['axial load', 'moment'], sep = ' ')
        dataframe = merge(deformation,momentframe, left_index=True, right_index=True)
        dataframe['moment'] = dataframe['moment']*sign
        # print( dataframe)

        dataframe.plot(x='curvature', y='moment')
        os.chdir(self.output_folder)
        plt.savefig(savename)
        plt.close(fig)
        os.chdir('../')#?????


    def mode_shapes_plot(self, mode_shapes_filename, savename='mode_shapes'):
        """ plots eigenvectors savename directory
        """
        os.chdir(self.input_folder)
        df = read_table(mode_shapes_filename, header=None, delim_whitespace = True)
        df = df.transpose()
        dofs = df.shape[0]
        zeroes = DataFrame(zeros((1, dofs)))
        df = zeroes.append(df)
        df.columns = ['dof ' + str(i) for i in range(dofs)]
        fig, ax = plt.subplots(1, dofs, sharey=True)
        fig.subplots_adjust(wspace=0)
        fig.suptitle('Modal shapes for ' + self.filename)
        for i in range(dofs):
            ax[i].plot(df['dof ' + str(i)], [x for x in range(dofs+1)])
            ax[i].set_title('Mode {}'.format(i+1))

        plt.savefig(savename + '.pdf')



    def capacity_curve_plot(self, drift_filename, column_forces_filename, savename='capacity_curve', normalize_weight = 1.0, normalize_drifts = 1.0, collapse_drift = 0.15, period = 1.0, elastic_slope_sensitivity = 2, drift_tol=1.0e-4, PDelta=False, delim_whitespace = True):
        """plots both capacity curves (Pdelta and linear) and approximates it with a
        piecewise bilinear function
        """

        os.chdir(self.input_folder)
        fig = plt.figure()
        # dataframe = read_csv(drift_filename, names = ['drift'])
        df = read_table(column_forces_filename, header=None, delim_whitespace = delim_whitespace)
        df['Regular'] = - 1.0/normalize_weight * df.sum(1)
        df['drift'] = read_table(drift_filename, names = ['drift'])/normalize_drifts

        if PDelta is False:
            # ax = plt.savefig(savename +'.pdf')
            ax = df.plot('drift', 'Regular')
            ax.set_xlabel(r"drift %")
            ax.set_ylabel(r"shear/weight %")
        else:
            df2 = read_table(PDelta, header=None, delim_whitespace = delim_whitespace)
            df['PDelta'] = -1.0/normalize_weight * df2.sum(1)
            print(df['PDelta'])
            ax = df.plot(x='drift', y = ["Regular", "PDelta"])
            ax.set_xlabel(r"distorsión de azotea")
            ax.set_ylabel(r"cortante/peso")
            # plt.savefig(savename +'_PDelta.pdf')

        drift = df.drift
        shear = df.Regular
        shearPD = df.PDelta
        du = collapse_drift

        try:
            du_index = drift[abs(drift - du) <= drift_tol].index[0]
        except IndexError as err:
            du_index = 0
            print("Specify a less strict drift tolerance (1.0e-3 recommended)", err)

        fu = shear.iloc[du_index]
        fuPD = shearPD.iloc[du_index]

        Ke = shear.iloc[0 + elastic_slope_sensitivity]/drift.iloc[0 + elastic_slope_sensitivity]

        KePD = shearPD.iloc[0 + elastic_slope_sensitivity]/drift.iloc[0 + elastic_slope_sensitivity]

        print(Ke, KePD, fu, fuPD)
        dy_index = argmin([abs(trapz([0, Ke*d, fu], [0, d, du]) - trapz(shear,drift)) for d in linspace(0, du, 100)])
        dy = linspace(0, du, 100)[dy_index]

        print("dy = {:.3f}".format(dy))

        fy = Ke*dy
        dyPD_index = argmin([abs(trapz([0, KePD*d, fuPD], [0, d, du]) - trapz(shearPD,drift)) for d in linspace(0, du, 100)])
        dyPD = linspace(0, du, 100)[dyPD_index]
        fyPD = KePD*dyPD
        self.fyPD = fyPD


        Ki = (fu-fy)/(du-dy)
        KiPD = (fuPD - fyPD)/(du - dy)
        ductility = du/dy
        ductilityPD = du/dyPD
        alpha = Ki/Ke
        theta_e = (Ke - KePD)/Ke
        theta_i = (Ki - KiPD)/Ke # this is what can be wrong
        period_aux = period*((1 - alpha)/(1 - alpha - theta_e + theta_i))**0.5
        theta_aux = (theta_i - alpha*theta_e)/(1 - alpha  - theta_e + theta_i)

        self.period = period_aux
        self.axial_load = theta_aux - alpha
        self.Ay = self.fyPD*9.81

        bilin = plt.plot([0, dy, du], [0, fy, fu], linestyle = '--', label="bilin. Regular")
        bilinPD = plt.plot([0, dyPD, du], [0, fyPD, fuPD],linestyle = '--',label ="bilin. PDelta")
        # bilin = plt.plot([0, dy, du], [0, fy, fu], linestyle = '--', label="$T_1$ = {:.3f} s \n".format(period) + r"$\alpha$ = {:.3f}".format(alpha) + '\n' +  r"$\theta_a - \alpha$ = {:.3f}".format(theta_aux -alpha))
        # bilinPD = plt.plot([0, dyPD, du], [0, fyPD, fuPD],linestyle = '--',label = r"$\theta_e$ = {:.3f}".format(theta_e) + "\n" + r"$\theta_i$ = {:.3f} ".format(theta_i) + "\n" + r"$ T_A$ = {:.3f} s".format(period_aux) +"\n" + r"$\theta_a$ = {:.3f}".format(theta_aux))

        legend = ax.legend(loc='lower center', shadow=True)
        frame = legend.get_frame()

        plt.xlim(0, 0.1)
        plt.savefig(savename+"_bilin"+".pdf")
        self.fig = fig
        plt.close(fig)
        plt.close("all")

    def collapse_capacity(self, savename, spectra_file, axial_loads=[0, 0.025, 0.05, 0.075, 0.10], extrapolate=[0.125, 0.150]):

        df = read_table(spectra_file, sep=",", skiprows=1, names = ["period"] + axial_loads)

        for axial_load in extrapolate:
            df.insert(len(df.columns), axial_load, "NaN")

        original_values = len(axial_loads)
        new = len(extrapolate)
        print(len(df.columns)-new, original_values)

        for row in range(len(df)):
            x = df.iloc[row, 1:original_values + 1].index
            y = df.iloc[row, 1:original_values + 1].values
            f = interp1d(x,y, fill_value = "extrapolate")
            df.iloc[row, original_values+1:] = f(extrapolate)

        periods = df.period
        axials = axial_loads + extrapolate
        func = interp2d(axials, periods, df.iloc[:, 1:])
        Ry = func(self.axial_load, self.period)
        cc = self.Ay*Ry
        self.cc = cc
        self.fig.suptitle(r"$T_A$ = {:.3f} s, $\theta_a - \alpha$ = {:.3f}, ".format(self.period, self.axial_load)+r"$R_y$ = {:.2f}, $A_y$ = {:.2f} g, $CC$ = {:.2f} g".format(Ry[0], self.Ay, cc[0]) )
        self.fig.savefig(savename + "_capacities"+".pdf")

        # df = df.iloc[:15]
        # df.plot(x="period", y=[0.10, 0.150])
        # plt.show()
        plt.close("all")
        plt.figure()
        ax2 = df.plot(x="period", y=axial_loads)
        ax2.set_xlabel("periodo [s]")
        ax2.set_ylabel(r"$A/A_y$")
        ax2.annotate(r"$\xi = 5$ %", (0.05, 0.45), textcoords='axes fraction', size=14)
        plt.savefig("collapse_capacity_spectra_soft_soil2.png")
        plt.savefig("collapse_capacity_spectra_soft_soil2.pdf")

    def idaplot(self, ida_output_dir, output_names = ['ida_performance_rec_' + str(rec+1) + '.out' for rec in range(2, )], collapse_drift= 0.15, collapse_im = 2.0):
        """ Plots ida results in a nice format
        requires all output files to have a similar 'output_name' and numbered 1..n
        inside 'ida_output_dir'

        RETURNS the median collapse intensity and the median collapse drift
        """
        ida_dict = {}
        sns.set_palette(sns.color_palette("Blues_d"))
        # sns.set_palette(sns.color_palette('Reds_d'))
        fig = plt.figure()
        plt.xlabel('Interstorey drift (rad)')
        plt.ylabel('Spectral first-mode pseudo-acceleration (g)')
        plt.xlim(0, collapse_drift)
        plt.ylim(0, collapse_im)
        name_list = ['drift', 'intensity']
        os.chdir(self.input_folder + ida_output_dir)

        for name in output_names:
            performance_record = read_table(name, delim_whitespace=True, skiprows=1, header=None, names=name_list)
            print(performance_record)
            x = performance_record['drift'].tolist()
            y = performance_record['intensity'].tolist()
            x.insert(0,0)
            y.insert(0,0)
            plt.plot(x, y , linewidth=0.7, alpha=0.7, marker='.')

        #     # xnew = np.linspace(0, x[-1])
        #     # print( len(x))
        #     # f1 = interp1d(x,y, kind='linear', bounds_error=False, fill_value = 'extrapolate')
        #     # f2 = interp1d(x,y, kind='linear', bounds_error=False, fill_value='extrapolate')
        #     # f3 = interp1d(x,y, kind='slinear')
        #     # f4 = interp1d(x,y, kind='quadratic')
        #     # try:
        #     # f2 = interp1d(x, y, kind='cubic')
        #     # plt.plot(xnew, f2(xnew), '-')
        #     # except Exception as e:
        #     #     print( e)
        #     #     print( 'analysis has too few points (0 or 1) or plot failed')
        #     # else:
        #     # plt.plot(x, y, '--')
        #     # plt.plot(x, y,'o', xnew, f2(xnew), '-', xnew, f3(xnew), 'b-', xnew, f4(xnew), '--')
        #     # plt.plot(xnew, f2(xnew))
        #     # IM = performance_record['Intensity']
        #     # EDP = performance_record['Displ Env']/self.storeys[0]
        #     # performance_record['edp'] = EDP
        #     intensities.append(max(y))
        #     idadict['record {}'.format(record+1)] = pd.DataFrame()
        #     idadict['spline record {}'.format(record+1)] = pd.DataFrame()
        #     idadict['record {}'.format(record+1)]['max drift'] = performance_record['max drift']
        #     idadict['record {}'.format(record+1)]['Sa/g'] = performance_record['intensity']
        #     # idadict['spline record {}'.format(record+1)]['max drift'] = xnew
        #     # idadict['spline record {}'.format(record+1)]['Sa/g'] = f2(xnew)
        #     # idadict[record]['max drift'] = performance_record['max drift']
        #     # idadict[record]['spline max drift'] =
        #     # idadict[record]['im'] = performance_record['Sa/g']
        #     # idadict[record]['im'] = IM
        #     # idadict[record]['true edp'] = idadict[record]['edp'].apply(lambda disp: 0 if disp > collapse_max drift else disp)
        #     # print( idadict[record])
        #     # trueEDP = idadict[record]['true edp']
        #     # performance_record.iloc[:,1:3]. plot(label = 'record {}'.format(record))
        #     # performance_record['stable Drift'].T.plot(label = 'record {}'.format(record))
        #     # plt.plot(performance_record['Drift'], performance_record['Int'], 'b--', performance_record['stable Drift'], performance_record['Int'], 'g-')
        #     # plt.plot(performance_record['stable Drift'], performance_record['smooth IM'])
        # idadict['intensties'] = intensities
        # idadict['median cc'] = np.median(intensities)
        # os.chdir(outputdir)
        # writer = csv.writer(open("idadict-{}.csv".format(self.name), "w"))
        # writer.writerow(['median cc', idadict['median cc']])
        # writer.writerow(['intensities Sa/g', intensities])
        # # for key, value in idadict.items():
        # #     writer.writerow([key, value])
        # self.idadict = idadict
        # plt.title(self.name + ' IDA curves for {} record(s), Median CC = {:.2f} g'.format(len(records), np.median(intensities)))
        # plt.savefig(self.name + '-IDA-records(1-{})'.format(len(records)) + '.png')
        # plt.clf()
        plt.savefig(self.filename + '_ida_records_{}.pdf'.format(len(output_names)))
        return ida_dict

def generic_plot(df_x, df_y, name_x = 'drift', name_y = 'moment', save_folder = os.getcwd(), save_name = 'generic_plot.pdf'):
    fig = plt.figure()
    plt.plot(df_x, df_y)
    plt.xlabel(name_x)
    plt.ylabel(name_y)
    os.chdir(save_folder)
    plt.savefig(save_name)
    plt.close(fig)
    os.chdir(os.getcwd())
    return


# =====================
# kaushik = Response(filename = "kaushik_manual", input_folder = dropbox + '/str/mod/frame/_results')
# kaushik.run_opensees(tcl_folder_location = "../frames/reverse_eng")
# save_as = 'kaushik2004_cc_infill'
# df = kaushik.capacity_curve_plot('kaushik_drift.out', 'shear_envelope_Linear.out', 'kaushik2017_capacity_curves', normalize_weight = 1640.0, PDelta = 'shear_envelope_PDelta.out', collapse_drift = 0.027, period = 1.063)
# cc = kaushik.collapse_capacity(save_as, dropbox + "/str/spectra/soft_soil_collapse_spectra.csv")
#------
#------
#------ 2019.
num_ida_records = 10

kaushik_tcl = Response('kaushik_tcl', dropbox + '/str/frame/results_kaushik/')
# pushover = kaushik_tcl.capacity_curve_plot('pushover/roof_drift.out', 'pushover/base_shear_springs.out', 'pushover/kaushik_tcl_cap_curves', 1640.0, 1.0, 0.03, 0.271, 2, 1.0e-3, PDelta='pushover/base_shear.out')
# ida_dir = 'ida/record_1/intensity_0.351/'
ida_dir = 'ida/record_98/intensity_0.157/'
ida_name = 'sct'

th1 = kaushik_tcl.time_history(ida_dir + 'drifts.out', ['time', 'drift1', 'drift2', 'drift3', 'drift4'], ida_dir,  'drift_kaushik_{}.pdf'.format(ida_name))
th2 = kaushik_tcl.time_history(ida_dir + 'columns_1.out', ['time', 'col1', 'col2', 'col3', 'col4'], ida_dir, 'columns_1_{}.pdf'.format(ida_name))
th3 = kaushik_tcl.time_history(ida_dir + 'roof_drift.out', ['time', 'roof drift'], ida_dir,  'roof_drift_kaushik_{}.pdf'.format(ida_name))

generic_plot(th1['drift1'], th2['col1'], 'drift', 'moment', dropbox+'/str/frame/results_kaushik/' + ida_dir, 'displacement_col1_{}.pdf'.format(ida_name))

# mod = kaushik_tcl.mode_shapes_plot('eigen/mode_shapes.out', 'eigen/modes_plot')
# cc = kaushik_tcl.collapse_capacity('kaushik_tcl_cc', dropbox+ '/str/spectra/soft_soil_collapse_spectra.csv')
# ida = kaushik_tcl.idaplot('ida/', ['ida_performance_rec_' + str(rec+1)+'.out' for rec in range(2,num_ida_records)])
ida = kaushik_tcl.idaplot('ida/', ['ida_performance_rec_98.out'], 0.1)

# saul = Response('saul', dropbox + '/str/frame/original_frame/')
# pushover2 = saul.capacity_curve_plot('_output/SPO_output/Disp-Roof1.out', '_output/SPO_output/Vbase1.out', 'saul_cap_curves', -1640, 0.5, 0.271, 2, PDelta = '_output/SPO_output/Vbase1_PD.out')
# ida2 = saul.idaplot('_output/IDA_output/', ['ida_performance_rec_' + str(rec+1)+'.out' for rec in range(2,num_ida_records)])

#       drift  intensity
# 0  0.002405        0.5
# 1  0.005775        0.7
# 2  0.023590        0.9
# 3  0.067010        1.1
# 4  0.121300        1.3
# 5  0.166800        1.5

#       drift  intensity
# 0  0.000663      0.151
# 1  0.001761      0.250
# 2  0.003299      0.350
# 3  0.005275      0.450
# 4  0.009928      0.550
# 5  0.026370      0.650
# 6  0.060620      0.750
# 7  0.094230      0.850
# 8  0.127400      0.950
# 9  0.157600      1.050

# ran it again at 0.5 with dSa = 0.2
# set spectra_unit_factor 981;# how many `spectra_units` are in a `g`. finds intensity of an oscillator under unscaled record. used to scale.
# set record_scale_factor 1;# 9.81 if record is in g's. else 1.
# 1  0.005271        0.7
# 2  0.017840        0.9
# 3  0.062300        1.1
# 4  0.113100        1.3
# 5  0.157600        1.5
# 6  0.231500        1.7
# >>>
# drift  intensity
# 0   0.000439        0.1
# 1   0.001318        0.2
# 2   0.002635        0.3
# 3   0.004392        0.4
# 4   0.007108        0.5
# 5   0.017840        0.6
# 6   0.043890        0.7
# 7   0.074050        0.8
# 8   0.113100        0.9
# 9   0.144100        1.0
# 10  0.176900        1.1
# 11  0.239600        1.2

# 0.07 1.2 <-
# backbone resorte del muro
# grafica del ida corregida
