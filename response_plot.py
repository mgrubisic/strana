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
print(cwd)

"""TODO: include units in plot automatically? better to plot without units! normalized.
- [ ] for record
- [ ]
"""

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


    def moment_curvature_plot(self, curvature_filename, moments_filename, savename = 'moment-rotation.pdf', sign=-1):
        fig = plt.figure()
        deformation = read_csv(curvature_filename, names = ['strain', 'curvature'], sep = ' ')

        momentframe = read_csv(moments_filename, names = ['axial load', 'moment'], sep = ' ')
        dataframe = merge(deformation,momentframe, left_index=True, right_index=True)
        dataframe['moment'] = dataframe['moment']*sign
        print( dataframe)

        dataframe.plot(x='curvature', y='moment')
        os.chdir(self.output_folder)
        plt.savefig(savename)
        plt.close(fig)
        os.chdir('../')
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



    def capacity_curve_plot(self, drift_filename, column_forces_filename, savename='drifts', normalize_weight = 1.0, collapse_drift = 0.02, period = 1.0, elastic_slope_sensitivity = 2, drift_tol=1.0e-4, PDelta=False, delim_whitespace = True):
        """plots both capacity curves (Pdelta and linear) and approximates it with a
        piecewise bilinear function
        """

        os.chdir(self.input_folder)
        fig = plt.figure()
        # dataframe = read_csv(drift_filename, names = ['drift'])
        df = read_table(column_forces_filename, header=None, delim_whitespace = delim_whitespace)
        df['Regular'] = - 1.0/normalize_weight * df.sum(1)
        df['drift'] = read_table(drift_filename, names = ['drift'])

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
            print("Specify a less strict drift tolerance (1.0e-3 recommended)", err)

        fu = shear.iloc[du_index]
        fuPD = shearPD.iloc[du_index]
        Ke = shear.iloc[0 + elastic_slope_sensitivity]/drift.iloc[0 + elastic_slope_sensitivity]
        KePD = shearPD.iloc[0 + elastic_slope_sensitivity]/drift.iloc[0 + elastic_slope_sensitivity]
        # print(Ke, KePD, fu, fuPD)
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

        plt.savefig(savename+"_bilin_"+".pdf")
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

# =====================
# kaushik = Response(filename = "kaushik_manual", input_folder = dropbox + '/str/mod/frame/_results')
# kaushik.run_opensees(tcl_folder_location = "../frames/reverse_eng")
# save_as = 'kaushik2004_cc_infill'
# df = kaushik.capacity_curve_plot('kaushik_drift.out', 'shear_envelope_Linear.out', 'kaushik2017_capacity_curves', normalize_weight = 1640.0, PDelta = 'shear_envelope_PDelta.out', collapse_drift = 0.027, period = 1.063)
# cc = kaushik.collapse_capacity(save_as, dropbox + "/str/spectra/soft_soil_collapse_spectra.csv")

kaushik_tcl = Response('kaushik_tcl', dropbox + '/str/frame/kaushik_results/')
save_as = 'kaushik_tcl_cc'
df = kaushik_tcl.capacity_curve_plot('pushover/roof_drift.out', 'pushover/base_shear_springs.out', 'pushover/kaushik_tcl_cap_curves', 1640.0, 0.03, 0.333, 2, 1.0e-3, PDelta='pushover/base_shear.out')
# mod = kaushik_tcl.mode_shapes_plot('eigen/mode_shapes.out', 'eigen/modes_plot')
cc = kaushik_tcl.collapse_capacity('kaushik_tcl_cc', dropbox+ '/str/spectra/soft_soil_collapse_spectra.csv')
