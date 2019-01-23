""" input ground motion and strata, output disps, vels, accs at surface
TODO: output solution fields at interfaces
TODO, make spectra out of multiple records
def spectra(record_array, frequencies):
"""
from pandas import read_csv
from numpy import exp, cos, sin, array, linspace, pi, fft, genfromtxt, cos, concatenate, append, convolve, clip
import seaborn as sns
import easyplot as ep
import matplotlib.pyplot as plt


class Record(object):
    def __init__(self, filename, steps, dt=0.01, kind='acc'):
        self.th = genfromtxt(filename)
        self.dt = dt
        self.steps = steps
        self.fft_steps = steps/2+1 if steps%2==0 else (steps/2)+1
        self.kind = kind
        self.duration = steps*dt
        self.timeseries = linspace(0, self.duration, steps)
        self.name = filename


class Soil(object):
    """initializes an elastic silty clay soil object"""
    def __init__(self, weight = 12.0, gravity = 9.81): 
        self.weight = weight
        self.density = weight/gravity
    def properties(self,thickness=1.0, shear_modulus=5000.0, damping=0.05):
        self.thickness = thickness
        self.shear_stiffness = shear_modulus
        self.x1 = thickness
        self.damping = damping
        self.complex_shear = shear_modulus*(1+2*damping**2+2j*damping*(1+damping**2)**.5)
        self.velocity = (self.complex_shear/self.density)**.5

def amplitudes(amplitude_in, amplitude_out, wavenumber, thickness, impedance):
    incident = 0.5*(amplitude_in*exp(1j*wavenumber*thickness)*(1+impedance) + amplitude_out*exp(-1j*wavenumber*thickness)*(1-impedance))
    reflected = 0.5*(amplitude_in*exp(1j*wavenumber*thickness)*(1-impedance) + amplitude_out*exp(-1j*wavenumber*thickness)*(1+impedance))
    return incident, reflected

def transfer_function(layers, frequencies_hz, name):
    """ returns the transfer function for unitary base amplitude for given layer configuration, top to bottom"""
    transference = []
    for frequency in frequencies_hz:
        incident = [1.0]
        reflected = [1.0]
        for m, layer in enumerate(layers):
            if m+1 == len(layers):
                continue
            impedance = layer.density*layer.velocity/layers[m+1].density/layers[m+1].velocity
            layer.wavenumber = 2*pi*frequency/layer.velocity
            incident.append(amplitudes(incident[m], reflected[m], layer.wavenumber, layer.thickness, impedance)[0])
            reflected.append(amplitudes(incident[m], reflected[m], layer.wavenumber, layer.thickness,impedance)[1])

        transference.append(2.0/abs(incident[-1].real + reflected[-1].real))
    ep.easyplot(frequencies_hz, transference, 'Hz', 'Transfer amplitude', name, save_as=name)
    return transference


def site_response(record, transference, clip_ampl = 50.0, max_freq=25.0, max_ampl = 50.0):
    """transfer function should be frequency compatible with record
    in the sense of SPACING!, length will be filled with zeroes
    TODO: strata name configuration
    """
    fourier = fft.rfft(record.th)
    ampl = array([(a.real**2+a.imag**2)**.5 for a in fourier])
    frequencies = linspace(0, 1.0/2/record.dt, len(ampl))
    transference = clip(transference, 0, clip_ampl)
    frequency_response =  fourier * transference
    time_response = fft.irfft(frequency_response)

    figure, (record_th, record_fft, transfer_plot, frequency_plot, time_plot) = plt.subplots(5, 1, figsize=(10, 15))
    figure.subplots_adjust(hspace = 0.4)
    figure.suptitle('Site response for record ' + record.name + ', strata configuration TODO')
    record_th.plot(record.timeseries, record.th, lw=0.5 , label='Record time history')
    record_th.set_xlabel('[s]')
    # record_th.xaxis.set_label_coords(1.05, -0.025)
    record_th.set_ylabel('Acceleration')
    record_th.legend()

    record_fft.plot(frequencies, ampl, lw=0.5, label='Record Fourier amplitude')
    record_fft.set_xlabel('[Hz]')
    record_fft.set_ylabel('Amplitude')
    record_fft.set_xlim([0, max_freq])
    record_fft.legend()

    record_fft.plot(frequencies, ampl, lw=0.5, label='Record Fourier amplitude')

    transfer_plot.plot(frequencies, transference, lw=0.5, label='Transfer function for strata')
    transfer_plot.set_xlabel('Hz')
    transfer_plot.set_ylabel('Amplitude')
    transfer_plot.set_xlim([0, max_freq])
    transfer_plot.set_ylim([0, max_ampl])
    transfer_plot.legend()

    frequency_plot.plot(frequencies, frequency_response, lw=0.5, label='Frequency response at surface')
    frequency_plot.set_xlabel('Hz')
    frequency_plot.set_ylabel('Amplitude')
    frequency_plot.set_xlim([0, max_freq])
    # frequency_plot.set_ylim([0, 10*clip_ampl])
    frequency_plot.legend()

    time_plot.plot(record.timeseries[1:], time_response, lw=0.5, label='Time response at surface')
    time_plot.set_xlabel('[s]')
    time_plot.set_ylabel('Accel')
    time_plot.legend()

    plt.savefig('site.response.record.7.pdf')
    plt.clf()
    return fourier, frequency_response
mu = 1.0*array([1.0, 4.0, 71.111])

w = 9.81*array([mu[0]/350**2, mu[1]/600**2, mu[2]/1500**2])

s1 = Soil(w[0])
s1.properties(5.0, mu[0], 0.05)

s2 = Soil(w[1])
s2.properties(5.0, mu[1], 0.0033)

s3 = Soil(w[2])
s3.properties(20.0, mu[2], 0.0013)
# # dt = 0.005
# damping = 0.15
# # steps = 7991
# nr = Record('recordset/01_NR94mu1.th', 2999, 0.01)
# freqs = linspace(0, 25.0, nr.fft_steps)
# # lp = Record('recordset/26_LP89ca2.th', steps = steps, dt = dt)
# s1 = Soil(18.9)
# s1.properties(50.0, shear_modulus = 50**2*18.9/9.81, damping = 0.15)
# # s2 = Soil(17.9)
# # s2.properties(38.0, shear_modulus = 300**2*18.9/9.81, damping = 0.03)
# roca = Soil(2.24*9.81)
# roca.properties(thickness=1e3, shear_modulus=1500.0**2*2.24, damping = 0.01)

# tf = transfer_function([s1, roca], freqs, 'elastic.bedrock4')
# fa = site_response(nr, tf)

cna = Record('recordset/eputf125.002', 4999, 0.02)
freqs = linspace(0, 15.0, cna.fft_steps)
tf = transfer_function([s1, s2,s3], freqs, 'ejemplo.cna2')
fou, fr = site_response(cna, tf)
