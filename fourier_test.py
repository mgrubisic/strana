import easyplot as ep
from numpy import pi, linspace, cos, sin, fft


freq1, freq2 = 5, 25
ampl = 2.0

bins = 200
dt = 0.01
duration = bins*dt


timeseries = linspace(0, duration, bins)
frequency = linspace(0, 1./dt/2, bins/2)
signal = ampl*sin(2*pi*freq1*timeseries) + 5*sin(2*pi*freq2*timeseries)

amplitudes = fft.rfft(signal, bins-1)

ep.easyplot(frequency, 2*abs(amplitudes)/bins, 'Hz', 'Amplitudes',  save_as='recordset/sinewave', xlim=[0, 30])
