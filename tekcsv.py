""" Read a tektronix folder and transform a Y-time oscilogram
    to a X-Y plot, using Channels 1 and 2.
"""

import math
import pandas as pd
import pylab
from scipy.optimize import curve_fit

def linear(x, a, b): # pylint: disable=invalid-name
    """ linear function equation for curve fitting
    """
    return a*x + b

class Channel():
    """ Extracts info from a channel CSV file
    """

    def get_param(self, param_name):
        """ Get a string parameter from the channel csv file
        """
        return self.csv_data[self.csv_data['parameters'] == param_name]['param_values'].values[0]

    def get_float_param(self, param_name):
        """ Get a float parameter from the channel csv file
        """
        return float(self.get_param(param_name))

    def get_int_param(self, param_name):
        """ Get an integer parameter from the channel csv file
        """
        return round(self.get_float_param(param_name))

    def __init__(self, csv_file):

        # read data to dataFrame
        self.csv_data = pd.read_csv(
            csv_file,
            sep=',',
            names=['parameters', 'param_values', 'empty', 'time', 'values', 'empty2']
        )

        # storing x and y values
        self.time = self.csv_data['time']
        self.value = self.csv_data['values']

        # storing aquisition parameters
        self.parameters = {}
        self.parameters['record_lenght'] = self.get_int_param('Record Length')
        self.parameters['sample_interval'] = self.get_float_param('Sample Interval')
        self.parameters['trigger_point'] = self.get_int_param('Trigger Point')

        self.parameters['source'] = self.get_param('Source')

        self.parameters['vertical_units'] = self.get_param('Vertical Units')
        self.parameters['vertical_scale'] = self.get_float_param('Vertical Scale')
        self.parameters['vertical_offset'] = self.get_float_param('Vertical Offset')

        self.parameters['horizontal_units'] = self.get_param('Horizontal Units')
        self.parameters['horizontal_scale'] = self.get_float_param('Horizontal Scale')

        self.parameters['pt_fmt'] = self.get_param('Pt Fmt')
        self.parameters['yzero'] = self.get_float_param('Yzero')

        self.parameters['probe_atten'] = self.get_float_param('Probe Atten')
        self.parameters['firmware_version'] = self.get_param('Firmware Version')

        # getting the rms value as the normalization value
        #self.norm = math.sqrt(sum(x**2 for x in self.value)/self.parameters['record_lenght'])
        self.norm = 1
        
        # saves the normalized y values
        self.value_normalized = self.value / self.norm

    def invert(self):
        """ Invert channel input (flips over y=0 axis).
        """
        self.value = -self.value
        self.value_normalized = -self.value_normalized

def analyze_aquisition(basefile):
    """ Giving basefile with ch1 data, plot Ch1 x Ch2
        as XY normalized plot, with linear fitting curve.
    """

    ch1 = Channel(basefile)
    ch2 = Channel(basefile.replace('CH1', 'CH2'))
    ch3 = Channel(basefile.replace('CH1', 'CH3'))

    #offset = syncronize(ch1.value_normalized,ch2.value_normalized)
    
    best_cov = 1e100
    for offset in range(-20,0):
        if offset>0:
            ch1.value_normalized2 = ch1.value_normalized.tolist()[:-offset]
            ch2.value_normalized2 = ch2.value_normalized.tolist()[offset:]
        elif offset<0:
            ch1.value_normalized2 = ch1.value_normalized.tolist()[-offset:]
            ch2.value_normalized2 = ch2.value_normalized.tolist()[:offset]

        pars, _cov = curve_fit(
            f=linear,
            xdata=ch1.value_normalized2,
            ydata=ch2.value_normalized2,
            p0=[1, 0],
            bounds=(-1e100, 1e100)
        )
        cov = _cov[0][0]**2 + _cov[0][1]**2 + _cov[1][0]**2 + _cov[1][1]**2
        if cov < best_cov:
            best_cov = cov
            best_pars = pars
            best_offset = offset
            
    
    ch1.value_normalized2 = (ch1.value_normalized/8.06).tolist()[-best_offset:]
    ch2.value_normalized2 = (ch2.value_normalized*(1.8+3.3)/3.3+2.5).tolist()[:best_offset]
    
    pars, _cov = curve_fit(
        f=linear,
        xdata=ch1.value_normalized2,
        ydata=ch2.value_normalized2,
        p0=[1, 0],
        bounds=(-1e100, 1e100)
    )
    lin_x = [-2.5/0.185, 2.5/0.185]
    lin_y = [x*0.185+2.5 for x in lin_x]
    #fit_y = [linear(x, pars[0], pars[1]) for x in lin_x]    
    

    pylab.figure()

    pylab.title("Saturação do sensor de corrente ACS712-05B (curva XY)")
    pylab.plot(ch1.value_normalized2, ch2.value_normalized2, label='Ganho do sensor')
    pylab.plot(lin_x, lin_y, label='ganho ideal')
    #pylab.plot(lin_x, fit_y, label='fitted {:.3f}*x{:+.3f}'.format(pars[0], pars[1]))
    pylab.xlabel("Corrente [A] Sobre a carga Resistiva de 8.06$\Omega$")
    pylab.ylabel("Tensão [V] de saída do sensor")
    
    pylab.yticks([ 0.37, 2.5,  4.76])
    pylab.xticks([-11.7, 0, 12.20])
    pylab.grid()

    pylab.legend()
    
    
    pylab.figure()
    pylab.title("Saturação do sensor de corrente ACS712-05B (curva YT)")
    pylab.plot(ch2.value_normalized/0.185*(3.3+1.8)/3.3, label = "Saída do sensor normalizada para a corrente de entrada")
    pylab.plot(ch1.value_normalized/8.06, label="Corrente [A] Sobre a carga Resistiva de 8.06$\Omega$")
    pylab.legend()
    
    pylab.yticks([-15,-11.7,0,12.2,15])
    pylab.xticks([])
    pylab.grid(axis='y')

def syncronize(xdata,ydata,testrange = 25):
    """ Find the index offset between datas so that we maximize their convolution
        (that is, they will be the most likely to each other)
    """
    def convolution(x,y):
        length = len(x)
        return sum(x[i] * y[i] for i in range(length))/length
    
        
    xdata = xdata.tolist()
    ydata = ydata.tolist()
    
    max_conv = convolution(xdata,ydata)
    max_i = 0
    print(f'i = {0}, conv = {max_conv}')
    
    for i in range(1,testrange):      
        print(i)
        # testing shorting ydata
        x = xdata[:-i]
        y = ydata[i:]
        conv = convolution(x,y)
        
        print(f'i = {i}, conv = {conv}')
        if conv > max_conv:
            max_conv = conv
            max_i = i
            
        # testing shorting xdata   
        x = xdata[i:]         
        y = ydata[:-i]
        conv = convolution(x,y)
        print(f'i = {-i}, conv = {conv}')
        if conv > max_conv:
            max_conv = conv
            max_i = -i
    
    print(f'max_i = {max_i}, max_conv = max_conv')
    return max_i       
    
        
    

for i in range(40,40+1):
    analyze_aquisition(f"DATA\\ALL00{i}\\F00{i}CH1.CSV")

pylab.show()
