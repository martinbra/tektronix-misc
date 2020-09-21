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
        self.norm = math.sqrt(sum(x**2 for x in self.value)/self.parameters['record_lenght'])

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

    pars, _cov = curve_fit(
        f=linear,
        xdata=ch1.value_normalized,
        ydata=ch2.value_normalized,
        p0=[1, 0],
        bounds=(-1e100, 1e100)
    )

    lin_x = [-1.5, 1.5]
    lin_y = [linear(x, 1, 0) for x in lin_x]
    fit_y = [linear(x, pars[0], pars[1]) for x in lin_x]

    pylab.figure()

    pylab.title(f"{basefile} {str(max(ch1.value))} Vmax")
    pylab.plot(ch1.value_normalized, ch2.value_normalized)
    pylab.plot(lin_x, lin_y, label='ideal')
    pylab.plot(lin_x, fit_y, label='fitted {:.3f}*x{:+.3f}'.format(pars[0], pars[1]))

    pylab.legend()

for i in range(29, 38+1):
    analyze_aquisition(f"DATA\\ALL00{i}\\F00{i}CH1.CSV")

pylab.show()
