import pylab
import numpy as np

acs_gain = 0.185 #V/A
acs_suply = 5 #V
acs_offset = acs_suply/2 #V

offset_sat = 0.3

current_min = -acs_offset/acs_gain
current_max = (acs_suply-acs_offset)/acs_gain

current_min_sat = (offset_sat-acs_offset)/acs_gain
current_max_sat = (acs_suply-offset_sat-acs_offset)/acs_gain

current_rms_max = min(current_max_sat,-current_min_sat)/2**0.5

print(current_min_sat, current_max_sat, current_rms_max)

acs_current = np.array(np.arange(current_min,current_max,0.1))
acs_output_teoretical = acs_current * acs_gain + acs_offset
acs_output_saturated = np.array([min(max(0.4,x),4.6) for x in acs_output_teoretical])

pylab.plot(acs_current,acs_output_teoretical)
pylab.plot(acs_current,acs_output_saturated)
pylab.show()

