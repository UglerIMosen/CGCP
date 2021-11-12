# coding: utf-8
import datetime
if datetime.datetime.today().weekday() == 4:
    #this is a really stupid animation. Happy friday
    from additional_stuff.GcFunc import GcFridaySpecial
    GcFridaySpecial()
else:
    from additional_stuff.GcFunc import GcLogo
    GcLogo()

import sys
import os
from pprint import pprint  # pprint == pretty print

# Python version check
if sys.version_info.major < 3: #courtesy of Keni 
    print("No! Use Python 3")
    sys.exit(1)

"""
#MANUAL SETTINGS
""" 
from additional_stuff.GcFunc import get_path
path_tool = get_path()
path = path_tool.ask_path()

if 'manual_settings.py' in os.listdir(path):
    sys.path.append(path)
    from manual_settings import MANUAL_SETTINGS
    print(' ')
    print('Using fitting information from data path') 
else:
    from additional_stuff.manual_settings import MANUAL_SETTINGS
    print(' ')
    print('Using default fitting information')
    
info = MANUAL_SETTINGS()

from additional_stuff.GcFunc import data_loading_tools
load_tools = data_loading_tools(info)

from additional_stuff.GcFunc import GC_plots
gc_plots = GC_plots(info)

print(' ')
print('All settings are configured:')
print('* Path of data: '+path)
print(' ')

"""
#Load Functions and settings
"""
import numpy as np
from matplotlib import pyplot as plt
import pickle
import os

#plot-design
import matplotlib as mpl
font = {'family' : 'serif',#'Arial',
        #'weight' : 'bold',
        'size'   : 18}
mpl.rc('font', **font)

#Figure size
from pylab import rcParams
rcParams['figure.figsize'] = 10, 7

print(' ')
print('All functions are loaded succesfully.')
print(' ')

folder = load_tools.make_numbered_folder(path+'/'+str(datetime.date.today()))
print('Data is saved to "'+folder+'".')
print(' ')

from additional_stuff.settings_generator import settings_generator
func = settings_generator(folder,info)
func.write_file()

"""
#Load and plot raw GC data
"""
# Load Sequence. 
sequence = load_tools.Sequence(path)
gc_plots.raw_data_heat_plot(sequence) #fig 1 and 2
print('Raw data has been plotted succesfully')
print(' ')

"""
#LOAD SETUP DATA
"""
##import reactor temperature data recorded in LabView
rel_time, reactor_temperature, reactor_pressure, massflow_O2, massflow_H2, massflow_CO2, massflow_CO, massflow_Ar, massflow_total = load_tools.load_setup_data(path)
gc_plots.MFC_plot(rel_time, reactor_temperature, reactor_pressure, massflow_O2, massflow_H2, massflow_CO2, massflow_CO, massflow_Ar) #fig 3

"""
#Data treatment
"""
from additional_stuff.GcFunc import integration_tools
int_tools = integration_tools(info)

print('----------------------------------------------------')
print(' ')
print('Start integrating of GC spectra. Progress:')

injectiontimes_sec = []
injection_pressure = []
for injection in sequence:
    cache = injection['injection time']
    injectiontimes_sec.append(cache)
    cache = injection['line pressure']
    injection_pressure.append(cache)
injectiontimes_sec = np.array(injectiontimes_sec)
injection_pressure = np.array(injection_pressure)

FID_res = {}
FID_err = {}
for peak in gc_plots.info.fit_info['FID']:
    FID_res[peak] = []
    FID_err[peak] = []
for j, injection in enumerate(sequence):
    print(f"{50*(j-gc_plots.info.skip_start)/len(injectiontimes_sec):.1f} %\r",end="")
    spectrum_area = int_tools.integrate_spectrum(injection,'FID')
    for peak in gc_plots.info.fit_info['FID']:
        FID_res[peak].append(spectrum_area[0][peak])
        FID_err[peak].append(spectrum_area[1][peak])
        
TCD_res = {}
TCD_err = {}
for peak in gc_plots.info.fit_info['TCD']:
    TCD_res[peak]=[]
    TCD_err[peak]=[]
for j, injection in enumerate(sequence):
    print(f"{(50*(j-gc_plots.info.skip_start)/len(injectiontimes_sec))+50:.1f} %\r",end="")
    spectrum_area = int_tools.integrate_spectrum(injection, 'TCD')
    for peak in gc_plots.info.fit_info['TCD']:
        TCD_res[peak].append(spectrum_area[0][peak])
        TCD_err[peak].append(spectrum_area[1][peak])

gc_plots.integrated_vs_time(injectiontimes_sec,FID_res,FID_err,'FID','GC FID',rel_time,reactor_temperature,reactor_pressure) #fig 4
gc_plots.integrated_vs_time(injectiontimes_sec,TCD_res,TCD_err,'TCD','GC TCD',rel_time,reactor_temperature,reactor_pressure) #fig 5

#Write data to .txt file
int_tools.write_to_txt(folder+'/FID_raw_data', FID_res, dimension='time', data_x=injectiontimes_sec)
int_tools.write_to_txt(folder+'/TCD_raw_data', TCD_res, dimension='time', data_x=injectiontimes_sec)
int_tools.write_to_txt(folder+'/FID_raw_data_error', FID_err, dimension='time', data_x=injectiontimes_sec)
int_tools.write_to_txt(folder+'/TCD_raw_data_error', TCD_err, dimension='time', data_x=injectiontimes_sec)

print('Integration of GC spectra is done')
print(' ')

"""
#Averaging Labview-data for synchronization with GC
"""
rel_time_GCmeaned = []
temperature_GCmeaned = []
pressure_GCmeaned = []
massflow_H2_GCmeaned = []
massflow_CO2_GCmeaned = []
massflow_CO_GCmeaned = []
massflow_Ar_GCmeaned = []
massflow_O2_GCmeaned = []
massflow_total_GCmeaned = []

cache_time = []
cache_pressure = []
cache_temperature = []
cache_H2 = []
cache_CO2 = []
cache_CO = []
cache_Ar = []
cache_O2 = []
cache_total = []

snif_start = 200 #how long before the injection
snif_end = 50 #how long before the injection

n=0
for time,pressure,temperature,H2_flow,CO2_flow,CO_flow,Ar_flow,O2_flow,total_flow in zip(rel_time,
                                                                                reactor_pressure,
                                                                                reactor_temperature,
                                                                                massflow_H2,
                                                                                massflow_CO2,
                                                                                massflow_CO,
                                                                                massflow_Ar,
                                                                                massflow_O2,
                                                                                massflow_total):
    if (injectiontimes_sec[n]-snif_start <= time <= injectiontimes_sec[n]-snif_end): #Snif 
        cache_time.append(time)
        cache_pressure.append(pressure)
        cache_temperature.append(temperature)
        cache_H2.append(H2_flow)
        cache_CO2.append(CO2_flow)
        cache_CO.append(CO_flow)
        cache_Ar.append(Ar_flow)
        cache_O2.append(O2_flow)
        cache_total.append(total_flow)
    elif time > snif_end+injectiontimes_sec[n]:#injectiontimes_sec[0]-snif_start+injectiontimes_sec[n]:
        rel_time_GCmeaned.append(injectiontimes_sec[n])
        pressure_GCmeaned.append(np.mean(cache_pressure))
        temperature_GCmeaned.append(np.mean(cache_temperature))
        massflow_H2_GCmeaned.append(np.mean(cache_H2))
        massflow_CO2_GCmeaned.append(np.mean(cache_CO2))
        massflow_CO_GCmeaned.append(np.mean(cache_CO))
        massflow_Ar_GCmeaned.append(np.mean(cache_Ar))
        massflow_O2_GCmeaned.append(np.mean(cache_O2))
        massflow_total_GCmeaned.append(np.mean(cache_total))
        cache_time = []
        cache_pressure = []
        cache_temperature = []
        cache_H2 = []
        cache_CO2 = []
        cache_CO = []
        cache_Ar = []
        cache_O2 = []
        cache_total = []
        n = n+1
    if n == len(injectiontimes_sec):
        break

rel_time_GCmeaned = np.array(rel_time_GCmeaned)
temperature_GCmeaned = np.array(temperature_GCmeaned)
pressure_GCmeaned = np.array(pressure_GCmeaned)
massflow_H2_GCmeaned = np.array(massflow_H2_GCmeaned)
massflow_CO2_GCmeaned = np.array(massflow_CO2_GCmeaned)
massflow_CO_GCmeaned = np.array(massflow_CO_GCmeaned)
massflow_Ar_GCmeaned = np.array(massflow_Ar_GCmeaned)
massflow_O2_GCmeaned = np.array(massflow_O2_GCmeaned)
massflow_total_GCmeaned = np.array(massflow_total_GCmeaned)

"""
#Concentrations
"""
#Ar, H2, CO2, CO and CH4 have been calibrated in March 2018
#MeOH and DME have been calibrated by Irek in 2013. Conversion from his calibration to the present one has been made.
#
#Temperature is taken as constant and thus is automatically included in the conversion-constants for detector to concentration. 
#Pressure at the GC-injection must be known. It is assumed to be the same as the inlet pressure.

GC_conversionP_to_Perc = gc_plots.info.GC_conversion_to_Perc

concentrations = {}
concentrations['time'] = rel_time_GCmeaned
concentrations['temperature'] = temperature_GCmeaned
concentrations['pressure'] = pressure_GCmeaned
concentrations['injection_pressure'] = injection_pressure
concentrations['total_massflow'] = massflow_total_GCmeaned
concentrations['MFC1'] = massflow_O2_GCmeaned
concentrations['MFC2'] = massflow_H2_GCmeaned
concentrations['MFC3'] = massflow_CO2_GCmeaned
concentrations['MFC4'] = massflow_CO_GCmeaned
concentrations['MFC5'] = massflow_Ar_GCmeaned

if len(injectiontimes_sec) != len(concentrations['time']):
    print('Warning: Uneven amount of GC- and MFC-data')
end_index = len(concentrations['time'])

"""
Fixing injection pressure
"""
if injection_pressure[0] != None:
    broken = np.where(injection_pressure < 0)[0]
    if len(broken) != 0:
        print('------------------------------------------')
        print('Issues with injection pressure encountered:')
        if len(broken) == len(injection_pressure):
            print('No injection pressure to use. Setting it to 1 bar')
            injection_pressure = np.ones(len(injection_pressure))
        else:
            print('Some injection pressures are faulty. Basing injection pressure on what values are available.')
            print('The following injections had a faulty injection pressure:')
            print(broken)
            for index in [*broken,*np.flip(broken)]:
                if injection_pressure[index] < 0 and index != 0 and index != len(injection_pressure)-1:
                    if injection_pressure[index-1] > 0 and injection_pressure[index+1] > 0:
                        injection_pressure[index] = (injection_pressure[index-1]+injection_pressure[index+1])/2
                    elif injection_pressure[index-1] > 0:
                        injection_pressure[index] = injection_pressure[index-1]
                    elif injection_pressure[index+1] > 0:
                        injection_pressure[index] = injection_pressure[index+1]
                elif injection_pressure[index] < 0 and index == 0:
                    injection_pressure[index] = injection_pressure[index+1]
                elif injection_pressure[index] < 0 and index == len(injection_pressure)-1:
                    injection_pressure[index] = injection_pressure[index-1]
            
    
for gas in gc_plots.info.fit_info['FID']:
    try:
        if injection_pressure[0] == None:
            concentrations[gas] = GC_conversionP_to_Perc['FID'][gas]*np.array(FID_res[gas][:end_index])/pressure_GCmeaned
        else:
            concentrations[gas] = GC_conversionP_to_Perc['FID'][gas]*np.array(FID_res[gas][:end_index])/injection_pressure
    except:
        print('Following gas is not calibrated: '+gas)
for gas in gc_plots.info.fit_info['TCD']:
    try:
        if injection_pressure[0] == None:
            concentrations[gas] = GC_conversionP_to_Perc['TCD'][gas]*np.array(TCD_res[gas][:end_index])/pressure_GCmeaned
        else:
            concentrations[gas] = GC_conversionP_to_Perc['TCD'][gas]*np.array(TCD_res[gas][:end_index])/injection_pressure
    except:
        print('Following gas is not calibrated: '+gas)

#Write data to .txt file
int_tools.write_to_txt(folder+'/concentration_vs_time',concentrations)
pickle.dump(concentrations, open(folder+'/concentration_vs_time'+".p", "wb"))

"""
#Concentration Figures
"""
gas_list = []
for gas in [*gc_plots.info.fit_info['TCD'].keys(),*gc_plots.info.fit_info['FID'].keys()]:
    if gas in list(GC_conversionP_to_Perc['TCD'].keys())+list(GC_conversionP_to_Perc['FID'].keys()):
        gas_list.append(gas)

gc_plots.log_concentration_vs_time(concentrations['time'],concentrations,gas_list,rel_time,reactor_temperature) #fig 6
gc_plots.log_concentration_vs_time(concentrations['time'],concentrations,gc_plots.info.linear_plotted_gasses,rel_time,reactor_temperature) #fig 7

"""
#Temperature-dependant data
"""
concentration_vs_temp = {}
#first all values are meaned over temperature and put into synchronized lists

for gas in gas_list:#list(GC_conversionP_to_Perc['TCD'].keys())+list(GC_conversionP_to_Perc['FID'].keys()):
    concentration_vs_temp['temperature']        = []
    concentration_vs_temp['temperature_std']    = []
    concentration_vs_temp['total_massflow']     = []
    concentration_vs_temp['total_massflow_std'] = []
    concentration_vs_temp['MFC1']               = []
    concentration_vs_temp['MFC2']               = []
    concentration_vs_temp['MFC3']               = []
    concentration_vs_temp['MFC4']               = []
    concentration_vs_temp['MFC5']               = []
    concentration_vs_temp['MFC1_std']           = []
    concentration_vs_temp['MFC2_std']           = []
    concentration_vs_temp['MFC3_std']           = []
    concentration_vs_temp['MFC4_std']           = []
    concentration_vs_temp['MFC5_std']           = []
    concentration_vs_temp['pressure']           = []
    concentration_vs_temp['pressure_std']       = []
    concentration_vs_temp[gas]                  = []
    concentration_vs_temp[gas+'_herr']          = [] #high value
    concentration_vs_temp[gas+'_lerr']          = [] #low value

    
    cache      = [concentrations[gas][0]]
    cache_1    = [temperature_GCmeaned[0]]
    cache_2    = [massflow_total_GCmeaned[0]]
    cache_3    = [pressure_GCmeaned[0]]
    cache_mfc1 = [massflow_O2_GCmeaned[0]]
    cache_mfc2 = [massflow_H2_GCmeaned[0]]
    cache_mfc3 = [massflow_CO2_GCmeaned[0]]
    cache_mfc4 = [massflow_CO_GCmeaned[0]]
    cache_mfc5 = [massflow_Ar_GCmeaned[0]]    
    for index, value in enumerate(np.diff(temperature_GCmeaned[1:len(temperature_GCmeaned)])):
        if abs(value) > 4: #when temperature changes less than 4 degrees
            if len(cache) > 1:
                cache   = cache[1:len(cache)]
            if len(cache_1) > 1: #removes first entry
                cache_1 = cache_1[1:len(cache_1)]
            if len(cache_2) > 1:
                cache_2 = cache_2[1:len(cache_2)]
            if len(cache_3) > 1:
                cache_3 = cache_3[1:len(cache_3)]
            if len(cache_mfc1) > 1:
                cache_mfc1 = cache_mfc1[1:len(cache_mfc1)]
            if len(cache_mfc2) > 1:
                cache_mfc2 = cache_mfc2[1:len(cache_mfc2)]
            if len(cache_mfc3) > 1:
                cache_mfc3 = cache_mfc3[1:len(cache_mfc3)]
            if len(cache_mfc4) > 1:
                cache_mfc4 = cache_mfc4[1:len(cache_mfc4)]
            if len(cache_mfc5) > 1:
                cache_mfc5 = cache_mfc5[1:len(cache_mfc5)]
            if len(cache) > 0:
                concentration_vs_temp['temperature'].append(np.mean(cache_1))
                concentration_vs_temp['temperature_std'].append(np.std(cache_1))
                concentration_vs_temp['total_massflow'].append(np.mean(cache_2))
                concentration_vs_temp['total_massflow_std'].append(np.std(cache_2))
                concentration_vs_temp['pressure'].append(np.mean(cache_3))
                concentration_vs_temp['pressure_std'].append(np.std(cache_3))
                concentration_vs_temp['MFC1'].append(np.mean(cache_mfc1))
                concentration_vs_temp['MFC1_std'].append(np.std(cache_mfc1))
                concentration_vs_temp['MFC2'].append(np.mean(cache_mfc2))
                concentration_vs_temp['MFC2_std'].append(np.std(cache_mfc2))                
                concentration_vs_temp['MFC3'].append(np.mean(cache_mfc3))
                concentration_vs_temp['MFC3_std'].append(np.std(cache_mfc3))                
                concentration_vs_temp['MFC4'].append(np.mean(cache_mfc4))
                concentration_vs_temp['MFC4_std'].append(np.std(cache_mfc4))                
                concentration_vs_temp['MFC5'].append(np.mean(cache_mfc5))
                concentration_vs_temp['MFC5_std'].append(np.std(cache_mfc5))                
                concentration_vs_temp[gas].append(np.mean(cache))
                concentration_vs_temp[gas+'_herr'].append(max(cache)-np.mean(cache))
                concentration_vs_temp[gas+'_lerr'].append(np.mean(cache)-min(cache))
            cache = []
            cache_1 = []
            cache_2 = []
            cache_3 = []
            cache_mfc1 = []
            cache_mfc2 = []
            cache_mfc3 = []
            cache_mfc4 = []
            cache_mfc5 = []
        else:
            cache.append(concentrations[gas][index+1]) #first index was written before loop
            cache_1.append(temperature_GCmeaned[index+1])
            cache_2.append(massflow_total_GCmeaned[index+1])
            cache_3.append(pressure_GCmeaned[index+1])
            cache_mfc1.append(massflow_O2_GCmeaned[index+1])
            cache_mfc2.append(massflow_H2_GCmeaned[index+1])
            cache_mfc3.append(massflow_CO2_GCmeaned[index+1])
            cache_mfc4.append(massflow_CO_GCmeaned[index+1])
            cache_mfc5.append(massflow_Ar_GCmeaned[index+1])

#Write data to .txt file
int_tools.write_to_txt(folder+'/concentration_vs_temp',concentration_vs_temp)
pickle.dump(concentration_vs_temp, open(folder+'/concentration_vs_temp'+".p", "wb"))

gc_plots.lin_concentration_vs_temp(concentration_vs_temp) #fig 8

"""
Inspection Tools #This whole section can be excluded
"""
#This next section could probably be cleaned up a lot. It holds an interactive plot of all the fitted spectra together with the time-change,
#so the raw data can be examined much faster. One note of caution: When changing the spectrum to be viewed, the program actually re-fits the
#background, which means that if you have a lot of gasses, it could be tiresome to interact with. I have no idea of the limit, but it is there.

from matplotlib.widgets import Slider

diff_time_hrs = np.mean(np.diff(injectiontimes_sec))/(60*60)
###
#FID
###

def integrated_data_plot(result_dict,error_dict,detector_key,title,ylabel,placement=111):
    f = fig.add_subplot(placement)
    for entry in result_dict:
        f.plot(injectiontimes_sec/(60*60),result_dict[entry],ms=1, label=gc_plots.info.fit_info[detector_key][entry]['name'],color=gc_plots.info.fit_info[detector_key][entry]['color'])
        f.errorbar(injectiontimes_sec/(60*60),result_dict[entry],yerr=error_dict[entry],fmt='o',capsize=1, color=gc_plots.info.fit_info[detector_key][entry]['color'])
    f.set_xlabel('time [hrs]')
    f.set_ylabel(ylabel)
    f.tick_params('y')
    f.axes.set_title(title)
    g = f.twinx()
    g.plot(rel_time/(60*60),reactor_temperature,'r--')
    g.set_ylabel('Temperature [celsius]',color='r')
    g.tick_params('y',colors='r')
    f.legend(frameon=0)
    
fig = plt.figure()
integrated_data_plot(FID_res,FID_err,'FID','GC FID','Integrated FID [Wb]',placement=211)
spectrum = (sequence[0]['retention time']/60, sequence[0]['FID'])
f = fig.add_subplot(212)
GC_line, = f.plot(spectrum[0],spectrum[1],'-',markersize=1.1)

for peak in gc_plots.info.fit_info['FID']:
    fitting_results = int_tools.fit_to_spectrum(spectrum, peak, 'FID')
    if 'background' in gc_plots.info.fit_info['FID'][peak]:
        globals()['bg_plot%s' %peak], = f.plot(fitting_results[0],int_tools.background_func[gc_plots.info.fit_info['FID'][peak]['background']](fitting_results[0],fitting_results[2],*fitting_results[3]))
    else:
        globals()['bg_plot%s' %peak], = f.plot(fitting_results[0],int_tools.errorfunc(fitting_results[0],fitting_results[2],*fitting_results[3]))
f.set_ylabel('Intensity [V]')
f.axes.set_ylim([-0.1, 10])

GC_injection_slider = Slider(plt.axes([0.16, 0.035, 0.705, 0.03]), 'Time (hrs)', injectiontimes_sec[0]/(60*60), injectiontimes_sec[-1]/(60*60),valinit=injectiontimes_sec[0]/(60*60),valstep=diff_time_hrs)
Axis_slider = Slider(plt.axes([0.16, 0.0, 0.705, 0.03]), 'Axis', -0.1, 10,valinit=1,valstep=0.001)

def update_injection(val):
    cache = int(round((val)*(1/np.mean(np.diff(injectiontimes_sec/(60*60))))+1))
    if int(cache) >= len(sequence):
        cache = len(sequence)-1
    spectrum = (sequence[int(cache)]['retention time']/60, sequence[int(cache)]['FID'])
    GC_line.set_xdata(spectrum[0])
    GC_line.set_ydata(spectrum[1])
    for peak in gc_plots.info.fit_info['FID']:
        try:
            fitting_results = int_tools.fit_to_spectrum(spectrum, peak, 'FID')
            globals()['bg_plot%s' %peak].set_xdata(fitting_results[0])
            if 'background' in gc_plots.info.fit_info['FID'][peak]:
                globals()['bg_plot%s' %peak].set_ydata(int_tools.background_func[gc_plots.info.fit_info['FID'][peak]['background']](fitting_results[0],fitting_results[2],*fitting_results[3]))
            else:
                globals()['bg_plot%s' %peak].set_ydata(int_tools.errorfunc(fitting_results[0],fitting_results[2],*fitting_results[3]))
        except:
            pass
    fig.canvas.draw_idle()

def update_axis(val):
    f.axes.set_ylim([-0.1, val])
    fig.canvas.draw_idle()

GC_injection_slider.on_changed(update_injection)
Axis_slider.on_changed(update_axis)

###
#TCD
###

def integrated_data_plot(result_dict,error_dict,detector_key,title,ylabel,placement=111):
    f_2 = fig_2.add_subplot(placement)
    for entry in result_dict:
        f_2.plot(injectiontimes_sec/(60*60),result_dict[entry],ms=1, label=gc_plots.info.fit_info[detector_key][entry]['name'],color=gc_plots.info.fit_info[detector_key][entry]['color'])
        f_2.errorbar(injectiontimes_sec/(60*60),result_dict[entry],yerr=error_dict[entry],fmt='o',capsize=1, color=gc_plots.info.fit_info[detector_key][entry]['color'])
    f_2.set_xlabel('time [hrs]')
    f_2.set_ylabel(ylabel)
    f_2.tick_params('y')
    f_2.axes.set_title(title)
    g_2 = f_2.twinx()
    g_2.plot(rel_time/(60*60),reactor_temperature,'r--')
    g_2.set_ylabel('Temperature [celsius]',color='r')
    g_2.tick_params('y',colors='r')
    f_2.legend(frameon=0)

fig_2 = plt.figure()
integrated_data_plot(TCD_res,TCD_err,'TCD','GC TCD','Integrated TCD [Wb]',placement=211)
spectrum = spectrum = (sequence[0]['retention time']/60, sequence[0]['TCD'])   
f_2 = fig_2.add_subplot(212)
GC_line_2, = f_2.plot(spectrum[0],spectrum[1],'-',markersize=1.1)

for peak in gc_plots.info.fit_info['TCD']:
    fitting_results = int_tools.fit_to_spectrum(spectrum, peak, 'TCD')
    if 'background' in gc_plots.info.fit_info['TCD'][peak]:
        globals()['bg_plot%s' %peak], = f_2.plot(fitting_results[0],int_tools.background_func[gc_plots.info.fit_info['TCD'][peak]['background']](fitting_results[0],fitting_results[2],*fitting_results[3]))
    else:
        globals()['bg_plot%s' %peak], = f_2.plot(fitting_results[0],int_tools.errorfunc(fitting_results[0],fitting_results[2],*fitting_results[3]))
    
f_2.set_ylabel('Intensity [V]')
f_2.axes.set_ylim([-1, 10])

GC_injection_slider_2 = Slider(plt.axes([0.16, 0.035, 0.705, 0.03]), 'Time (hrs)', injectiontimes_sec[0]/(60*60), injectiontimes_sec[-1]/(60*60),valinit=injectiontimes_sec[0]/(60*60),valstep=diff_time_hrs)
Axis_slider_2 = Slider(plt.axes([0.16, 0.0, 0.705, 0.03]), 'Axis', -1, 10,valinit=1,valstep=0.001)

def update_injection_2(val):
    cache = int(round((val)*(1/np.mean(np.diff(injectiontimes_sec/(60*60))))+1))
    if int(cache) >= len(sequence):
        cache = len(sequence)-1
    spectrum = (sequence[int(cache)]['retention time']/60, sequence[int(cache)]['TCD'])
    GC_line_2.set_xdata(spectrum[0])
    GC_line_2.set_ydata(spectrum[1])
    for peak in gc_plots.info.fit_info['TCD']:
        try:
            fitting_results = int_tools.fit_to_spectrum(spectrum, peak, 'TCD')
            globals()['bg_plot%s' %peak].set_xdata(fitting_results[0])
            if 'background' in gc_plots.info.fit_info['TCD'][peak]:
                globals()['bg_plot%s' %peak].set_ydata(int_tools.background_func[gc_plots.info.fit_info['TCD'][peak]['background']](fitting_results[0],fitting_results[2],*fitting_results[3]))
            else:
                globals()['bg_plot%s' %peak].set_ydata(int_tools.errorfunc(fitting_results[0],fitting_results[2],*fitting_results[3]))
        except:
            pass
    fig_2.canvas.draw_idle()

def update_axis_2(val):
    f_2.axes.set_ylim([-1, val])
    fig_2.canvas.draw_idle()

GC_injection_slider_2.on_changed(update_injection_2)
Axis_slider_2.on_changed(update_axis_2)

"""
#END of script
"""

FID = 'FID'
TCD = 'TCD'
print('Program has finished')
print(' ')
print('-------------------GENERAL INFO-----------------')
print(' ')
print('* GC data is contained in the dictionaries:')
print('  FID_res; FID_err; TCD_res; TCD_err;')
print(' ')
print('* Gas-concentration vs time and vs temperature is contained in the dictionaries:')
print('  concentrations; concentration_vs_temp;')
print(' ')
print('* use "dict.keys()" in order to get all the entries')
print('  use "dict[key]", in order to write out the entry')
print(' ')
print('------------------------------------------------')

plt.show()

#import code
#code.interact(local=locals())
