import os
import time
import re
import numpy as np
import sys
version = '2.1.0' #github version

def what_version():
    return version

#plot-design
import matplotlib as mpl
font = {'family' : 'serif',#'Arial',
        #'weight' : 'bold',
        'size'   : 18}
mpl.rc('font', **font)

from matplotlib import pyplot as plt
from scipy.special import erf
from scipy.optimize import curve_fit

import tkinter as tk
import tkinter.ttk as ttk
from tkinter.filedialog import askdirectory

"""
GC Tools
"""

class get_path(object):

    def ask_path(self):
        root = tk.Tk()
        style = ttk.Style(root) 
        style.theme_use("clam") 
        root.geometry("500x500")
        path = askdirectory(parent=root) 
        if sys.platform != 'darwin':
            root.destroy()
        return path

class data_loading_tools(object):

    def __init__(self,info):
        self.info = info

    def make_numbered_folder(self,folder_name):
        #will try to make a folder of said name, but extend with a number, if the folder already exists
        folder_made = False
        end = ''
        n=0
        while not folder_made:
            end='_'+str(n)
            try:
                directory = folder_name+end
                os.mkdir(directory)
                folder_made = True
                return directory
            except:
                n=n+1
    

    def read_injection_time(self,path):
        #read the injection time of the GC data file
        file=open(path,'r')
        lines=file.readlines()
        mask = ['#injection rel time: ' in i for i in lines]
        filtered_lines = np.array(lines)[mask]
        val = re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?",filtered_lines[0])[0]
        val = float(val)
        return val

    def read_line_pressure(self,path):
        #read the injection time of the GC data file
        file=open(path,'r')
        lines=file.readlines()
        mask = ['#line pressure: ' in i for i in lines]

        if list(np.array(lines)[mask]) == []:
            return 'unknown'

        else:
            filtered_lines = np.array(lines)[mask]
            val = re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?",filtered_lines[0])[0]
            val = float(val)
            return val
    
    def read_line_back_pressure(self,path):
        #read the injection time of the GC data file
        file=open(path,'r')
        lines=file.readlines()
        mask = ['#post GC pressure: ' in i for i in lines]

        if list(np.array(lines)[mask]) == []:
            return 'unknown'

        else:
            filtered_lines = np.array(lines)[mask]
            val = re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?",filtered_lines[0])[0]
            val = float(val)
            return val

    def fix_line_pressure(self,data):
        injection_pressure = data
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
        return injection_pressure
        
    def load_data_sheet(self,path, outcomment='#'):
        #This can open both .csv or tabulated data.
        #The data is stored in a numpy ndarray typed list, with dimensions matching the data
        file=open(path,'r')
        lines=file.readlines()
        mask = [outcomment not in i for i in lines]
        filtered_lines = np.array(lines)[mask]
        data_out = np.zeros(shape=(len(re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?",filtered_lines[0])),len(filtered_lines)))
        for row,line in enumerate(filtered_lines):
            values = re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?",line)
            for coloumn,value in enumerate(values):
                try:
                    data_out[coloumn,row]=float(value)
                except:
                    print('In text: '+str(filtered_lines)+'; of line: '+str(row)+'of the file: '+path+'; an error occurred')
        return data_out

    def Sequence(self,path):
        dir_list = os.listdir(path)
        liszt = sorted([j for j in dir_list if j[0:13] == 'GC_injection_' and j[-5] != ')' 
                                                 and not j[13:15].isdecimal()]) #first 9 injections
        cache_list = sorted([j for j in dir_list if j[0:13] == 'GC_injection_' 
                                                 and j[13:15].isdecimal() 
                                                 and not j[13:16].isdecimal()]) #next 90 injections
        liszt.extend(cache_list)
        cache_list = sorted([j for j in dir_list if j[0:13] == 'GC_injection_' and j[-5] != ')' 
                                                 and j[13:15].isdecimal() 
                                                 and j[13:16].isdecimal() 
                                                 and not j[13:17].isdecimal()]) #next 900 injections
        liszt.extend(cache_list)
        cache_list = sorted([j for j in dir_list if j[0:13] == 'GC_injection_' and j[-5] != ')'
                                                 and j[13:15].isdecimal() 
                                                 and j[13:16].isdecimal() 
                                                 and j[13:17].isdecimal()
                                                 and not j[13:18].isdecimal()]) #next 9000 injections, adds up to 2500 hours of experiment
        liszt.extend(cache_list)
    
        data_collection = []
        for file in liszt:
            data = self.load_data_sheet(path+'/'+file)
        
            data_dict                   = {}
            data_dict['injection time'] = self.read_injection_time(path+'/'+file)
            data_dict['line pressure']  = self.read_line_pressure(path+'/'+file)
            data_dict['line back pressure']  = self.read_line_back_pressure(path+'/'+file)
            data_dict['retention time'] = data[0]
            data_dict['FID']            = data[1]
            data_dict['TCD']            = data[2]
        
            data_collection.append(data_dict)
        
        
        
        data_collection = data_collection[self.info.skip_start:len(data_collection)-self.info.skip_end]
        return data_collection

    def load_setup_data(self,path):
        setup_data = self.load_data_sheet(path+'/'+[j for j in os.listdir(path) if j[0:10] == 'setup_data' and j[-4:] == '.txt'][0])

        #vectors: relative time(s), temperature before bed(celsius), pressure(bar), gas(H2), gas(CO2), gas(CO), gas(Ar), gas(O2), temperature after bed(celsius)
        rel_time = setup_data[0] #units in seconds
        reactor_temperature = setup_data[1] #units in celsius
        #try:
        #    reactor_pressure = (setup_data[7]+setup_data[8])*0.5 #units in bar
        #except:
        #    reactor_pressure = setup_data[7]
        reactor_pressure = setup_data[7] #front pressure is much closer to the actual reactor-pressure
        massflow_O2  = self.info.O2_MFC_factor*setup_data[2]+self.info.O2_MFC_offset #O2
        massflow_H2  = self.info.H2_MFC_factor*setup_data[3]+self.info.H2_MFC_offset #H2
        massflow_CO2 = self.info.CO2_MFC_factor*setup_data[4]+self.info.CO2_MFC_offset #CO2
        massflow_CO  = self.info.CO_MFC_factor*setup_data[5]+self.info.CO_MFC_offset #CO
        massflow_Ar  = self.info.Ar_MFC_factor*setup_data[6]+self.info.Ar_MFC_offset #Ar
        massflow_total = massflow_H2+massflow_CO2+massflow_CO+massflow_Ar+massflow_O2
        return rel_time, reactor_temperature, reactor_pressure, massflow_O2, massflow_H2, massflow_CO2, massflow_CO, massflow_Ar, massflow_total

class integration_tools(object):

    def __init__(self, info):
        self.info = info
        self.background_func = {'error': self.errorfunc,
                                'linear': self.linearfunc        
        }
    
    def errorfunc(self, time, center, size, height, slope):#step function
        step_speed = 10000 #making the erf-function an almost step-function
        return height+size*erf(step_speed*(time-center))

    def expfunc(self, time, center, size, height, slope):
        return height+np.exp(-slope*(time-size))-np.exp(0)#constrained exponential function

    def linearfunc(self, time, center, size, height, slope):
        return slope*(time-center)+height

    def fit_background(self, time_data, func_data, center, bkg_func = 'error'):
        time_data=np.array(time_data)
        func_data=np.array(func_data)
        if bkg_func == 'error':
            def background_function_1(time, size, height):
                return self.errorfunc(time, center, size, height, 0)
            best_fit = curve_fit(background_function_1,time_data,func_data)
            best_fit = [*best_fit[0], 0]
        elif bkg_func == 'linear':
            def background_function_2(time, height, slope):
                return self.linearfunc(time, center, 0, height, slope)
            best_fit = curve_fit(background_function_2,time_data,func_data)
            best_fit = [0, best_fit[0][0], best_fit[0][1]]
        else:
            def background_function_3(time, size, height):
                return self.errorfunc(time, center, size, height, 0)
            best_fit = curve_fit(background_function_3,time_data,func_data)
            best_fit = [*best_fit[0], 0]
        return best_fit

    def fit_exp_background(self, time_data, func_data, center):
        #this function is not used
        time_data=np.array(time_data)
        func_data=np.array(func_data)
        best_fit = curve_fit(self.expfunc,time_data,func_data,p0=[center,func_data[0],0.005],bounds=([time_data[0], 0, 0],[time_data[len(time_data)-1],5000.,2.]))#curvefit constrained exp function
        return best_fit

    def fit_to_spectrum(self, spectrum, peak, detector):
        mask = np.logical_and(spectrum[0]>self.info.fit_info[detector][peak]['start'],spectrum[0]<self.info.fit_info[detector][peak]['end'])
        peak_time = spectrum[0][mask]
        peak_count = spectrum[1][mask]
        if 'background_range' in self.info.fit_info[detector][peak]:
            background_mask_1 = np.logical_and(peak_time>self.info.fit_info[detector][peak]['start'],peak_time<self.info.fit_info[detector][peak]['start']+self.info.fit_info[detector][peak]['background_range'])
            background_mask_2 = np.logical_and(peak_time>self.info.fit_info[detector][peak]['end']-self.info.fit_info[detector][peak]['background_range'],peak_time<self.info.fit_info[detector][peak]['end'])
        else:
            background_mask_1 = np.logical_and(peak_time>self.info.fit_info[detector][peak]['start'],peak_time<self.info.fit_info[detector][peak]['start']+self.info.fit_info['settings'][detector]['background_range'])
            background_mask_2 = np.logical_and(peak_time>self.info.fit_info[detector][peak]['end']-self.info.fit_info['settings'][detector]['background_range'],peak_time<self.info.fit_info[detector][peak]['end'])
        peak_background_time  = np.append(peak_time[background_mask_1],peak_time[background_mask_2])
        peak_background_count = np.append(peak_count[background_mask_1],peak_count[background_mask_2])
        
        erf_center = [i for i, j in enumerate(peak_count) if j == np.amax(peak_count)] #the step of the error-function is placed at peak-maximum
        erf_center = peak_time[erf_center[0]] #value instead of a vector
        
        if 'background' in self.info.fit_info[detector][peak]:
            background_fit = self.fit_background(peak_background_time, peak_background_count, erf_center, bkg_func = self.info.fit_info[detector][peak]['background'])
        else:
            background_fit = self.fit_background(peak_background_time, peak_background_count, erf_center)
        return peak_time, peak_count, erf_center, background_fit
            
    def integrate_spectrum(self, injection, detector): 
        spectrum = (np.array(injection['retention time'])/60, np.array(injection[detector]))#fit background and integrate spectrum
                
        results_array = {}
        error_array = {}
    
        for peak in self.info.fit_info[detector]:
            fitting_results = self.fit_to_spectrum(spectrum, peak, detector)
            if 'background' in self.info.fit_info[detector][peak]:
                if self.info.fit_info[detector][peak]['background'] == 'error':
                    renorm_peak_count = np.array(fitting_results[1])-self.errorfunc(fitting_results[0],fitting_results[2],*fitting_results[3])     
                    area_error = (fitting_results[0][-1]-fitting_results[2])*fitting_results[3][0]*60 #60 is for seconds
                elif self.info.fit_info[detector][peak]['background'] == 'linear':
                    renorm_peak_count = np.array(fitting_results[1])-self.linearfunc(fitting_results[0],fitting_results[2],*fitting_results[3])     
                    area_error = ( ( 2*(fitting_results[0][-1]-fitting_results[2]) )**2 ) * fitting_results[3][-1]/2
                else:
                    renorm_peak_count = np.array(fitting_results[1])-self.errorfunc(fitting_results[0],fitting_results[2],*fitting_results[3])     
                    area_error = (fitting_results[0][-1]-fitting_results[2])*fitting_results[3][0]*60 #60 is for seconds
            else:
                renorm_peak_count = np.array(fitting_results[1])-self.errorfunc(fitting_results[0],fitting_results[2],*fitting_results[3])     
                area_error = (fitting_results[0][-1]-fitting_results[2])*fitting_results[3][0]*60 #60 is for seconds
                
            area = (renorm_peak_count*np.mean(np.diff(fitting_results[0]))*60).sum() #60 is for seconds
        
            results_array[peak] = abs(area)
            error_array[peak] = abs(area_error)
        
        for peak in self.info.fit_info[detector]:
            if 'mother_peak' in self.info.fit_info[detector][peak]:
                for motherpeak in self.info.fit_info[detector][peak]['mother_peak']:
                    results_array[motherpeak] = results_array[motherpeak]-results_array[peak]
        else:
            return results_array, error_array

    def write_to_txt(self,filename, data, dimension='', data_x = [None]):
        if filename[-4:len(filename)] == '.txt':
            file = open(filename,'w+')
        else: 
            file = open(filename+'.txt','w+')
        if dimension == '':
            str_cache = '#'
        else:
            str_cache = '#'+dimension+self.info.delimiter
        for key in list(data.keys()):
            str_cache = str_cache+str(key)+self.info.delimiter

        file.write(str_cache+'\n')
        
        if data_x[0] == None:
            str_cache = ''
            for n in range(0,len(data[key])):
                for key in list(data.keys()):
                    str_cache = str_cache+str(data[key][n])+self.info.delimiter
                file.write(str_cache+'\n')
                str_cache = ''
        else:
            for t,n in zip(data_x,range(0,len(data[key]))):
                str_cache = str(t)+self.info.delimiter
                for key in list(data.keys()):
                    str_cache = str_cache+str(data[key][n])+self.info.delimiter
                file.write(str_cache+'\n')
        file.close()

"""
all the plots
"""

class GC_plots(object):
    
    def __init__(self, info):
        self.info = info
    
        self.fig_info = {} #Make a dict for graph-colors and such
        for detector in ['TCD','FID']:
            for gas in list(self.info.fit_info[detector].keys()):
                self.fig_info[gas] = {}
                self.fig_info[gas]['name'] = gas
                self.fig_info[gas]['color'] = self.info.fit_info[detector][gas]['color']

    def smooth(self,f,pts):
        extent = np.ones(pts)/pts
        smoothed_f = np.convolve(f, extent, mode='same')
        return smoothed_f
    
    def raw_data_heat_plot(self, sequence):
        #make time-axis for plot of raw data (not retention time)
        first_injection = sequence[0+self.info.skip_start]['injection time']
        last_injection = sequence[-1-self.info.skip_end]['injection time']
        injection_spacing = (last_injection-first_injection)/len(sequence) #average injection spacing, this might not be true

        #setup of raw data
        #spectrum_length = len(sequence[0]['retention time'])-500 #hardcoded to exclude the last 10 points, since the spectra-length vary this value will be the reference for the length, and generally the sampling is high, so a GC-spectrum could be 10000-100000 points long
        spectrum_length = min([len(i['FID']) for i in sequence])-1
        TCD_raw = []
        FID_raw = []
        for injection in sequence: 
            TCD_raw.append(injection['TCD'][0:spectrum_length])
            FID_raw.append(injection['FID'][0:spectrum_length])
        TCD_raw = np.vstack(TCD_raw)
        FID_raw = np.vstack(FID_raw)
        
        retention_time = sequence[0]['retention time'][0:spectrum_length]

        def temp_func(data,title): #instead of writing out all the fancy plot, a function has been made (saving 7 lines of code)
            data = np.array(data)
            data[1] = data[1]+0.3
            fig, f = plt.subplots()
            img = f.imshow(abs(data),
                           aspect= 'auto',
                           extent= [sequence[0]['retention time'][0]/60,
                                    sequence[0]['retention time'][spectrum_length]/60,
                                    last_injection/3600,
                                    first_injection/3600],
                           norm  = mpl.colors.LogNorm(),
                           cmap  = 'hot')
            f.axes.set_title(title)
            f.set_xlabel('retention time [min]')
            f.set_ylabel('Experiment time [hrs]')
            fig.colorbar(img, ax=f)
            plt.tight_layout()
            plt.draw()

        temp_func(FID_raw,'FID raw data')
        temp_func(TCD_raw,'TCD raw data')

    def MFC_plot(self, rel_time, reactor_temperature, reactor_pressure, massflow_O2, massflow_H2, massflow_CO2, massflow_CO, massflow_Ar):
        smooth_parameter=15
        
        fig, f = plt.subplots()
        f.plot(rel_time/(3600),self.smooth(massflow_H2,smooth_parameter),label='H$_2$')
        f.plot(rel_time/(3600),self.smooth(massflow_CO2,smooth_parameter),label='CO$_2$')
        f.plot(rel_time/(3600),self.smooth(massflow_CO,smooth_parameter),label='CO')
        f.plot(rel_time/(3600),self.smooth(massflow_Ar,smooth_parameter),label='Ar')
        f.plot(rel_time/(3600),self.smooth(massflow_O2,smooth_parameter),label='O$_2$')
        f.legend(frameon=0)
        f.set_xlabel('Time [hours]')
        f.set_ylabel('Flow [ml/min]')
        f.set_title('MFC flows')
        
        g = f.twinx()
        g.plot(rel_time/(3600),reactor_temperature,'r--')
        g.set_ylabel('Temperature [celsius]',color='r')
        g.tick_params('y',colors='r')

        h = f.twinx()
        h.spines["right"].set_position(("axes",1.2))
        h.set_frame_on(True)
        h.patch.set_visible(False)
        for sp in h.spines.values():
            sp.set_visible(False)
        h.spines["right"].set_visible(True)
        h.plot(rel_time/(3600),self.smooth(reactor_pressure,smooth_parameter),'-o',color='b',markersize=0.5)
        h.set_ylabel('Pressure [bar]',color='b')
        h.tick_params('y',colors='b')
        
        plt.tight_layout()
        plt.draw()

    def integrated_vs_time(self, injectiontimes_sec, result_dict,error_dict,detector_key,title,rel_time,reactor_temperature,reactor_pressure):
        fig, f = plt.subplots()
        for entry in result_dict:
            f.plot(injectiontimes_sec/3600,
                   result_dict[entry],
                   ms   = 1, 
                   label= self.info.fit_info[detector_key][entry]['name'],
                   color= self.info.fit_info[detector_key][entry]['color'])
            f.errorbar(injectiontimes_sec/3600,
                       result_dict[entry],
                       yerr   = error_dict[entry],
                       fmt    = 'o',
                       capsize= 1, 
                       color  = self.info.fit_info[detector_key][entry]['color'])
        f.set_xlabel('time [hrs]')
        f.set_ylabel('Integrated detector [Wb]')
        f.tick_params('y')
        f.axes.set_title(title)
        f.legend(frameon=0)

        g = f.twinx()
        g.plot(rel_time/(3600),reactor_temperature,'r--')
        g.set_ylabel('Temperature [celsius]',color='r')
        g.tick_params('y',colors='r')

        #h = f.twinx()
        #h.spines["right"].set_position(("axes",1.2))
        #h.set_frame_on(True)
        #h.patch.set_visible(False)
        #for sp in h.spines.values():
        #    sp.set_visible(False)
        #h.spines["right"].set_visible(True)
        #h.plot(rel_time/(3600),reactor_pressure,'-o',color='b',markersize=0.5)
        #h.set_ylabel('Pressure [bar]',color='b')
        #h.tick_params('y',colors='b')

        plt.tight_layout()
        plt.draw()
    
    def log_concentration_vs_time(self,injectiontimes_sec,concentrations,gas_list,rel_time,reactor_temperature):
        break_statement = False
        for key in gas_list:
            if key not in concentrations.keys():
                print('KeyWarning: '+key+' not found in calibrated data. Check the "linear_plotted_gasses" array in "manual_settings.py".')
                break_statement = True
        if break_statement:
            print('PlotWarning: Errors found in "manual_settings.py". Skipping the log plot calibrated data.')
            return

        fig, f = plt.subplots()
        for result in gas_list:
            f.semilogy(injectiontimes_sec/(60*60),
                       concentrations[result],ms=1, 
                       label = self.fig_info[result]['name'],
                       color = self.fig_info[result]['color'])
            f.semilogy(injectiontimes_sec/(60*60),
                       concentrations[result],
                       'v',
                       ms    = 3, 
                       color = self.fig_info[result]['color'])
        f.set_xlabel('Time [hours]')
        f.set_ylabel('Volumetric Concentration [%]')
        f.tick_params('y')
        f.legend(loc=0,frameon=0)
        g = f.twinx()
        g.plot(rel_time/(60*60),reactor_temperature,'r--')
        g.set_ylabel('Temperature [celsius]',color='r')
        g.tick_params('y',colors='r')
        plt.tight_layout()
        plt.draw()
        
    def lin_concentration_vs_temp(self,concentration_vs_temp):
        collection_of_markers = ['o','*','s','^','v','<','>','o','*','s','^','v','<','>','o','*','s','^','v','<','>','o','*','s','^','v','<','>']
        cool_downs = np.where(np.diff(concentration_vs_temp['temperature']) <= -6)[0]
        if len(cool_downs) > 3 or len(cool_downs) == 0:
            cool_downs = [len(concentration_vs_temp['temperature'])]
        else:
            cool_downs = [*cool_downs, len(concentration_vs_temp['temperature'])]

        break_statement = False
        for key in self.info.linear_plotted_gasses:
            if key not in concentration_vs_temp.keys():
                print('KeyWarning: '+key+' not found in calibrated data. Check the "linear_plotted_gasses" array in "manual_settings.py".')
                break_statement = True
        if break_statement:
            print('PlotWarning: Errors found in "manual_settings.py". Skipping the linear concentration vs temperature plot.')
            return
        
        fig, f = plt.subplots()        
        for start_index,end_index in zip([-1,*cool_downs[:-1]],cool_downs):
            for gas,marker in zip(self.info.linear_plotted_gasses,collection_of_markers):
                f.errorbar(concentration_vs_temp['temperature'][start_index+1:end_index+1],
                           concentration_vs_temp[gas][start_index+1:end_index+1],
                           xerr = concentration_vs_temp['temperature_std'][start_index+1:end_index+1],
                           yerr = [concentration_vs_temp[gas+'_herr'][start_index+1:end_index+1],
                                   concentration_vs_temp[gas+'_lerr'][start_index+1:end_index+1]],
                           fmt  = marker,capsize=2,
                           color= self.fig_info[gas]['color'])
                f.plot(concentration_vs_temp['temperature'][start_index+1:end_index+1],
                       concentration_vs_temp[gas][start_index+1:end_index+1],
                       label = self.fig_info[gas]['name'],
                       color = self.fig_info[gas]['color'])
        f.set_xlabel('Temperature [celsius]')
        f.set_ylabel('Volumetric Concentration [%]')
        f.legend(loc=0,frameon=0)
        plt.tight_layout()
        plt.draw()
"""
stupid animation
"""

def CowGuin():
    if os.name == 'nt':
        clear = lambda: os.system('cls')
    else:
        clear = lambda: os.system('clear')
    time.sleep(0.9)
    clear()
    print(' ')
    print('    ')
    print('   ')
    print('   ')
    print('     ')
    print('   ')
    print('        .--.             ^__^')
    print('       |o_o |            (oo)\_______      ')
    print('       |:_/ |            (__)\       )\    ')
    print('      //   \ \               ||----w | *   ')
    print('     (|     | )              ||     ||     ')
    print('    /´\_   _/`\ ')
    print('    \___)=(___/          * ')
    print('                        \|/ ')
    time.sleep(0.8)
    clear()
    print(' ')
    print('  ______________________')
    print('< Nice code for your GC! >      ')
    print(' ------------------------     ')
    print('   \                            ')
    print('    \                          ')
    print('        .--.             ^__^')
    print('       |o_o |            (oo)\_______     ')
    print('       |:_/ |            (__)\       )\    ')
    print('      //   \ \               ||----w | *  ')
    print('     (|     | )              ||     ||    ')
    print('    /´\_   _/`\ ')
    print('    \___)=(___/          * ')
    print('                        \|/ ')
    time.sleep(1.1)
    clear()
    print(' ')
    print('    ')
    print('   ')
    print('   ')
    print('     ')
    print('   ')
    print('        .--.             ^__^')
    print('       |o_o |            (oo)\_______      ')
    print('       |:_/ |            (__)\       )\    ')
    print('      //   \ \               ||----w | *   ')
    print('     (|     | )              ||     ||     ')
    print('    /´\_   _/`\ ')
    print('    \___)=(___/          * ')
    print('                        \|/ ')
    time.sleep(0.8)
    clear()
    print(' ')
    print('    ')
    print('   ')
    print('   ')
    print('     ')
    print('   ')
    print('        .--.             ^__^           /       ')
    print('       |o_o |            (oo)\_______  /_´_       ')
    print('       |:_/ |            (__)\       )\       ')
    print('      //   \ \               ||----w | *      ')
    print('     (|     | )              ||     ||        ')
    print('    /´\_   _/`\ ')
    print('    \___)=(___/          * ')
    print('                        \|/ ')
    time.sleep(0.8)
    clear()
    print(' ')
    print('  ')
    print('                                 ______________________')
    print('                                < Beware of methane... >')
    print('                                 ----------------------')
    print('                                    /')
    print('        .--.             ^__^           ')
    print('       |o_o |            (oo)\_______         ')
    print('       |:_/ |            (__)\       )\       ')
    print('      //   \ \               ||----w | *      ')
    print('     (|     | )              ||     ||        ')
    print('    /´\_   _/`\ ')
    print('    \___)=(___/          * ')
    print('                        \|/ ')
    time.sleep(1.1)
    clear()

def GcFridayLogo():
    print(' FRIDAYFRIDAYFRIDAYFRIDAYFRIDAYFRIDAYFRIDAYFRIDAYFRIDAYFRIDAYFRIDAYFRIDAYFRIDAYFRIDAYFRIDAYFRIDAY ')
    print('     ______      ______      ______       __       _______      _______  _________   _______      ')
    print('   .´ ___  |   .´ ___  |    |_   __ \    /  \     |_   __ \    /  ___  ||_   ___  | |_   __ \     ')
    print('  / .´   \_|  / .´   \_|      | |__) |  / /\ \      | |__) |  |  (__ \_|  | |_  \_|   | |__) |    ')
    print('  | |    ____ | |             |  ___/  / ____ \     |  __ /    `.___`-.   |  _|  _    |  __ /     ')
    print('  \ `.___]  _|\ `.___.´\     _| |_   _/ /    \ \_  _| |  \ \_ |`\____) | _| |___/ |  _| |  \ \_   ')
    print('   `._____.´   `._____.´    |_____| |____|  |____||____| |___||_______.´|_________| |____| |___|  ')
    print(' FRIDAYFRIDAYFRIDAYFRIDAYFRIDAYFRIDAYFRIDAYFRIDAYFRIDAYFRIDAYFRIDAYFRIDAYFRIDAYFRIDAYFRIDAYFRIDAY ')  
    print(' For use at DTU Physics                                                    GC Parser '+version+' (2022) ')
 
def GcLogo():
    print('----------------------------------------------------------------------------------------------')
    print('   ______      ______      ______       __       _______      _______  _________   _______    ')
    print(' .´ ___  |   .´ ___  |    |_   __ \    /  \     |_   __ \    /  ___  ||_   ___  | |_   __ \   ')
    print('/ .´   \_|  / .´   \_|      | |__) |  / /\ \      | |__) |  |  (__ \_|  | |_  \_|   | |__) |  ')
    print('| |    ____ | |             |  ___/  / ____ \     |  __ /    `.___`-.   |  _|  _    |  __ /   ')
    print('\ `.___]  _|\ `.___.´\     _| |_   _/ /    \ \_  _| |  \ \_ |`\____) | _| |___/ |  _| |  \ \_ ')
    print(' `._____.´   `._____.´    |_____| |____|  |____||____| |___||_______.´|_________| |____| |___|')
    print('----------------------------------------------------------------------------------------------')  
    print('For use at DTU Physics                                                  GC Parser '+version+' (2022)')

def GcFridaySpecial():
    if os.name == 'nt':
        clear = lambda: os.system('cls')
    else:
        clear = lambda: os.system('clear')
    question = input('Is it friday?(y/n): ')
    if question != 'n':
        if question != 'y':
            print("Wait, what? I'll take that as an yes.")
        time.sleep(1.1)
        clear()
        GcFridayLogo()
        time.sleep(1.1)
        clear()
        GcFridayLogo()
        CowGuin()
        clear()
        GcFridayLogo()    
    elif question == 'n':
        print('.. hmm... My calendar says so...')
        question = input('Are you sure?(y/n): ')
        if question == 'y':
            clear()
            GcLogo()
        elif question != 'y':
            if question != 'n':
                print("You don't sound sure to me.")
            else:
                print('Yes, I thought so to.')
            time.sleep(1.1)
            clear()
            GcFridayLogo()
            time.sleep(1.1)
            clear()
            GcFridayLogo()
            CowGuin()
            clear()
            GcFridayLogo()
