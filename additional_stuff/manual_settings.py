"""
This script holds the information needed in 'CGCP' in order to treat data.
"""

"""
MANUEL SETTINGS SECTION
"""

class MANUAL_SETTINGS(object):

    def __init__(self):
        #settings
        self.skip_start = 1 #set to 1, and it will skip the first injection. The injection-times are unaltered
        self.skip_end   = 1   #set to 1, and it will skip the last injection
        self.delimiter  = '    ' #tabulated data

        #plot settings
        self.linear_plotted_gasses = ['CO', 'MeOH', 'CH4']#['CH4','C2H4','C2H6','C3H6','C3H8','C4Hx','MeOH','EtOH'] #possible gasses: H2, CO, CO2, Ar, CH4, MeOH, DME
        
        #MFC offset factors
        self.Ar_MFC_factor = 1.
        self.Ar_MFC_offset = 0.
        self.H2_MFC_factor = 1.
        self.H2_MFC_offset = 0.
        self.CO2_MFC_factor = 1.
        self.CO2_MFC_offset = 0.
        self.CO_MFC_factor = 1.
        self.CO_MFC_offset = 0.
        self.O2_MFC_factor = 1.
        self.O2_MFC_offset = 0.

        #Detector conversion factor
        
        self.GC_conversion_to_Perc = {}

        self.GC_conversion_to_Perc['calibration pressure'] = 1.0
        #If a pressure regulator is used, then a value in bars according to the pressure at the GC, should be entered (1.2 bars is usually a good guess)
        #If the pressure regulator is not used (it is not tightened), then leave the value as "None"

        #Info for CALIBRATION FACTORS. 
        #Following are calibrated on setup: H2, Ar, CO2, CO and CH4
        #calibration factors for other gasses can be found in [1] and [2]. 
        #[1] is the original calibration of the setup, and no other error should be included.
        #[2] is an estimate following the "effective carbon number", which should be used with caution, and with errors around 17%
        #
        #[1]: I. Sharafutdinov. "Investigations into low pressure methanol synthesis". PhD thesis, DTU March 2013
        #
        #[2]: K. Schofield. "The enigmatic mechanism of the flame ionization detector: Its overlooked implications for fossil fuel combustion modeling". In: Progress in Energy and Combustion Science 34 (2008), pp. 330-350.


        self.GC_conversion_to_Perc['TCD'] = {'H2'  : 314.0,#+-0.1,
                                             'Ar'  : 4.555,#+-0.001
                                             'CO2' : 3.701,#3.57,#+-0.001,
                                             'CO'  : 4.373}#4.28}#+-0.001}
                                             #'Methane': 6.2826}
        
        self.GC_conversion_to_Perc['FID'] = {'CH4' : 0.5158, #+-0.0031, or 1% error if used for other factors. Reponse is always 1
                                             'C2H4': 0.5158/1.97, #+-17%, response: 1.97 [2]
                                             'C2H6': 0.5158/1.99, #+-17%, response: 1.99 [2]
                                             'C3H6': 0.5158/3.00, #+-17%, response: 3 [2]
                                             'C3H8': 0.5158/3.00, #+-17%, response: 3 [2]
                                             'C4Hx': 0.5158/4.00,
                                             'DME' : 0.5158/1.10938, #+-1%; response: 1.10938 [1]
                                             'MeOH': 0.5158/0.628319, #+-1%; response: 0.628319 [1]
                                             'EtOH': 0.5158/1.59} #+-17%, response: 1.59 [2]

        
        #INFO FOR BACKGROUND FITTING AND INTEGRATION
        #
        #Names ('keys') and ordering of the peaks within the fit_info is abitrary, since fit_info is a dictionary.
        #The disignation 'P1', 'P2' or 'H2' will be the designated name of that gas in the output-data. Gasses that
        #are calibrated, are the only ones that are written to data. These are: H2, CO, CO2, Ar, CH4, MeOH, DME
        #
        #New peaks can be added by simply 'copy-paste' a section, containing a peak. It is important to have different
        #entries (P1,P2 etc.), else earlier entries with same name will be overwritten when integrating spectra.
        #
        #Peaks that need to be excluded, can be out-commented by # or sectional wise, like 'P2' in TCD.
        #
        #Compound-names are written in a latex-syntax.
        #
        #Colors are used in the plots
        
        self.fit_info = {}
        self.fit_info['settings'] = {}
        self.fit_info['settings']['FID'] = {}
        self.fit_info['settings']['FID']['background_range'] = 0.01 #minutes
        self.fit_info['settings']['TCD'] = {}
        self.fit_info['settings']['TCD']['background_range'] = 0.01 #minutes

        self.fit_info['FID'] = {}

        self.fit_info['FID']['CH4'] =  {'start'       : 0.74,
                                        'end'         : 1.17, 
                                        'name'        : 'CH$_4$/0.9 min',
                                        'color'       : 'orange',
                                        'mother_peak' : []} 

        self.fit_info['FID']['BKG'] =  {'start'       : 1.16,
                                        'end'         : 1.6, 
                                        'name'        : 'CO/CO$_2$ background',
                                        'color'       : 'black',
                                        'mother_peak' : []} 

        self.fit_info['FID']['C2H4'] =  {'start'      : 1.6,
                                        'end'         : 2.0, 
                                        'name'        : 'C$_2$H$_4$/1.7 min',
                                        'color'       : 'lightcoral',
                                        'mother_peak' : []} 

        self.fit_info['FID']['C2H6'] =  {'start'      : 2.0,
                                        'end'         : 2.5, 
                                        'name'        : 'C$_2$H$_6$/2.2 min',
                                        'color'       : 'red',
                                        'mother_peak' : []} 

        self.fit_info['FID']['C3H6'] =  {'start'      : 7.9,
                                        'end'         : 8.16, 
                                        'name'        : 'C$_3$H$_6$/8.0 min',
                                        'color'       : 'dodgerblue',
                                        'mother_peak' : []} 

        self.fit_info['FID']['C3H8'] =  {'start'      : 8.16,
                                        'end'         : 8.4, 
                                        'name'        : 'C$_3$H$_8$/8.2 min',
                                        'color'       : 'darkblue',
                                        'mother_peak' : []} 

        self.fit_info['FID']['DME'] =  {'start'       : 8.7,
                                        'end'         : 9.2, 
                                        'name'        : 'DME/8.8 min',
                                        'color'       : 'lightgreen',
                                        'mother_peak' : []} 

        self.fit_info['FID']['MeOH'] = {'start'       : 9.15,
                                        'end'         : 9.5, 
                                        'name'        : 'MeOH/9.3 min',
                                        'color'       : 'lightgreen',
                                        'mother_peak' : []} 

        self.fit_info['FID']['C4Hx'] = {'start'       : 9.6,
                                        'end'         : 10.2, 
                                        'name'        : 'C$_4$H$_x$/9.3 min',
                                        'color'       : 'purple',
                                        'mother_peak' : []} 

        self.fit_info['FID']['EtOH'] = {'start'       : 10.5,
                                        'end'         : 10.8, 
                                        'name'        : 'EtOH/10.6 min',
                                        'color'       : 'darkgreen',
                                        'mother_peak' : []} 

        self.fit_info['TCD'] = {}

        #self.fit_info['TCD']['Methane'] =  {'start'       : 2.9,
        #                                'end'         : 3.2,
        #                                'name'        : 'Methane/3.1 min',
        #                                'color'       : 'tab:orange',
        #                                'mother_peak' : [],
        #                                'background'        : 'error', #not mandatory
        #                                'background_range'  : 0.01} #not mandatory
        
        self.fit_info['TCD']['CO2'] =  {'start'       : 4.0,
                                        'end'         : 5.6,
                                        'name'        : 'CO$_2$/3.8 min',
                                        'color'       : 'tab:blue',
                                        'mother_peak' : []}

        self.fit_info['TCD']['H2'] =   {'start'       : 7.2,
                                        'end'         : 7.65,
                                        'name'        : 'H$_2$/7.4 min',
                                        'color'       : 'c',
                                        'mother_peak' : []}

        self.fit_info['TCD']['Ar'] =   {'start'       : 7.6,
                                        'end'         : 8.1,
                                        'name'        : 'Ar/7.8 min',
                                        'color'       : 'indigo',
                                        'mother_peak' : []}

#        self.fit_info['TCD']['O2'] =   {'start'       : 7.7,
#                                        'end'         : 8.0,
#                                        'name'        : 'O$_2$/7.8 min',
#                                        'color'       : 'indigo',
#                                        'mother_peak' : []}

        self.fit_info['TCD']['N2'] =   {'start'       : 8.0,
                                        'end'         : 8.6,
                                        'name'        : 'N$_2$/8.2 min',
                                        'color'       : 'm',
                                        'mother_peak' : []}

        self.fit_info['TCD']['CO'] =   {'start'             : 9.25,
                                        'end'               : 9.45,
                                        'name'              : 'CO/9.2 min',
                                        'color'             : 'saddlebrown',
                                        'mother_peak'       : [],
                                        'background'        : 'linear',
                                        'background_range'  : 0.01
                                        }