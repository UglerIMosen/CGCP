import os

class settings_generator(object):

    def __init__(self,path,manual_settings):
        self.path = path
        self.output = manual_settings
    
    def write_file(self):
        try:
            self.file = open(self.path+'/manual_settings.py','w+')
        except:
            os.mkdir(self.path)
            self.file = open(self.path+'/manual_settings.py','w+')
        self.file.write('"""\nThis autogenerated script holds the information needed for "CGCG" in order to treat data.\n"""\n\n"""\nMANUAL SETTINGS SECTION\n"""\n\n')
        self.file.write('class MANUAL_SETTINGS(object)\n')
        self.file.write('    def __init__(self):\n        #settings\n')
        self.file.write('        self.skip_start = '+str(self.output.skip_start)+' #set to 1, and it will skip the first injection. The injection-times are unaltered\n')
        self.file.write('        self.skip_end = '+str(self.output.skip_end)+' #set to 1, and it will skip the last injection\n')
        self.file.write('        self.delimiter = "'+self.output.delimiter+'"\n\n')
        self.file.write('        self.linear_plotted_gasses = [')
        for gas in self.output.linear_plotted_gasses:
            self.file.write('"'+gas+'"')
            if gas != self.output.linear_plotted_gasses[-1]:
                self.file.write(', ')
            else:
                self.file.write(']\n\n')
        self.file.write('        #MFC offset factors\n')
        self.file.write('        self.Ar_MFC_factor = '+str(self.output.Ar_MFC_factor)+'\n')
        self.file.write('        self.Ar_MFC_offset = '+str(self.output.Ar_MFC_offset)+'\n')
        self.file.write('        self.H2_MFC_factor = '+str(self.output.H2_MFC_factor)+'\n')
        self.file.write('        self.H2_MFC_offset = '+str(self.output.H2_MFC_offset)+'\n')
        self.file.write('        self.CO2_MFC_factor = '+str(self.output.Ar_MFC_factor)+'\n')
        self.file.write('        self.CO2_MFC_offset = '+str(self.output.Ar_MFC_factor)+'\n')
        self.file.write('        self.CO_MFC_factor = '+str(self.output.Ar_MFC_factor)+'\n')
        self.file.write('        self.CO_MFC_offset = '+str(self.output.Ar_MFC_factor)+'\n')
        self.file.write('        self.O2_MFC_factor = '+str(self.output.Ar_MFC_factor)+'\n')
        self.file.write('        self.O2_MFC_offset = '+str(self.output.Ar_MFC_factor)+'\n\n')
        self.file.write('        #Detector conversion factor\n\n')
        self.file.write('        self.GC_conversion_to_Perc = {}\n\n')
        self.file.write('        self.GC_conversion_to_Perc["calibration pressure"] = '+str(self.output.GC_conversion_to_Perc['calibration pressure'])+'\n\n')
        self.file.write('        #Info for CALIBRATION FACTORS. \n')
        self.file.write('        #Following are calibrated on setup: H2, Ar, CO2, CO and CH4\n')
        self.file.write('        #calibration factors for other gasses can be found in [1] and [2]. \n')
        self.file.write('        #[1] is the original calibration of the setup, and no other error should be included.\n')
        self.file.write('        #[2] is an estimate following the "effective carbon number", which should be used with caution, and with errors around 17%\n')
        self.file.write('        #\n')
        self.file.write('        #[1]: I. Sharafutdinov. "Investigations into low pressure methanol synthesis". PhD thesis, DTU March 2013\n')
        self.file.write('        #\n')
        self.file.write('        #[2]: K. Schofield. "The enigmatic mechanism of the flame ionization detector: Its overlooked implications for fossil fuel combustion modeling". In: Progress in Energy and Combustion Science 34 (2008), pp. 330-350.\n\n')
        self.file.write('        self.GC_conversion_to_Perc["TCD"] = {\n')
        for gas in list(self.output.GC_conversion_to_Perc["TCD"]):
            self.file.write('                                            "'+gas+'" : '+str(round(self.output.GC_conversion_to_Perc["TCD"][gas],5)))
            if gas != list(self.output.GC_conversion_to_Perc["TCD"])[-1]:
                self.file.write(',\n')
            else:
                self.file.write('}\n\n')
        self.file.write('        self.GC_conversion_to_Perc["FID"] = {\n')
        for gas in list(self.output.GC_conversion_to_Perc["FID"]):
            self.file.write('                                            "'+gas+'" : '+str(round(self.output.GC_conversion_to_Perc["FID"][gas],5)))
            if gas != list(self.output.GC_conversion_to_Perc["FID"])[-1]:
                self.file.write(',\n')
            else:
                self.file.write('}\n\n')
        self.file.write('        #INFO FOR BACKGROUND FITTING AND INTEGRATION\n')
        self.file.write('        #\n')
        self.file.write('        #Names ("keys") and ordering of the peaks within the fit_info is abitrary, since fit_info is a dictionary.\n')
        self.file.write('        #The disignation "P1", "P2" or "H2" will be the designated name of that gas in the output-data. Gasses that\n')
        self.file.write('        #are calibrated, are the only ones that are written to data. These are: H2, CO, CO2, Ar, CH4, MeOH, DME\n')
        self.file.write('        #\n')
        self.file.write('        #New peaks can be added by simply "copy-paste" a section, containing a peak. It is important to have different\n')
        self.file.write('        #entries (P1,P2 etc.), else earlier entries with same name will be overwritten when integrating spectra.\n')
        self.file.write('        #\n')
        self.file.write('        #Peaks that need to be excluded, can be out-commented by # or sectional wise, like "P2" in TCD.\n')
        self.file.write('        #\n')
        self.file.write('        #Compound-names are written in a latex-syntax.\n')
        self.file.write('        #\n')
        self.file.write('        #Colors are used in the plots\n\n')
        self.file.write('        self.fit_info = {}\n')
        self.file.write('        self.fit_info["settings"] = {}\n')
        self.file.write('        self.fit_info["settings"]["FID"] = {}\n')
        self.file.write('        self.fit_info["settings"]["FID"]["background_range"] = '+str(self.output.fit_info["settings"]["FID"]["background_range"])+' #minutes\n')
        self.file.write('        self.fit_info["settings"]["TCD"] = {}\n')
        self.file.write('        self.fit_info["settings"]["TCD"]["background_range"] = '+str(self.output.fit_info["settings"]["TCD"]["background_range"])+' #minutes\n\n')        
        for detector in ['FID','TCD']:
            self.file.write('        self.fit_info["'+detector+'"] = {}\n\n')
            for peak in list(self.output.fit_info[detector]):
                self.file.write('        self.fit_info["'+detector+'"]["'+peak+'"] = {\n')
                for item in list(self.output.fit_info[detector][peak]):
                    if isinstance(self.output.fit_info[detector][peak][item],int) or isinstance(self.output.fit_info[detector][peak][item],float) or item == 'mother_peak':
                        self.file.write('                                        "'+item+'" : '+str(self.output.fit_info[detector][peak][item]))
                    else:
                        self.file.write('                                        "'+item+'" : "'+str(self.output.fit_info[detector][peak][item])+'"')
                    if item != list(self.output.fit_info[detector][peak])[-1]:
                        self.file.write(',\n')
                    else:
                        self.file.write('}\n\n')
        
        self.file.close()

if __name__ == '__main__':
    from manual_settings import MANUAL_SETTINGS
    manual_settings = MANUAL_SETTINGS()
    func = settings_generator('dummy',manual_settings)
    func.write_file()