
----------------------------------------------------------------------------------------------
   ______      ______      ______       __       _______      _______  _________   _______    
 .' ___  |   .' ___  |    |_   __ \    /  \     |_   __ \    /  ___  ||_   ___  | |_   __ \   
/ .'   \_|  / .'   \_|      | |__) |  / /\ \      | |__) |  |  (__ \_|  | |_  \_|   | |__) |  
| |    ____ | |             |  ___/  / ____ \     |  __ /    `.___`-.   |  _|  _    |  __ /   
\ `.___]  _|\ `.___.'\     _| |_   _/ /    \ \_  _| |  \ \_ |`\____) | _| |___/ |  _| |  \ \_ 
 `._____.'   `._____.'    |_____| |____|  |____||____| |___||_______.'|_________| |____| |___|
----------------------------------------------------------------------------------------------

For use at DTU Physics                                                    GC Parser 2.0.0 (2021)

CGCP - Continuous Gas Chromatograph Parser

(+) Initial notes
    This readme file has been written for the updated version of CGCP, on 5th of March 2021.
    This particular version of CGCP is given with scripting for synchronizing GC data simultaneously 
with gas-flow, temperature and pressure data from a test reactor. This will not work for any 
other setup than the plug-flow reactor placed in HetCat lab. If this code is needed for another 
setup, please contact Thomas (thoe@fysik.dtu.dk)
    This code has been written to work well, launched from a commando-prompt (so not spyder or 
any sort). Wether it works well on Spyder is up to the user, and (s)he is free to rewrite the 
code for this purpose if needed. 
    The logo is from before version 1.1, where the very general name "GC parser" was used.

(+) Content
- A main script called "CGCP.py". This is the actual script to run.
- A folder called "additional_stuff", which contains:
	- The "Readme.txt" file, you are currently reading
	- A script called "manual_settings.py" containing most hard-coded settings, so coding in the
		main-script can be avoided.
	- A script called "GcFunc.py" which contains a lot of functions that the main-script uses.
	- An article containing relative FID-factors. This is need to have, when new hydrocarbons 
		are seen in the gas-flow.
- A folder called "test_data". This can be used to illustrate the capabilities of this program. This data has not been tested with versions of the program later than 1.5.0

(+) Requirements
- Python 3.6+ distribution with the following packages: numpy, matplotlib, scipy, unicodecsv, Tkinter
    - An easy way to find out, which packages are missing, is to simply run the script and get the errors.
- That your data is ordered in the right way. The folder of "test_data" can be explored for that reason.

(+) Running
- 	Make sure that "CGCP.py" is in a folder, and that there is a folder called "additional_stuff" that contains "GcFunc.py" and "manual_settings.py".
	- 	Alternatively you can place "additional_stuff" in your python-path, to make these functions available always.

- 	Set the values in "additional_stuff/manual_settings.py", which is a script that is executed from the main script. Most importantly but not exclusively the following:
		line 18, which contains a list of gasses of interest.
		line 51 to 62, contains calibration factors for the GC-detectors with expected uncertainties.
		line 83 and 85, the range of background before and after a peak. If it's to loong, peaks might overlap. If it's to short, the fitting might be unstable. It's better to go for less than for more.
		line 80 and out, the fitting parameters. Pay good attetion to this!

-    A new peak-entry into the fitting-dictionary is made up from the following template:

    self.fit_info['TCD']['N2'] =   {'start'       : 8.0, #integration start
                                    'end'         : 8.6, #integration stop
                                    'name'        : 'N$_2$/8.2 min', #name in plots. Takes latex-syntax
                                    'color'       : 'm', #color to use in matplotlib
                                    'mother_peak' : []} #name of another overlapping peak, that should have the integration of this peak subtracted

-    Additionally, you can add the following optional lines within the dictionary entry, {}:

                                    'background'        : 'linear', #background to use. "linear" and "error" is implemented, and "error" is the standard step-function background
                                    'background_range'  : 0.01 #the amount of time in each end of the integration-range, that the background should be fitted to

- 	Run "CGCP.py" with python from a commando-prompt (fx. Anaconda Prompt, or any prompt 
	that can launch python 3.6+). The syntax would be "python CGCP.py" or "python3 CGCP.py", if
	your current path is that, which contains "CGCP.py".

-	A gui appears for choosing the directory containing the data. Choose your head-directory containing 
	all your data.

-   	If a folder is already present in the work directory named by today's date, CGCP will make a new
	folder extended with "_#", where # is the number of times, this happened.

- 	Figures will open in windows themselves. The matplotlib GUI can be used to zoom in on the figures 
	or save them in a number of formats. Zooming in the gui is used by drawing square to zoom in, and
	use left-arrow to zoom out.

-   	The last generated figures contain interactive scalebars for surveying the raw data. Check that 
	your data is integrated with a proper background. If not, either move or extend the background-
	limits in the "manual_settings.py" script.

-	The generated data will be written in the directory named by the today's date. There
    	you'll find .txt files and .p files. ".p"-files are files used with pickle. They are loaded
    	within python using "pickle.load(open, "full_path_to_file","rb")", and saved to a variable. 
    	Use the variable-name followed by ".keys()" to inspect the dictionary, which contains all the data. 
    	ALL DATA IS SAVED IN PERCENTAGE! 

(+) Calibration of gasses
-	Calibration factors of a variety of gasses is already present in "manuel_settings.py".
- 	Calibration factors on different gasses can be updated in the "manual_settings.py". For more info, 
	on how to calibrate GC-detectors, read the appendices of the master-thesis "Synthesis and 
	characterization of novel catalysts: Heusler alloy nanoparticles" by Thomas Smitshuysen. This should
	be available from DTU Findit.
- 	Relative calibration factors vs methane for the FID-detector can be found in the attached article.
	- not all alkanes/alkenes/alkynes are present. Either find another article or use the rule of thumbs: the 		relative response to methane goes as 1 over every CHx-group. Thus CH4 has factor 1, C2H4 and C2H6 has 		1/2, C8H18 has 1/8.
-	There is also the possibility of changing factors and off-sets for data achieved from mass-flow 
	controllers (MFC's). This should be used, when gasses used in the experiment doesn't have the right
	flow-factors. Brooks and/or Bronkhorst have tables available for conversion-factors.

(+) Description of what the program does
- Loads GC, mass-flow, temperature and pressure data.
- numerically integrate each GC-spectrum according to the settings in the "manuel_settings.py"
	- at this point, data is written to a .txt file with time-resolved uncalibrated integrated GC-data
- synchronise the GC, mass-flow, temperature and pressure data, as the time-resolution is different
	- at this point, data is written to a .txt file with time-resolved data
- Tests in how many steps the temperature is changing with the following requirement: temperature is only changing, when the change is larger than 4 degrees celsius
- Order the time-resolved data in temperature-steps with two requirements: temperature has to have changed at least 4 degrees, and the first data-point after a temperature-change is never included.
	- at this point, data is written to a .txt file with temperature-resolved data
- Generate plots, where integrated data and raw GC-data can be viewed simultaneously
- all plots are printed to the user. These can be saved, zoomed or interacted with.
	- A good habit is NOT to use auto-generated figures, when making illustrations for a text. It is much better 	to find the data that is to be illustrated, and make the figure accordingly (and manually).

(+) Wish-list for future features of the program
- A GUI for changing the "manuel_settings.py" and possibility of running the main-features from the program within this GUI.
- Peak-finding algorithm to propose peaks
- real fitting of the peak-curvature and deconvolution of overlapping peaks

(+) Known unresolved bugs
- When running in OSX, the data-selection-window will not disappear and behaves weirdly odd. This has no planned solution at this time, because it is a bug within the Tkinter-package when running OSX.
- active figures doesn't work properly when using Spyder. An expert in spyder could tell you, how to set the settings properly in spyder in order to allow the figures that are generated in this script

(+) Help and questions
Write to thoe@fysik.dtu.dk or find me in 307:044 at DTU. 

(+) Thanks
To Kenneth Nielsen, one of the great contributors on earth. Especially towards proper python-programming.

=========================================================================================

Version history of gc_parser.py

0.1.0
- Load temperature and gc-data
- full raw data rainbow plot
- FID and TCD-signal plot based on area integration from Agilent software

0.2.0 (Faulty)
- fitting using "skewed gaussian”. By curtesy of Simon.

0.2.1
- fitting using "exponential modelled gaussian”.

0.2.2
- added constraints of fit.
- added better start-guesses for fit.

0.3.0
- numerical integration.
- background fitted with erf-function (or erfc)
- multiple spectra could be treated in sequence.
- FID and TCD-signal plot based on numerical integration

0.4.0
- final written code
- functions and fit_info moved to separate script
- Added constrained exponential background, for when peaks in spectrum were placed in close vicinity. “Constrained” in the sense that it has a lesser degree of freedom but is more numerical stable.

1.0
- identical in many ways to version 0.4.0.
- included GC-calibrated values of FID and TCD detector for plotting real concentrations

1.1
- fixed time
- stamps for GC-data, so they are loaded from metadata instead of being generated
- changed temperature-dependant measurements, so point taken at temperature
- changes beyond 4 degrees are excluded

Version history of CGCP.py

1.2
- Major rewriting to shorten the code.
- fixed time for real, compared to 1.1
- changed loaded functions to actual imported functions
- moved often changed settings to a separate script

1.3
- It is unclear, what changed for this version

1.4.0
- Wrote in a GUI for choosing your data
- Data are now also saved in .txt-files
- Included oxygen from the labview data
- Reduced the amount of information in "manual_settings"
- Updated the logo with year and version

1.4.1
- Data integrated from the FID and TCD are now written in .txt as well
- Solved error-message issue with the interactive slider-plots
- Moved the setting of "skipping data" into the "manual_settings" again

1.5.0
- Cleaned up the scripting. a lot. Now python-classes are used! :D
- Written the code for differently formatted data: this is a result of the loss of a very expensive license-key, written on stupid cardboard.
- This is the version released for our new raspberry-pi powered agilent GC

1.5.1
- cleaned some minor mistakes of the 1.5-release.
- Wrote the code locally, so a new version was needed to keep track.

1.5.2
- another clean-up of more bugs
- improved stability, for when data doesn't fit well together
- improved some of the plotting

1.5.3
- New calibration factors for FID and TCD has been included properly
- Better error handling

1.5.4
- a few bugs. nothing big.
- added the pressure-readout from the newly mounted pressurecontroller.

1.5.5
- resolved a few bugs
- better read-out of data to .txt

1.5.6
- made the script compatible with all data later than version 1.5

1.6.0
- This version is motivated by a conversation with Mikkel Stensberg
- Added the possibility of choosing backgrounds for the integration
	- alongside the step-function there is now also a linear background
	- this comes with the possibility of setting specific fitting-ranges for each peak
- additional fitting-error data wrote out to .txt files
- solved wrong units on raw data
- cleaned up the code. A lot

2.0.0
- Copy of 1.6.0, but made a github repository instead.
