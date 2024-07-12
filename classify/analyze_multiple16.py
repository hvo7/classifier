"""
Script Name: analyze_multiple16.py -- modified from make_spectrum.py
Description:  Program to read either 1 day's data or full year, accumulate and record classification data (latitude, 
     time to Fhv transition, solar angle, FWHM, t90), write to summary fFle, and make plots 
     (time profile, energy spectra, latitude, solar activity).
Note:  Horizontal and vertical plot limits are set (and can be adjusted) in set_plot_limit
Version  calculates GOES avvg and variance for +/- 1 hour aroune trigger time.
Recalculates time from hv on/off transition -- 6/17/24

Usage: python ./analyze_multiple14.py [fileid(yyyymmdd) or "Dec2023" for example to analyze multiple days ] chan_type(PHA/PI) plot_duration (+/-, optional)
goes
Current version uses old version of get_phdata for fileid > 20230719 and Yuta's new version for fileid < 20230720. (5/10/2023)
analyze_multple2 deletes individual plot files after they are merged together.

Incorporates Henry's GOES solar activity data routine for +/- 1 day.

Read event file and prepare list of previous classifications for each trigger number
If len(arg[1]) != 8, then read summary file; for each line, read fileid and call analyzexx to do normal analysis
if len(arg[1]) == 8, call analyzexx to do normal analysis only
    For example, Dec2023 calls summary_file_Dec2023.txt to analyze data for 12/15, 12/18, 12/31/2023:
          summary_file_20240430_2348
          remote_file_path = /mnt/cgbmdata2/cgbm_heasarc_archive/cgbmarchL3HEA_corrected/obs/2023/20230110/auxil/cgbm_20230110.att
          local_file_path = ./cgbml3_data/cgbm_20230110.attcgbm_20230110.att  
          1  20231215  [1386668718, 1386688282, 1386693962, 1386698214][755948846.3019712, 755968411.1089857, 755974091.3486664, 755978343.5158461]
          1  20231218  [1386944516, 1386950104, 1386955813][756224650.259189, 756230238.4677762, 756235946.6890068]
          1  20231231  [1388093748, 1388094024, 1388094280][757373905.3055743, 757374181.3115886, 757374437.3316015]

          Job duration = -1273.84 mins = -76431 secs.
    Note: First 3 lines are blank, end list of data files with a blank line.

For each trigger, write Prev_classification based on trigger number. (Trigger numbers may differ from list at Yoshidalab.mydns.jp by 1.)

Call solar only if trig_num = 0
Status as of 3/28/24:
     Correctly handles multiple triggers per day
     Needs to work on alternate plot input values, Norris backgound fitting
Program (as of 3/21/2024): 
     Moves all Global statements to start of main routine
     Includes t90
     Corrects fitting parameter printouts
     Writes name of program to output file
     Incorporate Norris fitting function
     determines ref_MET time
     calculates time to closest hv turn on/off
     plots counts vs time 
     fits background vs time 
     plots energy spectra for 6 detectors
     Hardness ratio
     determine latitude
     determine angle to Sun
     FWHM, t90 of time distribution
     Allow for multiple triggers on a single day 
     checks for triggers on a specified day
     chi-sq of fit
     spectra_y_or_n removed, effectively set to "y"
     Save plots and output data in summary files
     Include likelihood calculations 
     Incorporates Henry's solar activity program.
     Version incorporates find_triggers routi‪ne to find MET time.
     This script determines time difference between trigger time and hv turn on/off.
     Looks for hv turn-on/turn-off (HXM1L time change) > threshold counts (def = 2000) between times -2000 and 2000 sec from trigger.
     Prints out time difference.
     User has option to plot spectra or not
Still:
     Check with Yuta about hardness ratio energy ranges; need to account for background
     Clean up redundant ouput and redundant plot1
     Clean up outputs and pdfmerger
     Add weighting for likelihoods
     Add lots of comments
     Choose L or 1-L for each parameter
     Program frequently hangs at solar call. If I run again, it usually works fine.
     Code only prints out counts around hv turn on/off for first trigger of the day. 

Author: Mike Cherry
Date: 5/10/2024
Version: 2
"""

import sys
import datetime
from astropy.table import Table, vstack
import numpy as np
import cgutil as cg
import matplotlib.pyplot as plt
#from pypdf import PdfMerger
from pypdf import PdfWriter
# from matplotlib.backends.backend_pdf import PdfPages
from dotenv import load_dotenv
import os
from statistics import mean
import math
import pandas as pd
from scipy.optimize import curve_fit

from astropy.table import Table
from astropy.time import Time
from astropy.time import TimeDelta
import astropy.units as u
from scipy.spatial.transform import Rotation as R
from astropy.coordinates import SkyCoord, UnitSphericalRepresentation
from astropy.coordinates import get_sun
import cartopy.crs as ccrs

#Input line
args = sys.argv
if len(args) < 3 :
    print("Usage: python ./analyze_multiple2.py fileid(yyyymmdd) chan_type(PHA/PI) spectrum duration (+/-, optional)")
    sys.exit(1)

else:
    fileid = str(args[1])
    chan_type = str(args[2])
 
if len(args) == 4:
    tmin = -int(args[3])
    tmax = int(args[3])

else:
    tmin= -100
    tmax = 100

###############################################################
###############################################################

def main (fileid, chan_type, tmin, tmax):

# Read event file and prepare list of previous classifications for each trigger number
# Read in prior event classifications
    global prior_trig_id, prior_trig_class
    event_class = open ("Event_classification.txt", "r")  	
    lines = len(event_class.readlines())  		# Determine length of file
    event_class.seek(0)					# Rewind to beginning
    content = []
    category = []  
    prior_trig_id = []
    prior_trig_class = []
    for i in range (0, lines):	
        data = event_class.readline().strip()		# Read one line at a time, strip out EOF	
        if (not data) or (len(data) == 0):
            break
        content.append (data)			
        category = content[i].split()     # Category = data for a single line split into individual variables
        if category[3] == 'T':
            prior_trig_id.append (category[1])		# prior_trig_id = trigger id
            prior_trig_class.append (category[5])	# prior_trig_class = previously assigned event classification       
        else:
            prior_trig_id.append (category[1])		# prior_trig_id = trigger id
            prior_trig_class.append ('Unknown')		# prior_trig_class = previously assigned event classification       
    event_class.close()

    """
    Prepare looping -- one day or many
    if len(fileid) == 8, call analyzexx to do normal analysis only
    If len(fileid) != 8, then read run_list; for each line, read fileid and call analyzexx to do normal analysis
    For example, Dec2023 calls run_list_Dec2023.txt to analyze data for 12/15, 12/18, 12/31/2023:
          run_list_20240430_2348
          remote_file_path = /mnt/cgbmdata2/cgbm_heasarc_archive/cgbmarchL3HEA_corrected/obs/2023/20230110/auxil/cgbm_20230110.att
          local_file_path = ./cgbml3_data/cgbm_20230110.attcgbm_20230110.att  
          1  20231215  [1386668718, 1386688282, 1386693962, 1386698214][755948846.3019712, 755968411.1089857, 755974091.3486664, 755978343.5158461]
          1  20231218  [1386944516, 1386950104, 1386955813][756224650.259189, 756230238.4677762, 756235946.6890068]
          1  20231231  [1388093748, 1388094024, 1388094280][757373905.3055743, 757374181.3115886, 757374437.3316015]

          Job duration = -1273.84 mins = -76431 secs.
    Note: First 3 lines are blank, end list of data files with a blank line.
    """
# Analyze list of events
# if len(fileid) == 8, call analyze39 to do normal analysis on single day only
    if len(fileid) == 8:		
        analyze39 (fileid, chan_type, tmin, tmax)
        sys.exit()

# if len(fileid) != 8, call analyze39 to do normal analysis on all days in "summary_file"+fileid+".txt"
    if len(fileid) != 8:			
#        summary_file_list = open ("summary_file_"+fileid+".txt", "r")
        run_list = open ("run_list_"+fileid+".txt", "r")
        content = []
        category = []  
        for i in range (0,3):				# Skip 3 header lines
            data = run_list.readline().strip()	
        i_event = 0
        while i_event > -1:
            data = run_list.readline().strip()	
            content.append (data)			
            category = content[i_event].split()     # Category = data for a single line split into individual variables
            if len(category) == 0:
                i_event = -1
                break
            fileid = category [1]
            i_event += 1
            try:
                analyze39 (fileid, chan_type, tmin, tmax)
            except:
                continue
        run_list.close()

###############################################################

def analyze39 (fileid, chan_type, tmin, tmax):
    # Input example
    #fileid = '20231110'
    #chan_type = 'PHA'
    #ref_MET = 752966385.210472

#    global tmin, tmax
#    global summary_file, summary_file_name
    global data_types, event_classes, max_rows, max_bins
    global xL, Ngtx

# Set up local disk directory environment
    load_dotenv()
    data_dir = os.getenv("LOCAL_DATA_DIR")

# Set up figure saving to file
    plt.rcParams["figure.figsize"] = [7.00, 3.50] 
    plt.rcParams["figure.autolayout"] = True

# Account for format change in mid-2023
    if int(fileid) > 20230725:		
        ph_table = get_phdata(fileid, data_dir, chan_type)
    else:
        ph_table = get_phdata_modified(fileid, data_dir, chan_type)

################################
    """
# Read in likelihood tables -- old version
    likelihood_fileid = '20230615'	# file name for likelihood data; this can be changed if necessary
    input_file_name = "ranking_" + likelihood_fileid + ".txt"

#    global data_types, event_classes, max_rows

    max_rows = 1000		# Max number of rows; this can be changed if necessary

    data_types = ['Latitude','Time_to_transition','Solar_angle','Gaussian_FWHM','Gaussian_Amplitude', 'Gaussian_Chi-sq','Norris_t90',
        'Norris_Amplitude', 'Norris_Chi-sq', 'Hardness_ratio', 'GOES', 'GOES_variance', 'GOES_avg_1hr', 'GOES_var_1hr']				# Data types

    event_classes = ['GRB', 'Solar','Particle']				# Event classes
#    likelihood_list = [lat, transition, angle_to_Sun, FWHM, t90]

#  Create 3d arrays
#    global xL, Ngtx
    xL =   [[ ['0' for row in range(max_rows)] for col in range(len(event_classes))] for type in range(len(data_types))]	        #  x
    Ngtx =   [[ ['0' for row in range(max_rows)] for col in range(len(event_classes))] for type in range(len(data_types))]	#  N(>x)

    rows = ranking (max_rows, input_file_name)

    norm = [0] * len(event_classes)			# Normalize Ngtx to 1 at row_unmber = 0
    for itype in range (len(data_types)):
       for events in range (len(event_classes)):
           norm[events] =  (Ngtx [itype][events][0])
           for row_number in range (rows [itype]):
#                print ('$$$$$', itype, events, row_number, norm, Ngtx [itype][events][row_number] , type (Ngtx [itype][events][row_number]) )
                Ngtx [itype][events][row_number] = Ngtx [itype][events][row_number] / norm [events]

#    for itype in range (len(data_types)):
#        print (data_types [itype],':')  
#        for row_number in range (rows [itype]):
#            print (row_number, end = '     ')
#            for events in range (len(event_classes)):		
#                print (xL [itype][events][row_number], end = '  ')
#            print ('  ',end='   ')
#            for events in range (len(event_classes)):     
#                print (  "{:.3f}".format( Ngtx [itype][events][row_number] ), end = '  ')
#            print()
    """

# Read in likelihood tables -- new version

    global data_types, event_classes, midpoints, likelihood
    likelihood_fileid = '20240602'	# file name for likelihood data; this can be changed if necessary
    input_file_name = "ranking_" + likelihood_fileid + ".txt"
    #print(likelihood_fileid)
	
    data_types = ['Latitude','Time_to_transition','Solar_angle','Gaussian_FWHM','Gaussian_Amplitude', 'Gaussian_Chi-sq','Norris_t90',
        'Norris_Amplitude', 'Norris_Chi-sq', 'Hardness_ratio', 'GOES', 'GOES_variance', 'GOES_avg_1hr', 'GOES_var_1hr']				# Data types
    event_classes = ['GRB', 'Solar','Particle']				# Event classes


    lhood_input = open ( input_file_name, 'r')
    lines = len(lhood_input.readlines())  		# Determine length of file
    lhood_input.seek(0)					# Rewind to beginning
    max_bins = int((lines - 1)/len(data_types))	- 1	# Number of histogram bins must be the same for all data types

    likelihood = [[ [0 for row in range(max_bins)] for col in range(len(event_classes))] for type in range(len(data_types))]
    midpoints = [[ 0 for bin_num in range(max_bins)] for type in range(len(data_types))]

    data = lhood_input.readline().strip()			# Skip header in first line of input file				

    for type in range (len(data_types)):
        category = []  
        data = lhood_input.readline().strip()			# Read data type				
        category = data.split()	
        if category [0] != data_types [type]:
            print ('Wrong data type: ', category[0])
            sys.exit()
        for bin in range (max_bins):
            category = []
            data = lhood_input.readline().strip()		# Read data 				
            category = data.split()	
            midpoints [type] [bin] = category [0]
            likelihood [type] [0] [bin] = category [1]
            likelihood [type] [1] [bin] = category [2]
            likelihood [type] [2] [bin] = category [3]
#            print (midpoints [type] [bin], likelihood [type] [index] [bin] for index in range (3))

################################
# Read in event list with previous characterizations.
    max_rows = 2000
    event_file_name = "event_classification.txt"
#    event_id_rows, event_UT_rows, event_type_rows = []
    event_id_rows, event_UT_rows, event_type_rows = events_list (max_rows, event_file_name)

################################
# find MET time. How about if there are multiple files on a single day?
    trig_index = []
    trig_time = []
    triggers, trig_index, trig_time = determine_MET_time (fileid, data_dir)
#    print ('Result of finding triggers:')
#    print (triggers)
#    print (trig_index)
#    print (trig_time)
#    num = input()   

################################
# Step through triggers
    for trig_num in range (len(trig_index)):
        xyz = analyze_single_trigger (fileid, chan_type, data_dir, ph_table, trig_index, trig_time, trig_num, triggers,
        event_id_rows, event_UT_rows, event_type_rows, tmin, tmax)
    
    return

##########################################################################################
##########################################################################################

# Read in likelihood tables
def ranking (max_rows, input_file_name):
    """
Description:  Program to read in, categorize likelihood tables. Based on ranking4.py.
Program reads ranking_yyyymmdd.txt tab-delimited data file of likelihood data for 
    event classes, uses 3d list.
    """

    rows = [0] * len(data_types)

    input_file = open(input_file_name,"r")
    
    for line in range (0,4):                     # Read and skip over header lines
        header = input_file.readline().strip()  

# Create 3d arrays
#    global xL, Ngtx
#    xL =   [[ ['#' for row in range(max_rows)] for col in range(len(event_classes))] for type in range(len(data_types))]	#  x
#    Ngtx =   [[ ['#' for row in range(max_rows)] for col in range(len(event_classes))] for type in range(len(data_types))]	#  N(>x)

    read_more_data = 'y'
    while read_more_data == 'y':
        content = []    			# Content = Each line of Content is a raw data line from spreadsheet
        category = []  				# Category = data for a single line split into N(x), N(>x), etc
        data = input_file.readline().strip()	# Check for end of file
        if (not data) or (len(data) == 0):
#            print('\nEnd of data')
            break
        content.append (data)			
        category = content[0].split()
#        print ('\nNew data set: ',category[0])
        type = 'None'
        for itype in range(len(data_types)):		
            if category[0] == data_types[itype]:
                type = category[0]
                break
        if type == 'None':
            print("Data_type must be one of types listed as available in code: ", data_types)
            sys.exit(1)

        header = input_file.readline().strip()	# Read in potential event class (GRB, solar, particle) 
        content.append (header)
        event_class = content[1].split()
        if (event_class != event_classes):
            print ("Available event classes: ", event_class,'\n',event_classes)
            sys.exit(2)

        content = []  				# Read in data; first initialize content array to zero again
        row_number = 0				# Count row_number in range 0 to row_number - 1
        while True:					# Read in data line by line
            data = input_file.readline().strip()
            content.append (data)
            category = content[row_number].split()
            end_indicator = category[0]
            if (end_indicator == "***"):
                break 

            for j in range(len(event_classes)):					# GRB, solar, or particle ?
                xL [itype][j][row_number] = eval (category [4 + j])
#                Ngtx [itype][j][row_number] = category[7 + j]    		# Use integral distribution
                Ngtx [itype][j][row_number] = eval (category[13 + j])		# Use differentiated integral distribution

            row_number += 1				# row_number is now actual count of rows that have been read in
            rows [itype] = row_number

            if row_number > max_rows:
                print ('Maximum number of rows (',max_rows,' has been exceeded. max_rows must be increased.')
                sys.exit()
        continue
    input_file.close()
 
    return rows

###############################################################  
def events_list (max_rows, event_file_name):
    event_file = open(event_file_name,"r")

    event_id_rows = []
    event_UT_rows = []
    event_type_rows = []
    content = []    			# Content = Each line of Content is a raw data line from spreadsheet
    category = []  				# Category = data for a single line split into N(x), N(>x), etc
    for ievent_num in range (max_rows):
        data = event_file.readline().strip()	# Check for end of file
        if (not data) or (len(data) == 0):
            break		
        category = data.split()
        if len(category) < 6:
            data = event_file.readline().strip()	# Check for end of file
            if (not data) or (len(data) == 0):
                break		
            category = data.split()
        event_id_rows.append(category[1])
        event_UT_rows.append(category[2])
        if category[5] == 'LGRB' or category[5] == 'SGRB' or category[5] == 'PGRB':
            category[5] = event_classes[0]			# GRB
        if category[5] == 'SFLR':				# solar flare
            category[5] = event_classes[1]
        if category[5] == 'PART':				# Particle
            category[5] = event_classes[2]
        event_type_rows.append(category[5])
#        print (event_id_rows[ievent_num], event_UT_rows[ievent_num], event_type_rows[ievent_num])  
    event_file.close()						
#    print ('Length:  ',len(event_id_rows))
    return event_id_rows, event_UT_rows, event_type_rows

###############################################################
def determine_MET_time (fileid,data_dir):
#     Returns triggers, trig_index,trig_time

    hk_file = '%s/cgbm_%s.hk' % (data_dir, fileid)
    hk_table = Table.read (hk_file, hdu=1)

    # Timestamps of CGBM data are end of the bins.
    # Timestamps were changed to center of the bins
    hk_table['TIME'] = hk_table['TIME'] - cg.PERIODIC_HBIN
    hk_table['MDCTIME'] = hk_table['MDCTIME'] - cg.PERIODIC_HBIN

    # The transition of TRIG_STATUS (from 0 to 1) means onboard trigger happened.
    hk_table['TRIG_TRANSITION'] = 0
    hk_table['TRIG_TRANSITION'][0] = 0
    hk_table['TRIG_TRANSITION'][1:] = hk_table['TRIG_STATUS'].data[1:] - hk_table['TRIG_STATUS'].data[:-1]

    # I don't know why this is needed but TRIG_TRANSITION can be unexpected values.
    bad_mask = (hk_table['TRIG_TRANSITION'] == 0) | (hk_table['TRIG_TRANSITION'] == 1)
    hk_table['TRIG_TRANSITION'][~bad_mask] = 0

    trig_index = hk_table['TRIG_TRANSITION'] == 1

    # Extracting rows of which TRIG_TRANSITION == 1
    triggers = hk_table[trig_index]

    # Counting and Printing the number of triggers
    print('\n%d triggers in %s ' % (len(triggers), fileid))

    trig_index = []
    trig_time = []
    # Printing the TriggerID and Trigger time
    indexp = 0
    for row in triggers:
        indexp +=1
        print('TriggerID = %d   Trigger Time (MET) = %.6f   Trigger Time (MDC) = %.6f' 
            % (int(row['MDCTIME']), float(row['TIME']), float(row['MDCTIME'])))
        trig_index.append (int(row['MDCTIME']))
        trig_time.append (float(row['TIME']))
    print ()

    return triggers, trig_index,trig_time

###############################################################
def get_phdata(fileid, data_dir, chan_type):
#  Returns merged
    if not (chan_type == 'PHA' or chan_type == 'PI'):
        print('chan_type is not expected.')
        return None

    hxm1h_ph_file = '%s/cgbm_%s_hx1_hph.fits' % (data_dir, fileid)
    hxm1h_ph_table = Table.read(hxm1h_ph_file, hdu=1)
    hxm1l_ph_file = '%s/cgbm_%s_hx1_lph.fits' % (data_dir, fileid)
    hxm1l_ph_table = Table.read(hxm1l_ph_file, hdu=1)

    hxm2h_ph_file = '%s/cgbm_%s_hx2_hph.fits' % (data_dir, fileid)
    hxm2h_ph_table = Table.read(hxm2h_ph_file, hdu=1)
    hxm2l_ph_file = '%s/cgbm_%s_hx2_lph.fits' % (data_dir, fileid)
    hxm2l_ph_table = Table.read(hxm2l_ph_file, hdu=1)

    sgmh_ph_file = '%s/cgbm_%s_sgm_hph.fits' % (data_dir, fileid)
    sgmh_ph_table = Table.read(sgmh_ph_file, hdu=1)
    sgml_ph_file = '%s/cgbm_%s_sgm_lph.fits' % (data_dir, fileid)
    sgml_ph_table = Table.read(sgml_ph_file, hdu=1)

    # 'TIME' is Misson Elapsed Time (MET) same as Suzaku/MAXI
    # MET = 0 is 2000-01-01T00:00:00 UTC
    # 'MDCTIME' is time in MDC time
    merged = hxm1h_ph_table[['TIME', 'MDCTIME']]   
    merged = hxm1l_ph_table[['TIME', 'MDCTIME']]

    # '*_EXPOSURE's are exposure for each detector.
    merged['HXM1_EXPOSURE'] = hxm1h_ph_table['EXPOSURE']
    merged['HXM1_EXPOSURE'] = hxm1l_ph_table['EXPOSURE']

    merged['HXM2_EXPOSURE'] = hxm2h_ph_table['EXPOSURE']
    merged['SGM_EXPOSURE'] = sgmh_ph_table['EXPOSURE']
#    print (sgmh_ph_table['PHA'])
#    merged
#    print (merged['SGM_EXPOSURE'])
#    input()

    # '*_COUNTS's are counts per PH time binning (4s)
    merged['HXM1H_COUNTS'] = hxm1h_ph_table[chan_type]
    merged['HXM1L_COUNTS'] = hxm1l_ph_table[chan_type]
    merged['HXM2H_COUNTS'] = hxm2h_ph_table[chan_type]
    merged['HXM2L_COUNTS'] = hxm2l_ph_table[chan_type]
    merged['SGMH_COUNTS'] = sgmh_ph_table[chan_type]
    merged['SGML_COUNTS'] = sgml_ph_table[chan_type]

    merged['HXM1H_RATE'] = merged['HXM1H_COUNTS'] / np.tile(merged['HXM1_EXPOSURE'], (cg.PH_BINNUM_HIGH, 1)).T
    merged['HXM2H_RATE'] = merged['HXM2H_COUNTS'] / np.tile(merged['HXM2_EXPOSURE'], (cg.PH_BINNUM_HIGH, 1)).T
    merged['SGMH_RATE'] = merged['SGMH_COUNTS'] / np.tile(merged['SGM_EXPOSURE'], (cg.PH_BINNUM_HIGH, 1)).T
    merged['HXM1L_RATE'] = merged['HXM1L_COUNTS'] / np.tile(merged['HXM1_EXPOSURE'], (cg.PH_BINNUM_LOW, 1)).T
    merged['HXM2L_RATE'] = merged['HXM2L_COUNTS'] / np.tile(merged['HXM2_EXPOSURE'], (cg.PH_BINNUM_LOW, 1)).T
    merged['SGML_RATE'] = merged['SGML_COUNTS'] / np.tile(merged['SGM_EXPOSURE'], (cg.PH_BINNUM_LOW, 1)).T

    return merged
###############################################################

def get_phdata_modified(fileid, data_dir, chan_type):
    if not (chan_type == 'PHA' or chan_type == 'PI'):
        print('chan_type is not expected.')
        return None

    hxm1h_ph_file = '%s/cgbm_%s_hx1_hph.fits' % (data_dir, fileid)
    hxm1h_ph_table = Table.read(hxm1h_ph_file, hdu=1)
    hxm1l_ph_file = '%s/cgbm_%s_hx1_lph.fits' % (data_dir, fileid)
    hxm1l_ph_table = Table.read(hxm1l_ph_file, hdu=1)

    hxm2h_ph_file = '%s/cgbm_%s_hx2_hph.fits' % (data_dir, fileid)
    hxm2h_ph_table = Table.read(hxm2h_ph_file, hdu=1)
    hxm2l_ph_file = '%s/cgbm_%s_hx2_lph.fits' % (data_dir, fileid)
    hxm2l_ph_table = Table.read(hxm2l_ph_file, hdu=1)

    sgmh_ph_file = '%s/cgbm_%s_sgm_hph.fits' % (data_dir, fileid)
    sgmh_ph_table = Table.read(sgmh_ph_file, hdu=1)
    sgml_ph_file = '%s/cgbm_%s_sgm_lph.fits' % (data_dir, fileid)
    sgml_ph_table = Table.read(sgml_ph_file, hdu=1)

    # 'TIME' is Misson Elapsed Time (MET) same as Suzaku/MAXI
    # MET = 0 is 2000-01-01T00:00:00 UTC
    # 'MDCTIME' is time in MDC time
    merged = hxm1h_ph_table[['TIME', 'MDCTIME']]

    # '*_EXPOSURE's are exposure for each detector.
    merged['HXM1_EXPOSURE'] = hxm1h_ph_table['EXPOSURE']
    merged['HXM2_EXPOSURE'] = hxm2h_ph_table['EXPOSURE']
    merged['SGM_EXPOSURE'] = sgmh_ph_table['EXPOSURE']

    # '*_COUNTS's are counts per PH time binning (4s)
    chan_type = 'COUNTS_' + chan_type
    merged['HXM1H_COUNTS'] = hxm1h_ph_table[chan_type]
    merged['HXM1L_COUNTS'] = hxm1l_ph_table[chan_type]
    merged['HXM2H_COUNTS'] = hxm2h_ph_table[chan_type]
    merged['HXM2L_COUNTS'] = hxm2l_ph_table[chan_type]
    merged['SGMH_COUNTS'] = sgmh_ph_table[chan_type]
    merged['SGML_COUNTS'] = sgml_ph_table[chan_type]

    merged['HXM1H_RATE'] = merged['HXM1H_COUNTS'] / np.tile(merged['HXM1_EXPOSURE'], (cg.PH_BINNUM_HIGH, 1)).T
    merged['HXM2H_RATE'] = merged['HXM2H_COUNTS'] / np.tile(merged['HXM2_EXPOSURE'], (cg.PH_BINNUM_HIGH, 1)).T
    merged['SGMH_RATE'] = merged['SGMH_COUNTS'] / np.tile(merged['SGM_EXPOSURE'], (cg.PH_BINNUM_HIGH, 1)).T
    merged['HXM1L_RATE'] = merged['HXM1L_COUNTS'] / np.tile(merged['HXM1_EXPOSURE'], (cg.PH_BINNUM_LOW, 1)).T
    merged['HXM2L_RATE'] = merged['HXM2L_COUNTS'] / np.tile(merged['HXM2_EXPOSURE'], (cg.PH_BINNUM_LOW, 1)).T
    merged['SGML_RATE'] = merged['SGML_COUNTS'] / np.tile(merged['SGM_EXPOSURE'], (cg.PH_BINNUM_LOW, 1)).T

    return merged

###############################################################
def set_plot_limits (x, yarray):
    xmin = -100
    xmax = 100
    ymin = 1000000000
    ymax = 0
    for i in range (len(x)):
        if x [i] >= xmin :
            if x[i] > xmax:
                break
            yi = np.sum(yarray, axis = 1)[i]
            if yi < ymin:
                ymin = yi
            if yi > ymax:
                ymax = yi
    ymin = int(ymin/500)*500
    ymax = (int(ymax/500) + 1) * 500
    return xmin, xmax, ymin, ymax

###############################################################
# Routine to plot counts vs time     
def time_dist(x, y_HXM1L, y_HXM2L, y_SGML, ref_MET, fileid, trig_num):
    fig, cv = plt.subplots(3, 1, sharex=True, figsize=(8, 6))
    #    cv[0].set_xlim(-1800, 1800)
    xmin, xmax, ymin, ymax = set_plot_limits(x, y_HXM1L)
    cv[0].set_xlim(xmin,xmax)
    cv[0].set_ylim(ymin,ymax)
    cv[0].plot(x, np.sum(y_HXM1L, axis=1), drawstyle='steps-mid', label='HXM1')
    xmin, xmax, ymin, ymax = set_plot_limits(x, y_HXM2L)
    cv[1].set_ylim(ymin,ymax)
    cv[1].plot(x, np.sum(y_HXM2L, axis=1), drawstyle='steps-mid', label='HXM2')
    xmin, xmax, ymin, ymax = set_plot_limits(x, y_SGML)
    cv[2].set_ylim(ymin,ymax)
    cv[2].plot(x, np.sum(y_SGML, axis=1), drawstyle='steps-mid', label='SGM')
 
    id = 'Trigger #' + fileid + '.' + str(trig_num)
    cv[2].set_xlabel('Time [s] since %.6f -- %s' % (ref_MET, id), fontsize=14)
    plt.subplots_adjust(hspace=0, top=0.95, bottom=0.1, left=0.15, right=0.95)

    for ax in cv:
        ax.axvline(x=0, color='red', linestyle='--')
        ax.tick_params(labelsize=14)
        ax.set_ylabel('Counts/s', fontsize=14)
        plt.legend(loc='best')

#    print (' \nClose plot to continue.')
#    plt.savefig ('Plot1.pdf')
#    plt.show()

    return

###############################################################
def fitter(x_data, y_data, ref_MET, tmin, tmax):
# Version of  .                   for Gaussian function
#   Returns ampl, peak_pos, bkg, s_over_rootB, FWHM, chisq

    average = mean(np.sum(y_data, axis = 1))
    peak = max(np.sum(y_data, axis = 1))

# generate data over restricted range tmin to tmax defined as global variables at start of main routine
    x_data0 = []		# x_data, y_data are original data
    y_data0 = []		# x_data0, y_data0 are data defined in restricted range

    ij = 0
    for i in range (0,len(x_data)):
        if ((x_data[i] >= tmin) & (x_data[i] <= tmax)):
            x_data0.append (x_data[i])
            y_data0.append (np.sum(y_data, axis = 1)[i])
            ij += 1
#            print ('Gauss:  ', xdata0[i]. ydata0[i])

    params = [0.,0.,20., 2., .03, .0003]		# Fitting parameters
    a = ampl = peak - average 
    b = peak_pos = params[1] 
    c = FWHM = params[2] 
    d = bkg1 = average 
    e = bkg_linear = params[4]
    f = bkg_quad = params[5]

    y_fit = fit_function (x_data0, a, b, c, d, e, f)  

# curve_fit() function takes the test-function x-data0 and y-data0 as arguments and returns 
# the coefficients a - f in param and the estimated covariance of param in param_cov
 
    p0 = [ a, b, c, d, e, f ]				# Initial values of fitting parameters
    try:
        param, param_cov = curve_fit(fit_function, x_data0, y_data0, p0, sigma=None, absolute_sigma=True,maxfev = 100000)
    except:
        print ('No Gaussian fit possible (maxfev > 10000).')
        param = [-1, -1, -1, -1, -1, -1]
# Calculate chi-squared
    chisq = 0
    for i in range (0,len(x_data0)):   
        y_fit[i] = fit_function (x_data0[i], *param)
        chisq += (y_fit[i] - y_data0[i]) ** 2
    chisq = chisq / (len(x_data0)-1)

#     print("Fitting function coefficients:", param, "     Chi-squared: ", chisq)
#    print("Covariance of coefficients:")
#    print(param_cov)

    plt.plot(x_data0, y_data0, 'o', color ='red', label ="data")
    FWHM = 2.355 * round(param[2],1)
    ampl = param [0]
    peak_pos = param [1]
    bkg = d + e * peak_pos + f * peak_pos * peak_pos
    s_over_rootB = ampl/math.sqrt(bkg)
    FWHM = abs(round (FWHM, 1))

#Plot fitted function
    plt.plot(x_data0, fit_function(x_data0, *param), '--', color ='blue', label =f"Gaussian function:  FWHM = {FWHM} s")
    plt.legend()
#    plt.savefig ('Plot_Gaussian.pdf')
#    plt.show()

    return ampl, peak_pos, bkg, s_over_rootB, FWHM, chisq, d, e, f

###############################################################
# Gaussian fit function with quadratic background
def fit_function (x, a, b, c, d, e, f):
    return a * np.exp(-.5*((np.array(x)-b)/c)**2) + d + np.array(x) * e + np.square(x) * f 

###############################################################
def fitter_Norris(x_data, y_data, ref_MET, peak_pos1, FWHM, d0, e0, f0):
# Version of fitter for Norris function
#   Returns ampl, peak_pos, bkg, s_over_rootB, t90, chisq

    tmin0 = -4 * abs(FWHM)            # Perform Norris fit only over restricted range around Gaussian peak
    tmax0 = 4 * abs(FWHM)

# generate data over restricted range tmin0 to tmax0
    x_data0 = []
    y_data0 = []

    ij = 0
    peak_pos = peak_pos1
    peak = 0
    for i in range (0,len(x_data)):
        if x_data[i] > tmax0:
            break
        if x_data[i] >= tmin0 :
            x_data0.append (x_data[i])
            y_data0.append (np.sum(y_data, axis = 1)[i])
            bkgi = d0 + e0 *  x_data[i] + f0 * x_data[i] * x_data[i] 
            y_data0[ij] = y_data0[ij] - bkgi		# Subtract off background estimated for Gaussian
            if (abs(x_data[i]) - peak_pos) < 4 :
                bkg = bkgi
                peak = y_data0[ij]
            ij += 1

    tau1 = tau2 = FWHM  
    a = peak * np.exp(2)			# Differentiate Norris function, set to zero to determine amplitude
    b = peak_pos - FWHM   			# Initial estimate of starting point

#    params = [a, b, tau1, tau2] 

# curve_fit() function takes the test-function x-data and y-data as argument and returns 
# the coefficients a, b, tau1, tau2 in param and the estimated covariance of param in param_cov
 
    p0 = [ a, b, tau1, tau2, d0, e0, f0 ]	# Initial values of fit parameters
    try:
        param, param_cov = curve_fit(fit_Norris, x_data0, y_data0, p0, maxfev = 10000)
    except:
        print ('No Norris fit possible (maxfev > 10000).')    
        param = [-1,-1,-1,-1,-1,-1,-1]
    peak = param[0]
    b = param[1]
    param[2] = abs(param[2])
    param[3] = abs(param[3])
    tau1 = param[2]
    tau2 = param[3]
    e = param[4]
    f = param[5]
    g = param[6]

    y_Norris = fit_Norris (x_data0, peak, b, tau1, tau2, e, f, g)

# Calculate t90 from fit
    sumy1 = sumy = 0
    t0 = b
    t90 = 1000
    for i in range (0,len(x_data0)): 
        bkgi = e + f * x_data0[i] + g * x_data0[i] * x_data0[i]
        sumy += y_Norris[i] - bkgi			# Remove background, calculate total number of counts
    for i in range (0,len(x_data0)): 
        bkgi = e + f * x_data0[i] + g * x_data0[i] * x_data0[i]
        sumy1 += y_Norris[i] - bkgi
        if (sumy1 >= .05 * sumy) and (t0 == b) :
            t0 = x_data0[i]
        if sumy1 >= .95 * sumy:
            t90 = round (x_data0[i] - t0, 1)
            break
# Calculate chi-squared, maxy, peak_pos, bkg2
    chisq = 0   
    for i in range (0,len(x_data0)):   	
        bkgi = d0 + e0 * x_data0[i] + f0 * x_data0[i] * x_data0[i]
        y_Norris[i] += bkgi			# Add the Gaussian background back in 
        y_data0[i] += bkgi
        chisq += (y_Norris[i] - y_data0[i]) ** 2
#        print ('Norris   :', xdata0[i], ydata0[i])
    chisq = chisq / (len(x_data0)-1)

    maxy = peak * np.exp(-(2 * math.sqrt(tau2/tau1)))
    peak_pos = b + math.sqrt(tau1*tau2)
    bkg = e + d0 + (f + e0) * peak_pos + (g + f0) * peak_pos * peak_pos
    s_over_rootB = maxy/math.sqrt(max(bkg,0.01))
    plt.plot(x_data0, y_Norris, '+', color ='green', label =f"Norris function:  t90 = {t90} s")
    plt.legend(fontsize = '6',loc='upper left')
    plt.savefig ('Plot2.pdf')
#    plt.show()
#    plt.close()

    return maxy, peak_pos, bkg, s_over_rootB, t90, chisq

###############################################################
# Norris et al 2005 fit function with quadratic background
def fit_Norris (x, a, b, c, d, e, f, g):

    yNorris = e + np.array(x) * f + np.square(x) * g
    for i in range (0, len(x)):
        if x[i] > 1000:
            break
        if x[i] > b :
            yNorris[i] += a * np.exp(-(x[i]-b)/c - d/(x[i]-b))    

    return yNorris

#   return a * np.exp(-(np.array(x)-b)/c - d/(np.array(x)-b))  + e + np.array(x) * f + np.square(x) * g 

###############################################################
def get_th_data(data_dir, fileid):
    """
        This function loads six TH data files and merge them to one astropy Table.

        Args:
            data_dir (str): path to the data directory
            fileid (str): date of the data (yyyymmdd)

        Returns:
            Table: a table includes six TH data for each detector & gain
    """

    hx1_hth_file = '%s/cgbm_%s_hx1_hth.fits' % (data_dir, fileid)
    hx1_hth_table = Table.read(hx1_hth_file, hdu=1)
    hx2_hth_file = '%s/cgbm_%s_hx2_hth.fits' % (data_dir, fileid)
    hx2_hth_table = Table.read(hx2_hth_file, hdu=1)
    sgm_hth_file = '%s/cgbm_%s_sgm_hth.fits' % (data_dir, fileid)
    sgm_hth_table = Table.read(sgm_hth_file, hdu=1)
    hx1_lth_file = '%s/cgbm_%s_hx1_lth.fits' % (data_dir, fileid)
    hx1_lth_table = Table.read(hx1_lth_file, hdu=1)
    hx2_lth_file = '%s/cgbm_%s_hx2_lth.fits' % (data_dir, fileid)
    hx2_lth_table = Table.read(hx2_lth_file, hdu=1)
    sgm_lth_file = '%s/cgbm_%s_sgm_lth.fits' % (data_dir, fileid)
    sgm_lth_table = Table.read(sgm_lth_file, hdu=1)
    merge = Table()

    # Timestamps of CGBM data are end of the bins.
    # Timestamps were changed to center of the bins
    merge['TIME'] = hx1_hth_table['TIME'] - cg.TH_HBIN
    merge['MDCTIME'] = hx1_hth_table['MDCTIME'] - cg.TH_HBIN
    merge['HXM1_EXPOSURE'] = hx1_hth_table['EXPOSURE']
    merge['HXM1H_COUNTS'] = hx1_hth_table['COUNTS']
    merge['HXM1H_RATE'] = (
            hx1_hth_table['COUNTS']
            / np.array([
        hx1_hth_table['EXPOSURE'],
        hx1_hth_table['EXPOSURE'],
        hx1_hth_table['EXPOSURE'],
        hx1_hth_table['EXPOSURE']]).T)
    merge['HXM1L_COUNTS'] = hx1_lth_table['COUNTS']
    merge['HXM1L_RATE'] = (
            hx1_lth_table['COUNTS']
            / np.array([
        hx1_lth_table['EXPOSURE'],
        hx1_lth_table['EXPOSURE'],
        hx1_lth_table['EXPOSURE'],
        hx1_lth_table['EXPOSURE']]).T)

    merge['HXM2_EXPOSURE'] = hx2_hth_table['EXPOSURE']
    merge['HXM2H_COUNTS'] = hx2_hth_table['COUNTS']
    merge['HXM2H_RATE'] = (
            hx2_hth_table['COUNTS']
            / np.array([
        hx2_hth_table['EXPOSURE'],
        hx2_hth_table['EXPOSURE'],
        hx2_hth_table['EXPOSURE'],
        hx2_hth_table['EXPOSURE']]).T)
    merge['HXM2L_COUNTS'] = hx2_lth_table['COUNTS']
    merge['HXM2L_RATE'] = (
            hx1_lth_table['COUNTS']
            / np.array([
        hx2_lth_table['EXPOSURE'],
        hx2_lth_table['EXPOSURE'],
        hx2_lth_table['EXPOSURE'],
        hx2_lth_table['EXPOSURE']]).T)

    merge['SGM_EXPOSURE'] = sgm_hth_table['EXPOSURE']
    merge['SGMH_COUNTS'] = sgm_hth_table['COUNTS']
    merge['SGMH_RATE'] = (
            sgm_hth_table['COUNTS']
            / np.array([
        sgm_hth_table['EXPOSURE'],
        sgm_hth_table['EXPOSURE'],
        sgm_hth_table['EXPOSURE'],
        sgm_hth_table['EXPOSURE']]).T)
    merge['SGML_COUNTS'] = sgm_lth_table['COUNTS']
    merge['SGML_RATE'] = (
            sgm_lth_table['COUNTS']
            / np.array([
        sgm_lth_table['EXPOSURE'],
        sgm_lth_table['EXPOSURE'],
        sgm_lth_table['EXPOSURE'],
        sgm_lth_table['EXPOSURE']]).T)
    return merge

###############################################################
def plot_orb(data_dir, fileid, trigtime):
    """
        Routine copied from process_trigger to make an orbit plot and save as png file.

        Args:
            data_dir (str): path to the data directory
            fileid (str): date of the data (yyyymmdd)
            trigtime (float): trigger time in MET
        Returns:
            orb_table['LATITUDE'][diff_min_index], orb_table['LONGITUDE'][diff_min_index]

    """
    # mergin is a variable to set the time interval for the track of the pointing direction.
    mergin = 150.0
    orb_file = '%s/cgbm_%s.orb' % (data_dir, fileid)
    orb_table = Table.read(orb_file, hdu=1)
    orb_table['TIME'] = orb_table['TIME'] - cg.PERIODIC_HBIN
    neg_mask = orb_table['LONGITUDE'] > 180
    orb_table['LONGITUDE'][neg_mask] = orb_table['LONGITUDE'][neg_mask] - 360

    # Find the nearest time stamp
    orb_table['abs_diff'] = np.abs(orb_table['TIME'] - trigtime)
    diff_min_index = np.argmin(orb_table['abs_diff'])

    # Extracting row within triggertime +/- mergin
    around_mask = (trigtime - mergin <= orb_table['TIME']) & ( orb_table['TIME'] <= trigtime + mergin)
    around = orb_table[around_mask]

#    print ('Latitude = ',orb_table['LATITUDE'][diff_min_index],'   Longitude = ', orb_table['LONGITUDE'][diff_min_index])
    
    fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
    ax.set_global()
    ax.coastlines(resolution='110m', linewidth=1, edgecolor='black')
    plt.plot(
        around['LONGITUDE'],
        around['LATITUDE'],
        marker=',',
        ls='None')

    plt.plot(
        orb_table['LONGITUDE'][diff_min_index],
        orb_table['LATITUDE'][diff_min_index],
        color='purple',
        marker='*',
        ls='None'
    )

    ax.set_xticks(range(-180, 210, 30), crs=ccrs.PlateCarree())
    ax.set_yticks(range(-90, 120, 30), crs=ccrs.PlateCarree())
    ax.xaxis.set_major_formatter(plt.FixedFormatter(ax.get_xticks()))
    ax.yaxis.set_major_formatter(plt.FixedFormatter(ax.get_yticks()))
    plt.xlabel('Longitude [deg]')
    plt.ylabel('Latitude [deg]')
    plt.savefig ('Plot4.pdf')
#    plt.show()
    plt.close ('all')
    return orb_table['LATITUDE'][diff_min_index], orb_table['LONGITUDE'][diff_min_index]

###############################################################
def met2utc(trigtime_met):
    """
        This function get Time in astropy.time from trigger time in MET

        Args:
            trigtime_met (float): trigger time in MET
        Returns:
            utc (astropy.time.Time):
    """
    met = TimeDelta(trigtime_met, format="sec")
    utc = (met + cg.T0)
    return utc

###############################################################
def get_iss_att(trigtime_met, iat_table):
    """
        This function get ISS attitude data.

        Args:
            trigtime_met (float): trigger time in MET
            iat_table (astropy.table.Table): Table of orbit data
        Returns:
            iss2sky (scipy.spatial.transform.Rotation): Rotation matrix from ISS
            coordinate to the celestial coordinates.
    """

    iat_table['ZTIME'] = iat_table['TIME'] - trigtime_met

    # Find the nearest time stamp
    iat_table['ABS_TIME'] = np.abs(iat_table['TIME'] - trigtime_met)
    nearest_index = np.argmin(iat_table['ABS_TIME'])

    # Geting the quaternion
    qparam = iat_table['QPARAM'][nearest_index]
    iss2sky = R.from_quat(qparam)
#    Should I keep these printouts?
#    print('Input MET = %.6f' % (trigtime_met))
#    print('Nearest MET = %.6f' % (iat_table['TIME'][nearest_index]))
#    print('Time Diff. = %.6f' % (iat_table['TIME'][nearest_index]-trigtime_met))
#    print()
    return iss2sky

###############################################################
def get_cgbm_pointing(iss2sky, theta, phi, det):
    """
        This function get CGBM pointing directions.

        Args:
            iss2sky (scipy.spatial.transform.Rotation): Rotation matrix from ISS
            coordinate to the celestial coordinates.
            theta (float) : zenith angle [deg]
            phi (float) : azimuth [deg]
            det (str): detector name (HXM1/HXM2/SGM)
        Returns:
            sky_pos (astropy.coordinates.SkyCoord): Sky position of the direction
            which was represented by theta and phi in the detector coordinate.
    """

    det = det.upper()

    if theta < 0 or theta > 180:
        raise ValueError("theta must be in the [0, 180] deg range.")
    if phi < 0 or phi > 360:
        raise ValueError("phi must be in the [0, 360] deg range")

    # Convert position from spherical to cartesian coordinates.
    lon = phi * u.deg
    lat = (90.0 - theta) * u.deg
    pos_cgbm = UnitSphericalRepresentation(lon, lat).to_cartesian()

    # Rotation from the CGBM to the celestial coordinates
    cgbm2sky = iss2sky * cg.CAL2ISS * cg.CGBM2CAL[det]

    # Rotate the position
    sky_pos = cgbm2sky.apply(pos_cgbm.xyz)

    # Transfrom to SkyCoord with ra, dec
    sky_pos = SkyCoord(sky_pos[0], sky_pos[1], sky_pos[2],
                       representation_type='cartesian', frame='icrs')
    sky_pos.representation_type = 'spherical'
#    print('%s zenith (R.A., Dec.) = (%.3f, %.3f)' % (det , sky_pos.ra.deg, sky_pos.dec.deg))
    return sky_pos

###############################################################
def get_cgbm_theta_phi(iss2sky, ra, dec, det):
    """
        This function get the incident angle of the direction which is given as ra and dec.

        Args:
            iss2sky (scipy.spatial.transform.Rotation): Rotation matrix from ISS
            coordinate to the celestial coordinates.
            ra (float) : Right Ascension [deg]
            dec (float) : Declination [deg]
            det (str): detector name (HXM1/HXM2/SGM)
        Returns:
            pos_cgbm_deg ((float, float)): Incident angles (zenith angle [deg], azimuth [deg])
    """

    sky2iss = iss2sky.inv()
    det = det.upper()

    pos_sky = SkyCoord(ra, dec, unit='deg', frame='icrs')

    # Rotation from sky to CGBM coordinates

    sky2cgbm = cg.CAL2CGBM[det] * cg.ISS2CAL * sky2iss

    # `pos_sky.cartesian.xyz` yields a Quantity object, adding `.value`
    # gives `numpy.ndarray` (though this will also work with Quantity)
    pos_cgbm = sky2cgbm.apply(pos_sky.cartesian.xyz.value)
    theta = np.arccos(pos_cgbm[2])
    phi = np.arctan2(pos_cgbm[1], pos_cgbm[0])

    # np.arctan2 yields an angle in the [-pi, pi] range, so we need to
    # add 2*pi for negative values to get the [0, 2*pi] range

    if phi < 0:
        phi += 2 * np.pi
    pos_cgbm_deg = (np.rad2deg(theta), np.rad2deg(phi))
    return pos_cgbm_deg

###############################################################
def get_solar_incident_angle(iss2sky, trigtime_met, det):
    """
        This function get the incident angle of the direction which is given as ra and dec.

        Args:
            iss2sky (scipy.spatial.transform.Rotation): Rotation matrix from ISS
            coordinate to the celestial coordinates
            trigtime_met (float): trigger time in MET
            det (str): detector name (HXM1/HXM2/SGM)
        Returns:
            pos_cgbm_deg ((float, float)): Incident angles (zenith angle [deg], azimuth [deg]) of the sun
            sun_ra_rad, sun_dec_rad
    """

    trigtime_utc = met2utc(trigtime_met)

    # Get solar position
    sun = get_sun(trigtime_utc)
    solar_angle = get_cgbm_theta_phi(iss2sky, sun.ra.deg, sun.dec.deg, det)

#    print('%s solar angle (theta, phi) = (%.3f, %.3f)' % (det , solar_angle[0], solar_angle[1]))
#    print ('get_solar_incident_angle: ',sun.ra.deg, sun.dec.deg)

    sun_ra_rad = sun.ra.deg * np.pi/180.
    sun_dec_rad = sun.dec.deg * np.pi/180.
#    print ('get_solar_incident_angle in radians: ',sun_ra_rad, sun_dec_rad)

    return solar_angle, sun_ra_rad, sun_dec_rad

###############################################################
def get_angles(data_dir, fileid, trigtime_met):
    """
        This function calculates pointing directions and solar angles and print them.

        Args:
            data_dir (str): path to the data directory
            fileid (str): date of the data (yyyymmdd)
            trigtime_met (float): trigger time in MET
        Returns:
            ra_rad, sun_dec_rad, zenith_hxm.ra.deg, zenith_hxm.dec.deg, zenith_sgm.ra.deg, zenith_sgm.dec.deg, solar_hxm, solar_sgm
    """

    iat_file = '%s/cgbm_%s.iat' % (data_dir, fileid)
    iat_table = Table.read(iat_file, hdu=1)
    iat_table['TIME'] = iat_table['TIME'] - cg.PERIODIC_HBIN

    iss2sky = get_iss_att(trigtime_met, iat_table)
    zenith_hxm = get_cgbm_pointing(iss2sky, 0., 0., 'HXM')
    zenith_sgm = get_cgbm_pointing(iss2sky, 0., 0., 'SGM')

    solar_hxm, sun_ra_rad, sun_dec_rad = get_solar_incident_angle(iss2sky, trigtime_met, 'HXM')
    solar_sgm, sun_ra_rad, sun_dec_rad = get_solar_incident_angle(iss2sky, trigtime_met, 'SGM')    

    return sun_ra_rad, sun_dec_rad, zenith_hxm.ra.deg, zenith_hxm.dec.deg, zenith_sgm.ra.deg, zenith_sgm.dec.deg, solar_hxm, solar_sgm

###############################################################
def process_trigger(data_dir, fileid, trigid, trigtime_met, trigtime_mdc):
    """
        This function make plots and gets pointing directions & solar angles.

        Args:
            data_dir (str): path to the data directory
            fileid (str): date of the data (yyyymmdd)
            trigid (int): integer of the trigger time in MDC time
            trigtime_met (float): trigger time in MET
            trigtime_mdc (float): trigger time in MDC time
        Returns:
            sun_ra_rad, sun_dec_rad, hxm_ra_deg, hxm_dec_deg, sgm_ra_deg, sgm_dec_deg, solar_hxm, solar_sgm, lat,lon, UT

    """

    th_data = get_th_data(data_dir, fileid)

    # Create CGBM light curves
    # plt_lc_single_det(th_data, trigid, trigtime_met, 'HXM1')
    # plt_lc_single_det(th_data, trigid, trigtime_met, 'HXM2')
    # plt_lc_single_det(th_data, trigid, trigtime_met, 'SGM')

    # Get the ISS orbit
    lat, lon = plot_orb(data_dir, fileid, trigtime_met)
    trigtime_utc = met2utc(trigtime_met)

    # Printing the triggerID and trigger time
    # print('TriggerID: %d' % (trigid))
    # print('Trigger Time (MDC): %.6f' % (trigtime_mdc))
    # print('Trigger Time (MET): %.6f' % (trigtime_met))
    # print('Trigger Time (UT): %s' % (trigtime_utc.isot))

    UT = trigtime_utc.isot
    # Printing the pointing directions and solar angles
    sun_ra_rad, sun_dec_rad, hxm_ra_deg, hxm_dec_deg, sgm_ra_deg, sgm_dec_deg, solar_hxm, solar_sgm = get_angles(data_dir, fileid, trigtime_met)
    return sun_ra_rad, sun_dec_rad, hxm_ra_deg, hxm_dec_deg, sgm_ra_deg, sgm_dec_deg, solar_hxm, solar_sgm, lat,lon, UT

###############################################################
def Haversine(delta1, delta2, alpha1, alpha2):
# routine to calculate difference between two locations on sky
    arg1 = math.sin( .5*(delta2 - delta1))
    arg2 = math.sin( .5*(alpha2 - alpha1))
    arg = arg1 * arg1 + math.cos(delta1) * math.cos(delta2) *arg2 * arg2
    theta = 2. * math.asin( math.sqrt(arg))
# Astronomy cafe
#    theta = math.acos ( math.sin(delta1) * math.sin(delta2) + math.cos(delta1) * math.cos(delta2) * math.cos(alpha1 - alpha2) )
    return theta

###############################################################
def solar (fileid, time_UT, duration):
# Solar activity reader from Henry Vo, 2/28/2024

    # pip install bs4 
    # pip install selenium 
    # pip install pandas
    # pip install requests

    from selenium import webdriver
    from selenium.webdriver.common.by import By
    import requests
    from PIL import Image
    from io import BytesIO
    from datetime import datetime, time, timedelta

    # Create a webdriver object to automate chrome
    driver = webdriver.Chrome()
    # Open the solar activity webpage
    driver.get('https://www.polarlicht-vorhersage.de/goes-archive/')

    # Prompt the user for date, time, and duration of the graph
    # date = input('Enter the date of the desired event (yyyymmdd)')
    # fileid = '20231201'
    date = fileid
    # time_UT = input('Enter the time of the desire event (UT) (HHMM)')
    time_UT = "1200"
    # duration = input("Enter the desired time duration")
    duration = "12"

### Filter the graph to be in the desired time window
#### First, change the date to time window
# Find all the input dropboxes from the html
    elements = driver.find_elements(By.XPATH, '//input')

    min_date_element = elements[0] # First date dropbox html element
    min_time_element = elements[1] # First time dropbox html element
    max_date_element = elements[2] # Second date dropbox html element
    max_time_element = elements[3] # Second time dropbox html element

### Take the inputted date and convert to a usable date and time that can be subtracted and added by the time duration
    date = datetime.strptime(date, "%Y%m%d") 

    hour = time_UT[:2]
    minute = time_UT[2:]

    time_UT = time(int(hour), int(minute), 0)
    datetime_obj = datetime.combine(date, time_UT) #Create a datetime object with the date and time

    max_datetime = datetime_obj+timedelta(hours=int(duration)) # Add time duration +- to the datetime object
    min_datetime = datetime_obj-timedelta(hours=int(duration))

    max_date = max_datetime.strftime("%m/%d/%Y")
    min_date = min_datetime.strftime("%m/%d/%Y")
    max_time = max_datetime.strftime("%I:%M %p")
    min_time = min_datetime.strftime("%I:%M %p")

# Clear all the dropboxes of any information
    min_date_element.clear() 
    min_time_element.clear()
    max_date_element.clear()
    max_time_element.clear()

# Insert the desired values into the dropboxes
    min_date_element.send_keys(min_date)
    min_time_element.send_keys(min_time)
    max_date_element.send_keys(max_date)
    max_time_element.send_keys(max_time)

# Print the image
    graph_element = driver.find_elements(By.XPATH, '//img')[2]

# Take the link from the 'src' attribute in the html
    graph = graph_element.get_attribute('src')

# Open the link and display it as a .png
    response = requests.get(graph)
    image_date = response.content
    image = Image.open(BytesIO(image_date))
#    image.show()
    image.save('Plot_solar.pdf')
    image.close ()

# tmp0oy0f1dj.png

# Close out of the webpage
    driver.quit()

    return

###############################################################
# X-Ray Flux Data Retriever.py -- Henry Vo,5/21/24
# Import the libraries

def Xray_flux1 (fileid, UT):
    import os
    import requests
    import xarray as xr
    from bs4 import BeautifulSoup
    import netCDF4 as nc
    import datetime
    import numpy as np
    import sys

    print ('*** 5a Enter X-ray')
# Take the input date and sort into year, month, day

# Input the trigger date
#    trig_date = input("Enter the date with a format: yyyymmdd") + "T00:00:00Z"
    trig_date = fileid + "T00:00:00Z"

# Split the date into year, month, and day by indexing the date string
    year = trig_date[:4]
    month = trig_date[4:6]
    day = trig_date[6:8]

# Determine current directory
#    print(os.getcwd() + "\n")
    current_dir = os.getcwd()

# Determine which GOES contains the data
    if trig_date[:4] == '2022':
        if month in ('09,10,11,12'):
            goes = 'goes18'
        else:
            goes = 'goes17'
    elif trig_date[:4] in ('2024,2023'):
        goes = 'goes18'
    elif trig_date[:4] in ('2018,2019,2020,2021'): 
        goes = 'goes17'
    else:
        goes = 'goes16'

# Find the link of the year and the date suggested
    base_link = 'https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/' + goes + '/l2/data/xrsf-l2-flx1s_science/' + year +'/' + month + '/'

# Get the HTML content of the webpage
    response = requests.get(base_link)

# Parse the HTML content using BeautifulSoup
    soup = BeautifulSoup(response.content, 'html.parser')

# ## Look  through the <\a> tags in the  html and look for the link of the target date, then download the nc file

    a_links = soup.find_all('a')

# Loop through each <a> tag ang find the link to download the nc file
    target_link = ""
    for link in a_links[4:-1]:
        if len (link['href'].split('_')) > 3:
            if day == link['href'].split('_')[3][7:9]:
#                print ('*** day = ')
                target_file = link.get('href')
                target_link = base_link + link.get('href')
    if target_link == "":
        Average = Variance = GOES_avg_1hr = GOES_var_1hr = -1
        return Average, Variance, GOES_avg_1hr, GOES_var_1hr

#    print ('*** 5b, link   ')
#    print (' 5b, link again  ', target_link)

# Get the content from the target_link
    response = requests.get(target_link)
#    print ('*** 5b0')
# ## Write the nc file to a local directory named nc_files

# #### First check to see if a nc_files directory already exists, if not make one

# Set the base of the path where the nc files should be stored
#base = r'C:\Users\henry\OneDrive - Louisiana State University\Laptop Desktop\Research\Jupyter'
#    print ('*** 5b1')
    base = r'C:\Users\Owner\code\classifier\classify\Xray'
#    print ('*** 5b2')

# Join the base of the path with the file name 'nc_files'
    flux_data_path = os.path.join(base, 'nc_files')

# Make a nc_files directory if one doesnt exist. If it already exists, change the working directory to the nc_files directory
    if not os.path.exists(flux_data_path):
      os.makedirs(flux_data_path)
    else:
        os.chdir(flux_data_path)
#    print ('*** 5b3')

#### Now write the nc_file to this directory

# Generate a unique filename based on the current date and time
    timestamp = str(year +  month + day)
    nc_file_path = f"flux_data_{timestamp}.nc" # Set the name of the file as data_yyyymmdd 


    cwd_files = os.listdir('.') # Get all the files inside the nc_files directory
# If the file does not already exist in the nc_files folder, write an nc.file

    if nc_file_path not in cwd_files:
        with open(nc_file_path, "wb") as f:
            f.write(response.content)

# Open the NetCDF file
    nc_file = nc.Dataset(nc_file_path)

##     *****************************                Calculate the average counts and variance               ************************
#    print ('*** 5b4')

    flux = nc_file.variables['xrsa_flux'][:]
    Average = np.average(flux)
    Variance = np.sqrt(np.var(flux) / (86400-1))

    print ("UT:  ",UT, UT[11:13], UT[14:16])
    trigger_time_in_secs = int(UT[11:13]) * 3600 + int(UT[14:16]) * 60
    print (trigger_time_in_secs)
    start_time = trigger_time_in_secs - 3600
    end_time = trigger_time_in_secs + 3600
    GOES_avg_1hr = np.average(flux[start_time:end_time])
    GOES_var_1hr = np.sqrt(np.var(flux[start_time:end_time]) / (7200-1))
    print ("Fluxes:  ", Average, Variance, GOES_avg_1hr, GOES_var_1hr)

#    print("Average_Flux: ", Average, "      Variance: ", Variance)

## Now plot the data from the newly written nc_file 

# Change directory into the one where XFlux_plotter.py is located
#os.chdir(r'C:\Users\henry\OneDrive - Louisiana State University\Laptop Desktop\Research\Jupyter')
    os.chdir(r'C:\Users\Owner\code\classifier\classify\Xray')
    print ('*** 5c')


# Import the plot_values function from XFlux_plotter.py
# from XFlux_plotter import plot_xflux

# Call the function to plot the data
# plot_xflux(nc_file)

    """
XFlux_plotter.py

This function takes one argument:
1) A nc_file object from the netCDF4 library

The function takes the values from the nc file and plots it over a 24 hr window. 
    """
# def plot_xflux(nc_file):
#    print ('*** 5d')
    import netCDF4 as nc
    import numpy as np
    import matplotlib.pyplot as plt
    from datetime import datetime, timedelta

    # Extract which goes from the ncfile
    goes = nc_file.platform

    # Set the y variable to xrayflux
    flux = nc_file.variables['xrsa_flux'][::5]

    # Set the x variable to time
    start_time = datetime.strptime(nc_file.time_coverage_start, '%Y-%m-%dT%H:%M:%S.%fZ')

    # Create an array of 17280 datetime elements, starting from 0 to 86400-1 that represents every 5 seconds
    array_of_seconds = np.arange(0.0, 86400.0-1.0, 5)                                      

    time_index = [start_time + timedelta(seconds=i) for i in array_of_seconds]    

    # Plot the xrayflux(y) vs. time_index(x)
    plt.figure(figsize=(10,8)) # Set the size of the graph
    plt.plot(time_index, flux)  

    plt.yscale('log')
#    plt.title(goes.upper() + ' X-ray Flux vs. Time   (Trigger time: ' + UT + ')') # Set the title of the plot
    plt.title('GOES X-ray flux vs. time        Trigger time (UT): ' + UT)
    plt.ylim(1e-8, None)
    plt.yticks([1e-9, 1e-8, 1e-7,1e-6]) # Set the y axis tick markers

# Switch back to original directory to store plot fileand continue program execution
    os.chdir(current_dir)
    plt.savefig ("GOES_flux.pdf")

#    plt.show()
    plt.close ('all')
    if GOES_avg_1hr > 0 and GOES_var_1hr > 0:
        print ('***5e GOES values: ', Average, Variance, GOES_avg_1hr, GOES_var_1hr)
    else:
#        Average = Variance = GOES_avg_1hr = GOES_var_1hr = -1
        GOES_avg_1hr = GOES_var_1hr = -1
        print ('***5f GOES values: ', Average, Variance, GOES_avg_1hr, GOES_var_1hr)

    if Average == 'nan':		# Check that GOES values are not 'nan'
        Average = -1
    if Variance == 'nan':
        Variance = -1
    if GOES_avg_1hr == 'nan':
        GOES_avg_1hr = -1
    if GOES_var_1hr == 'nan':
        GOES_var_1hr = -1
 
    return Average, Variance, GOES_avg_1hr, GOES_var_1hr

###############################################################
def analyze_single_trigger (fileid, chan_type, data_dir, ph_table, trig_index, trig_time, trig_num, triggers,
event_id_rows, event_UT_rows, event_type_rows, tmin, tmax):
#    This is basically analyze8.py
    trigtime_met = ref_MET = trig_time[trig_num]  #  Multiple triggers?   # These three lines ensure compatibility with process_trigger
    trigtime_mdc = trig_index[trig_num]
    trigid = int(trigtime_mdc)

    ph_table['ZTIME'] = ph_table['TIME'] - cg.PH_HBIN - ref_MET

################################
# Find nearest hv turn on/off time
# Revised version -- Search over +-10800 s = 2700 time bins -- 6/17/24
    """
    threshold = 2000
    ir_on = ir_off = transition = turn_off = turn_on = -86400

    ir_end = len(ph_table['ZTIME'])
    for ir in range (1,ir_end):
        if ph_table['ZTIME'][ir] > -4000:			# Start search for transition at 4000 before trigger
            ir_start = ir + 1
            break
    for ir in range (ir_start,min(ir_start+2000,ir_end)):
        print (i2, ph_table['ZTIME'][ir])
#        xx1 = int(np.sum(ph_table['HXM1L_COUNTS'], axis=1)[ir])
#        xx2 = int(np.sum(ph_table['HXM1L_COUNTS'], axis=1)[ir-1])
        xx1 = int(np.sum(ph_table['SGML_COUNTS'], axis=1)[ir])
        xx2 = int(np.sum(ph_table['SGML_COUNTS'], axis=1)[ir-1])

        if (xx1-xx2) > threshold:				# Starting at ir = -2000, is difference between two successive records > 2000?
            turn_on = ph_table['ZTIME'][ir]
            ir_on = ir						# Identify index ir at which turn_on occurs
        if (xx2-xx1) > threshold:
            turn_off = ph_table['ZTIME'][ir]
            ir_off = ir						# Identify index ir at which turn_off occurs
        if turn_on > -86400 and turn_off > -86400 :		# Have I found both an up and a down transition within 2000 counts of trigger?
            break
    if abs(turn_on) > abs(turn_off):
        transition = int(turn_off)				# Transition = time of turn on/off, ir_off/ir_on = index
        ir = ir_off						#    - ==> turn on, + ==> turn off
    else:							#    -86400 means turn on at very start of event or never found a transition
        transition = int(turn_on)
        ir = ir_on
    """
    threshold = 1000
    imax= len(ph_table['ZTIME'])
#    print ('Start time: ', ph_table['ZTIME'][0], '  End time: ',ph_table['ZTIME'] [imax-1] )
    transition = 86400
    ir_transition = 0
    ir_start = 0
    ir_end = len(ph_table['ZTIME'])
#    for ir in range (ir_start, ir_end,5):
#        print ('@@@@',ir, ph_table['ZTIME'][ir], np.sum(ph_table['HXM1L_COUNTS'], axis=1)[ir])
    ir_start = -float(ph_table['ZTIME'][0])/4 - 2700		# Start 1/2 orbit prior to trigger, extend to 1/2 orbit after
    ir_start = int(ir_start)
    ir_end = ir_start + 5400
    for ir in range (ir_start, ir_end):
        xx1 = int(np.sum(ph_table['HXM1L_COUNTS'], axis=1)[ir])
        xx2 = int(np.sum(ph_table['HXM1L_COUNTS'], axis=1)[ir-1])
        if abs(xx1-xx2) > threshold:				# Starting at ir = -2700, is difference between two successive records > 2000?
            turn_on = int(ph_table['ZTIME'][ir])		# Tentative time at which transition occurs
            if abs(turn_on) < abs(transition):
                transition = turn_on				# Identify time at which closest turn_on/off occurs
                ir_transition = ir				# Identify index at which closest turn_on/off occurs

    print ('Time to turn on/off (sec) = ',transition)

#   Print counts vs time around hv transition
    ir_first = ir_transition - 10
    if (ir_first < ir_start):
        ir_first = ir_start
    ir_last = ir_first + 20
    if (ir_last > ir_end):
        ir_last = ir_end
    if ir_first <= ir_last :
        print ('Counts around hv turn on/off transition:')
        for i in range (ir_first, ir_last):
            print(' HXM1L counts near hv transition: ', ph_table['ZTIME'][i],np.sum(ph_table['HXM1L_COUNTS'], axis = 1)[i])
        print ()
#    print ('*** Finished 1')

################################
# Plot distribution of counts vs time 
    time_dist(ph_table['ZTIME'], ph_table['HXM1L_COUNTS'], ph_table['HXM2L_COUNTS'], ph_table['SGML_COUNTS'], ref_MET, fileid, trig_num)
   
################################
# Fit counts vs time plus background -- Gaussian
    ampl1, peak_pos1, bkg1, s_over_rootB1, FWHM, chisq1, d, e, f = fitter(ph_table['ZTIME'], ph_table['SGML_COUNTS'], ref_MET, tmin, tmax)
# Fit counts vs time plus background -- Norris function

    ampl2, peak_pos2, bkg2, s_over_rootB2, t90, chisq2 = fitter_Norris(ph_table['ZTIME'], 
        ph_table['SGML_COUNTS'], ref_MET, peak_pos1, FWHM, d, e, f)
#    print ('***** 2 Return from Norris fit')
    
################################
# Plot spectra vs energy
#   Use default start and stop time

#    tmin = float(input("Please input the start time:"))
#    tmax = float(input("Please input the end time:"))
#    print ()
#    tmin = -200
#    tmax = 200
    tmask = (tmin <= ph_table['ZTIME']) & (ph_table['ZTIME'] <= tmax)
    roi = ph_table[tmask]
    row_num = float(len(roi))

    detectors = ['HXM1H', 'HXM2H', 'SGMH', 'HXM1L', 'HXM2L', 'SGML']
    id = 0  
    sumplot = np.zeros (len(detectors))
    for detector in detectors:
        y = np.sum(roi[detector + '_RATE'], axis = 0)
        if (detector == 'HXM1H' or detector == 'HXM2H' or detector == 'SGMH'):      # Change x, z values for first 3 high-threshold detectors
            Emin = 20
            Emax = 100
            x = np.arange(cg.PH_BINNUM_HIGH)
            fig = plt.figure(figsize=(8, 6))
            plt.plot( x, y / cg.high_ph_width() / row_num, drawstyle='steps-mid')
        else:                                                                       # Change x, z values for last 3 low-threshold detectors
            Emin = 10
            Emax = 400
            x = np.arange(cg.PH_BINNUM_LOW)
            fig = plt.figure(figsize=(8, 6))
            plt.plot( x, y / cg.low_ph_width() / row_num, drawstyle='steps-mid')

#    Calculate hardness ratio between Emin and Emax separately for low- and high-threshold detectors
        for i in range (0, len(x)):        # Sum up counts in appropriate energy ranges
            if (x[i] > Emin) & (x[i] < Emax):
                sumplot[id] += y[i]
   
        id += 1
        plt.yscale('log')
        plt.xlabel('PH channel')
        plt.ylabel('Counts /s /ADC channel')
        plt.title(detector + ' spectrum')
        fig.savefig ('plot3'+str(id)+'.pdf')
#       plt.show()
        plt.close (fig)

#    print ('**** 3 Return from plotting spectra vs energy')

################################
#     Hardness ratios
    HXM1_ratio = sumplot[3]/sumplot[0]    
    HXM2_ratio = sumplot[4]/sumplot[1]
    SGM_ratio  = sumplot[5]/sumplot[2]
#    print('Hardness ratios --      HXM1: %5.1f     HXM2: %5.1f     SGM: %5.1f \n' % (HXM1_ratio, HXM2_ratio, SGM_ratio))
#    print ('****   4  Return from hardness ratio calculation')

################################
#     determine latitude
#    plot_orb(data_dir, fileid, ref_MET)     # Multiple trigger times?

################################
#    determine latitude and angle to Sun
    sun_ra_rad, sun_dec_rad, hxm_ra_deg, hxm_dec_deg, sgm_ra_deg, sgm_dec_deg, solar_hxm, solar_sgm, lat, lon, UT = process_trigger(
          data_dir, fileid, trigid, trigtime_met, trigtime_mdc)
    if (sun_dec_rad) < 0 :
        sun_dec_rad += 2.*math.pi

    hxm_ra_rad = hxm_ra_deg * math.pi/180.
    hxm_dec_rad = hxm_dec_deg * math.pi/180.
    sgm_ra_rad = sgm_ra_deg * math.pi/180.
    sgm_dec_rad = sgm_dec_deg * math.pi/180.

    alpha1 = sun_ra_rad
    alpha2 = hxm_ra_rad
    alpha3 = sgm_ra_rad
    dec1 = sun_dec_rad
    dec2 = hxm_dec_rad
    dec3 = sgm_dec_rad

# Calculation between 2 points with RA, dec using Haversine formula
    theta12 = Haversine (dec1, dec2, alpha1, alpha2) * 180./math.pi
    theta13 = Haversine (dec1, dec3, alpha1, alpha3) * 180./math.pi

# Calculation between 2 points with theta, phi using dot product of vectors
    theta1 = phi1 = 0.
    theta2 = solar_hxm[0] * math.pi/180.
    phi2 = solar_hxm[1] * math.pi/180.
    delta_theta_hxm = math.acos( math.sin(theta1) * math.cos(phi1) * math.sin(theta2) * math.cos(phi2) 
                          + math.sin(theta1) * math.sin(phi1) * math.sin(theta2) * math.sin(phi2)
                          + math.cos(theta1) * math.cos(theta2) ) * 180./math.pi
    theta3 = solar_sgm[0] * math.pi/180.
    phi3 = solar_sgm[1] * math.pi/180.
    delta_theta_sgm = math.acos( math.sin(theta1) * math.cos(phi1) * math.sin(theta3) * math.cos(phi3) 
                          + math.sin(theta1) * math.sin(phi1) * math.sin(theta3) * math.sin(phi3)
                          + math.cos(theta1) * math.cos(theta3) ) * 180./math.pi

################################
#     Solar activity

    global GOES_avg, GOES_var
    duration = "1200"
    GOES_avg = GOES_var = -1
#    print('*** 5 Enter solar', trig_num, GOES_avg, GOES_var)
    if trig_num == 0 or trig_num > 0:
        GOES_avg, GOES_var, GOES_avg_1hr, GOES_var_1hr = Xray_flux1 (fileid, UT)
        try: 
            solar (fileid, UT, duration)
        except:
            duration = '1200'
    print ('**** 6 Trig_num =', trig_num, ' Return from solar', GOES_avg, GOES_var, GOES_avg_1hr, GOES_var_1hr)
################################

    """
# As of 2/22/2004 -- Calculate likelihoods for -- old version
#        lat, transition, angle_to_Sun = .5 *(solar_hxm[0] + solar_sgm [0]), FWHM, t90
# To add categories, add them to data_types list in main routing and likelihood_list here
# data_types, x, Ngtx, event_classes are global variables
# data_types = ['Latitude','Time_to_transition','angle_to_Sun','FWHM','t90']	-- Copy of data_types from main routine
# event_classes = ['GRB', 'Solar','Particle']				-- Copy of event classes from main routine
# Ngtx [itype][events][row_number] = int (Ngtx [itype][events][row_number]) / norm [events]  --  Copy from main routine

    angle_to_Sun = .5 *(solar_hxm[0] + solar_sgm [0])
    likelihood_list = [lat, transition, angle_to_Sun, FWHM, ampl1, chisq1, t90, ampl2, chisq2, SGM_ratio, GOES_avg, GOES_var, 
        GOES_avg_1hr, GOES_var_1hr]

    LGRB = [0] * len(data_types)
    Lsolar =[0] * len(data_types)
    Lpart = [0] * len(data_types)
    for idata in range (len(data_types)):
        for ievent in range (len(event_classes)):
            total = 0
            if ievent == 0:
                LGRB [idata] = 0 
            if ievent == 1:
                Lsolar [idata] = 0
            if ievent == 2:
                Lpart [idata] = 0
            jrow = 0
            while (jrow < max_rows):
                jrow += 1
#                print ('lat = ', lat, likelihood_list [idata], type(likelihood_list [idata]))
#                print ('Ngtx ', idata, ievent, jrow, Ngtx [idata][ievent][jrow], type (Ngtx [idata][ievent][jrow]))  
#                print ('xL ', idata, ievent, jrow, xL [idata][ievent][jrow], type (xL [idata][ievent][jrow]))
#                print ('@@@@@$$$$$  ', idata, ievent, jrow, xL [idata][ievent][jrow], Ngtx [idata][ievent][jrow],
#                    type(xL [idata][ievent][jrow]), type(Ngtx [idata][ievent][jrow]) )
                if likelihood_list [idata] < float(xL [idata][ievent][jrow]) or Ngtx [idata][ievent][jrow] == 0:
                    jrow -= 1
                    if ievent == 0:
                        LGRB [idata] = Ngtx [idata][ievent][jrow] 
                    if ievent == 1:
                        Lsolar [idata] = Ngtx [idata][ievent][jrow] 
                    if ievent == 2:
                        Lpart [idata] = Ngtx [idata][ievent][jrow]  
                    total += Ngtx [idata][ievent][jrow]  
                    jrow = 10000 
#            continue
#        print ('***** 9c ', idata, ievent)
        LGRB [idata] = LGRB [idata] / total
        Lsolar [idata] = Lsolar [idata] / total
        Lpart [idata] = Lpart [idata] / total    

    sum_data_types = [0] * len(event_classes)						# Calculate Sum (L) and Product (ln L)
    sum_log_data_types = [1] * len(event_classes)
    count_sum = [0] * len(event_classes)
    for idata in range (len(data_types)):
        sum_data_types [1] += Lsolar [idata]
        sum_data_types [2] += Lpart [idata]
        if LGRB [idata] > 0:
            sum_data_types [0] += LGRB [idata]
            sum_log_data_types [0] = sum_log_data_types [0] + math.log (LGRB[idata])
            count_sum [0] += 1
         if Lsolar [idata] > 0:
            sum_data_types [1] += Lsolar [idata]
            sum_log_data_types [1] = sum_log_data_types [1] + math.log (Lsolar[idata])
            count_sum [1] += 1
         if Lpart [idata] > 0:
            sum_data_types [2] += Lpart [idata]
            sum_log_data_types [2] = sum_log_data_types [2] + math.log (Lpart[idata])
            count_sum [02 += 1
     
    for ievent in range (len(event_classes)):						# Normalize by number of non-zero likelihoods
#        sum_data_types [ievent] = sum_data_types [ievent] / len (data_types)              
        sum_data_types [ievent] = sum_data_types [ievent] / count_sum [ievent]
        sum_log_data_types [ievent] = - math.exp (sum_log_data_types [ievent] ** (1 / count_sum [ievent] ))

#    print ('*** 10  End likelihoods, enter print section, trig_num = ', trig_num)
    """

##################################################################    

# As of 6/5/2024 -- Calculate likelihoods for -- New version
#        lat, transition, angle_to_Sun = .5 *(solar_hxm[0] + solar_sgm [0]), FWHM, t90
# To add categories, add them to data_types list in main routing and likelihood_list here
# data_types, x, Ngtx, event_classes are global variables
# data_types = ['Latitude','Time_to_transition','angle_to_Sun','FWHM','t90']	-- Copy of data_types from main routine
# event_classes = ['GRB', 'Solar','Particle']				-- Copy of event classes from main routine

    angle_to_Sun = .5 *(solar_hxm[0] + solar_sgm [0])
    likelihood_list = [lat, transition, angle_to_Sun, FWHM, ampl1, chisq1, t90, ampl2, chisq2, SGM_ratio, GOES_avg, GOES_var,
        GOES_avg_1hr, GOES_var_1hr]
    L_ratio = [0] * len(event_classes)

    LGRB = [0] * len(data_types)
    Lsolar =[0] * len(data_types)
    Lpart = [0] * len(data_types)
#    for idata in range (len(data_types)):
    for idata in range (len(data_types)):
        x = float( likelihood_list [idata])
        if idata == 2: 				# Normalizations
            arg = math.sin (x *3.14159265/180.)
            x = arg * arg
        if idata == 10:				
            x = x * 1.0e8
        if idata == 11:				
            x = x * 1.0e12
        
        total = 0
        for iclass in range (len(event_classes)):
#            for bin in range (0,10):
#                print (likelihood_list [idata], midpoints [idata] [bin] )
#                print (likelihood [idata] [iclass] [bin])
            for bin in range (1, max_bins):
                L_ratio [iclass] = 0
#                if likelihood_list [idata] > float(midpoints [idata] [max_bins - 1]):
                if x > float(midpoints [idata] [max_bins - 1]):
                    break
#                if likelihood_list [idata] <= float(midpoints [idata] [bin]):
                if x <= float(midpoints [idata] [bin]):
#                    x = float( likelihood_list [idata])
                    x2 =float( midpoints [idata] [bin])
                    x1 = float( midpoints [idata] [bin-1])
                    if x < x1:
                        break
                    y2 = float( likelihood [idata] [iclass] [bin] )   
                    y1 = float( likelihood [idata] [iclass] [bin-1])  
                    L_ratio [iclass] = (y2 - y1) / (x2 - x1) * (x - x1) + y1
#                    if idata == 2:
#                        print (x, x1, x, x2, y1, y2, L_ratio [iclass])
                    break
            total += L_ratio [iclass]
        for iclass in range (len(event_classes)):
            if total != 0:
                L_ratio [iclass] = L_ratio [iclass] / total    # For each event, L_ratio is GRB/solar/particle likelihood for particular data class
#        print (total, L_ratio [0], L_ratio [1], L_ratio [2])
        LGRB [idata] = L_ratio [0] 
        Lsolar [idata] = L_ratio [1]
        Lpart [idata] = L_ratio [2] 
#        sys.exit()
#        print ('@@@@@',LGRB[idata], Lsolar[idata], Lpart[idata], total)

    sum_data_types = [0] * len(event_classes)						# Calculate Sum (L) and Product (ln L)
    sum_log_data_types = [1] * len(event_classes)
    count_sum = [0] * len(event_classes)
    for idata in range (len(data_types)):
        sum_data_types [1] += Lsolar [idata]
        sum_data_types [2] += Lpart [idata]
        if LGRB [idata] > 0:
            sum_data_types [0] += LGRB [idata]
            sum_log_data_types [0] = sum_log_data_types [0] + math.log (LGRB[idata])
            count_sum [0] += 1
        if Lsolar [idata] > 0:
            sum_data_types [1] += Lsolar [idata]
            sum_log_data_types [1] = sum_log_data_types [1] + math.log (Lsolar[idata])
            count_sum [1] += 1
        if Lpart [idata] > 0:
            sum_data_types [2] += Lpart [idata]
            sum_log_data_types [2] = sum_log_data_types [2] + math.log (Lpart[idata])
            count_sum [2] += 1
     
    for ievent in range (len(event_classes)):						# Normalize by number of non-zero likelihoods
#        sum_data_types [ievent] = sum_data_types [ievent] / len (data_types)              
        sum_data_types [ievent] = sum_data_types [ievent] / count_sum [ievent]
        sum_log_data_types [ievent] = math.exp (sum_log_data_types [ievent]) ** (1 / count_sum [ievent] )

################################
#     Print out results

#    if trig_num > 0:
#        summary_file.write (summary_file_name + '   written by   ' + os.path.basename(sys.argv[0]) + '   ' + str(now) + '\n\n') 
    print('\n %d of %d trigger(s) in %s ' % (trig_num + 1, len(triggers), fileid), end = '   --   ')   # print number of triggers

    print('TriggerID: %d     \nTrigger Time (MDC): %.6f     Trigger Time (MET): %.6f     Trigger Time (UT): %s' 
          % (trigid, trigtime_mdc, trigtime_met, UT))

    print ('\nLatitude = ',lat,'   Longitude = ', lon)
    print ('Time to hv turn on/off = ',transition)

    print('\nHardness ratios --      HXM1: %5.1f     HXM2: %5.1f     SGM: %5.1f \n' % (HXM1_ratio, HXM2_ratio, SGM_ratio))
    print ('Solar RA, dec (radians): ', sun_ra_rad, sun_dec_rad)
    print ('HXM RA, dec (radians)  :   ', hxm_ra_rad, hxm_dec_rad)
    print ('SGM RA, dec (radians)  :', sgm_ra_rad, sgm_dec_rad)
    print('HXM solar angle (theta, degrees) : %.3f' % (solar_hxm[0]))
    print('SGM solar angle (theta, degrees) : %.3f' % (solar_sgm[0]))

    print('\nGOES average activity = ', GOES_avg, '     Variance of GOES activity = ', GOES_var) 
    print('GOES avg +/1 hr = ', GOES_avg_1hr,'      Variance +/- 1 hr = ',GOES_var_1hr)

    print ('\nTime history fit parameters (Gaussian):\nAmplitude = %8d,  Peak pos = %6d,   Bkg = %8d,  S/root(B) = %5.1f,   FWHM = %5.1f,    chi-sq = %5.1f' 
        %(ampl1, peak_pos1, bkg1, s_over_rootB1, FWHM, chisq1) )
    print ('Time history fit parameters (Norris):\nAmplitude = %8d,  Peak pos = %6d,   Bkg = %8d,  S/root(B) = %5.1f,    t90 = %5.1f,    chi-sq = %5.1f' 
           % (ampl2, peak_pos2, bkg2, s_over_rootB2, t90, chisq2) )
    print ('\nLikelihood ratios:')
    space = ' '
    for idata in range (len(data_types)):
        num_spaces = 21 - len (data_types [idata])
        if idata < 5 or idata == 6 or idata == 7 or idata == 9:
            print (num_spaces*space, data_types [idata],":  {:9.2f}".format(likelihood_list[idata]),"     LGRB = {:9.2f}".format(LGRB[idata]),
                "   Lsolar = {:9.2f}".format(Lsolar[idata]), "     Lpart = {:9.2f}".format(Lpart [idata]))
        else:
            print (num_spaces*space, data_types [idata],":  {:.3e}".format(likelihood_list[idata]),"     LGRB = {:9.2f}".format(LGRB[idata]),
                "   Lsolar = {:9.2f}".format(Lsolar[idata]), "     Lpart = {:9.2f}".format(Lpart [idata]))
    print (18*space, 'Avg :             ', end = space)
    for ievent in range (len(event_classes)):
        print ("          %9.2f" % (sum_data_types [ievent]) , end = (ievent + 3) * space)
    print ('\n',14*space, 'Product :             ', end = space)
    for ievent in range (len(event_classes)):
        print ("          %9.2f" % (sum_log_data_types [ievent]) , end = (ievent + 3) * space)


#    print ('\n**************************')
    
 ################################   
# Determine previous classification from event_file (listing of previous events with classifications)
    """
    classification = 'Unknown'
    for i in range (len(event_id_rows)):
#        print(trigid-int(event_id_rows[i]), event_id_rows[i], event_UT_rows[i], event_type_rows[i])   
        if trigid == int(event_id_rows[i]):
#            print (i, event_id_rows[i], event_UT_rows[i], event_type_rows[i])   
            classification = event_type_rows[i]
            break
    """     
    classification = '????'   
    for i in range (len(prior_trig_id)):
#        print (trigid, type (trigid), prior_trig_id[i], type (prior_trig_id[i]))
#        if str(trigid) == prior_trig_id[i]:
        if abs(trigid - int(prior_trig_id[i])) < 2:   # Is trigid equal to or within 1 of prior_trig_id?
            classification = prior_trig_class[i]
            break
    print ('\n   Prior classification:    ',classification,'\n')
    

################################    
# Open and write to output summary file and values file

#    print(' **** 11   Start on summary file, trig_num = ', trig_num, len(triggers))
    if trig_num == 0:
        global summary_file, summary_file_name
        now = datetime.datetime.now()
        d1 = now.strftime("%Y%m%d_%H%M")
        summary_file_name = "Output_summary_file_" + fileid + '_' + d1+'.txt'
        summary_file = open (summary_file_name,"w")
        values_file = open ("values_file.txt","a")
        values_file.write (d1 + '   ' + UT + '   ' + str(trigid) + '\n')
#        values_file = open ("values_file.txt","a")
    
        line = 0
        while line < len(data_types): 
            values_file.write (data_types[line] + ' ')
            line += 1
        values_file.write ('\n')
        summary_file.write (summary_file_name + '   written by   ' + os.path.basename(sys.argv[0]) + '   ' + str(now) + '\n\n') 
    if trig_num > 0:
        summary_file = open (summary_file_name,"a")   
        values_file = open ("values_file.txt","a")

    summary_file.write ('#%d of %d trigger(s) in %s \n' % (trig_num + 1, len(triggers), fileid))   # print number of triggers
    summary_file.write ('   TriggerID: %d \n' % (trigid))
    summary_file.write ('   Trigger Time (MDC): %.6f ' % (trigtime_mdc))
    summary_file.write ('   Trigger Time (MET): %.6f ' % (trigtime_met))
    summary_file.write ('   Trigger Time (UT): %s \n\n' % (UT))

    summary_file.write ('   Latitude = %5.1f   Longitude = %5.1f \n' % (lat, lon))
    summary_file.write ('   Time to hv turn on/off = %d \n\n' % transition)

    summary_file.write ('   Hardness ratios -- HXM1: %5.1f     HXM2: %5.1f     SGM: %5.1f \n\n' % (HXM1_ratio, HXM2_ratio, SGM_ratio))
    summary_file.write ('   Solar RA, dec (radians): %7.3f   %7.3f \n' % (sun_ra_rad, sun_dec_rad))
    summary_file.write ('   HXM RA, dec (radians):   %7.3f   %7.3f   ' % (hxm_ra_rad, hxm_dec_rad))
    summary_file.write ('   SGM RA, dec (radians): %7.3f   %7.3f \n' % (sgm_ra_rad, sgm_dec_rad))
    summary_file.write ('   HXM solar angle (theta, degrees) = %7.3f   ' % (solar_hxm[0]))
    summary_file.write ('   SGM solar angle (theta, degrees) = %7.3f \n\n' % (solar_sgm[0]))

    summary_file.write ('\nGOES average activity %10.3f   Variance of GOES activity = %10.3f' % (GOES_avg, GOES_var) )
    summary_file.write ('\nGOES avg act +/- 1 hr %10.3f   Variance of GOES +/- 1 hr = %10.3f' % (GOES_avg_1hr, GOES_var_1hr) )

    summary_file.write ('   Time history fit parameters (Gaussian):\n      Amplitude = %8d,  Peak pos = %6d,   Bkg = %8d,  S/root(B) = %5.1f,   FWHM = %5.1f,    chi-sq = %5.1f \n' % (ampl1, peak_pos1, bkg1, s_over_rootB1, FWHM, chisq1) )
    summary_file.write ('   Time history fit parameters (Norris):\n      Amplitude = %8d,  Peak pos = %6d,   Bkg = %8d,  S/root(B) = %5.1f,    t90 = %5.1f,    chi-sq = %5.1f \n\n' % (ampl2, peak_pos2, bkg2, s_over_rootB2, t90, chisq2) )

    for idata in range (len(data_types)):
        summary_file.write (data_types [idata])
        num_spaces = 21 - len(data_types[idata])

        for j in range (num_spaces):
            summary_file.write (" ")
        if likelihood_list[idata] >= 100000 or likelihood_list[idata] < .01:
            summary_file.write ("%.3e" % (likelihood_list[idata]))
        else:
            summary_file.write ("%9.2f" % (likelihood_list[idata]))
        summary_file.write ("        LGRB = %9.2f   Lsolar = %9.2f    Lpart = %9.2f \n" % (LGRB[idata], Lsolar [idata], Lpart [idata]) )

        if idata > 9:
            if idata == 10 or idata == 12:
                likelihood_list[idata] = float ( likelihood_list[idata]) * 1.0e8
            if idata == 11 or idata == 13:
                likelihood_list[idata] = float ( likelihood_list[idata]) * 1.0e12
#            values_file.write (str(likelihood_list[idata]))
#            values_file.write ('      ')
#        else:
#            values_file.write (str(int(likelihood_list[idata]*1000)/1000) + '      ')
        values_file.write (str(int(likelihood_list[idata]*1000)/1000) + '      ')

    summary_file.write ('            Avg :                ')
    for ievent in range (len(event_classes)):                
        summary_file.write ("            %9.2f" % (sum_data_types [ievent]))
    summary_file.write ('\nPrior classification:     ' + classification + '\n\n**************************\n\n')
    values_file.write ('Prev. class = ' + classification + '\n')

# Combine individual plots into single file
#    pdfs = ['plot1.pdf', 'plot2.pdf', 'plot31.pdf', 'plot32.pdf', 'plot33.pdf', 'plot34.pdf', 'plot35.pdf', 'plot36.pdf', 'plot4.pdf']
# Why is plot1 redundant?
    if GOES_avg == -1:
        pdfs = ['Plot2.pdf','plot31.pdf', 'plot32.pdf', 'plot33.pdf', 'plot34.pdf', 'plot35.pdf', 'plot36.pdf', 
        'plot4.pdf'] 
    else:
       pdfs = ['Plot2.pdf','plot31.pdf', 'plot32.pdf', 'plot33.pdf', 'plot34.pdf', 'plot35.pdf', 'plot36.pdf', 
        'plot4.pdf', 'plot_solar.pdf', 'GOES_flux.pdf']
#    merger = PdfMerger()
    merger = PdfWriter()
    for pdf in pdfs:
        merger.append (pdf)
    filename = "summary_plot_file_" + fileid + '-'+str(trig_num)+'.pdf'

    merger.write (filename)

    merger.close()
    plt.close()

# Delete plot files after they are merged together.
    os.remove ("Plot2.pdf")
    os.remove ("plot31.pdf")
    os.remove ("plot32.pdf")
    os.remove ("plot33.pdf")
    os.remove ("plot34.pdf")
    os.remove ("plot35.pdf")
    os.remove ("plot36.pdf")
    os.remove ("plot4.pdf")
    if trig_num == len(triggers) - 1 :   	# After final trigger, close plots of geomagnetic activity
        os.remove ("plot_solar.pdf")
        if GOES_avg > -1:
            os.remove ("GOES_flux.pdf")

# After each trigger, close summary file, close values file
#    if trig_num == len(triggers) - 1 :
    summary_file.close()
    values_file.close()
#    print(' **** 12   Return from writing, trig_num = ', trig_num, len(triggers))
    return


        
if __name__ == "__main__":
    main(fileid, chan_type, tmin, tmax)




