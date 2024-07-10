61113#!/usr/bin/env python
# coding: utf-8

"""
Script Name: copy_multiple.py

This script get multiple days of CGBM data from the data server in  
    Aoyama Gakuin University. Modified version of Yuta Kawakubo's 
    copy_cgbml3_wo_rsync.py to copy multiple files. Version calls copy_Mike9 repeatedly
This version reads successive days of data, writes to local file, makes
    summary list of files successfully copied in summary_file_yyyymmdd_hhmm, where time is time at which code is run.
 
Example of usage: python ./copy_multiple1.py 

Author: Mike Cherry
Date: 2023-4-21
Version: 1.0

"""
import sys
import datetime
from dotenv import load_dotenv
import os
import paramiko
from astropy.table import Table
import cgutil as cg
from math import trunc

def main():

    load_dotenv()     # .env is required to set up the environments

    # Open summary file
    global summary_file_name, summary_file, psummary

    now = datetime.datetime.now()
    d1 = now.strftime("%Y%m%d_%H%M")
    summary_file_name = "summary_file_" + d1
    summary_file = open(summary_file_name,"w")
    print ('\nJob started at ',d1)
    d2 = now.strftime("%Y%m%d_%H%M%S")
    psummary = 0


    copy(20220103,1)	#	PART
    copy(20220106,1)	#	SGR
    copy(20220111,1)	#	SGR
    copy(20220112,1)	#	SGR
    copy(20220114,1)	#	SGR
    copy(20220115,1)	#	PART
    copy(20220119,1)	#	PART
    copy(20220120,1)	#	SFLR
    copy(20220124,1)	#	LGRB
    copy(20220126,1)	#	SFLR
    copy(20220130,1)	#	SFLR
    copy(20220205,1)	#	LGRB
    copy(20220207,1)	#	PART
    copy(20220209,1)	#	LGRB
    copy(20220210,1)	#	SGRB
    copy(20220219,1)	#	LGRB
    copy(20220222,1)	#	PART
    copy(20220224,1)	#	PART
    copy(20220306,1)	#	PART
    copy(20220308,1)	#	LGRB
    copy(20220310,1)	#	LGRB
    copy(20220312,1)	#	PART
    copy(20220314,1)	#	PART
    copy(20220320,1)	#	LGRB
    copy(20220322,1)	#	SFLR
    copy(20220328,1)	#	UNKN
    copy(20220329,1)	#	SFLR
    copy(20220402,1)	#	SFLR
    copy(20220403,1)	#	LGRB
    copy(20220408,1)	#	PART
    copy(20220410,1)	#	PART
    copy(20220412,1)	#	PART
    copy(20220414,1)	#	PART
    copy(20220415,1)	#	PART
    copy(20220418,1)	#	PART
    copy(20220419,1)	#	SFLR
    copy(20220420,1)	#	SFLR
    copy(20220421,1)	#	LGRB
    copy(20220422,1)	#	SFLR
    copy(20220423,1)	#	LGRB
    copy(20220425,1)	#	FLCT
    copy(20220427,1)	#	PART
    copy(20220428,1)	#	PART
    copy(20220430,1)	#	SFLR
    copy(20220504,1)	#	SFLR
    copy(20220504,1)	#	SFLR
    copy(20220505,1)	#	SFLR
    copy(20220510,1)	#	PART
    copy(20220511,1)	#	SFLR
    copy(20220514,1)	#	LGRB
    copy(20220519,1)	#	SFLR
    copy(20220520,1)	#	SGR
    copy(20220521,1)	#	LGRB
    copy(20220523,1)	#	SGR
    copy(20220525,1)	#	LGRB
    copy(20220526,1)	#	SFLR
    copy(20220527,1)	#	LGRB
    copy(20220603,1)	#	SGR
    copy(20220606,1)	#	LGRB
    copy(20220609,1)	#	LGRB
    copy(20220612,1)	#	LGRB
    copy(20220614,1)	#	PART
    copy(20220616,1)	#	SGR
    copy(20220618,1)	#	PART
    copy(20220620,1)	#	PART
    copy(20220624,1)	#	LGRB
    copy(20220627,1)	#	PART
    copy(20220628,1)	#	PART
    copy(20220702,1)	#	PART
    copy(20220703,1)	#	PART
    copy(20220704,1)	#	PART
    copy(20220705,1)	#	PART
    copy(20220708,1)	#	SFLR
    copy(20220711,1)	#	LGRB
    copy(20220712,1)	#	PART
    copy(20220715,1)	#	PART
    copy(20220719,1)	#	PART
    copy(20220720,1)	#	PART
    copy(20220722,1)	#	PART
    copy(20220723,1)	#	PART
    copy(20220728,1)	#	PART
    copy(20220731,1)	#	PART
    copy(20220810,1)	#	PART
    copy(20220812,1)	#	PART
    copy(20220814,1)	#	PART
    copy(20220817,1)	#	PART
    copy(20220818,1)	#	UNKN
    copy(20220819,1)	#	PART
    copy(20220821,1)	#	PART
    copy(20220822,1)	#	PART
    copy(20220824,1)	#	PART
    copy(20220825,1)	#	PART
    copy(20220827,1)	#	SFLR
    copy(20220829,1)	#	PGRB
    copy(20221029,1)	#	PART
    copy(20221109,1)	#	SGR
    copy(20221110,1)	#	PART
    copy(20221111,1)	#	UNKN
    copy(20221112,1)	#	UNKN
    copy(20221113,1)	#	SGR
    copy(20221114,1)	#	LGRB
    copy(20221115,1)	#	LGRB
    copy(20221120,1)	#	SGR
    copy(20221121,1)	#	LGRB
    copy(20221124,1)	#	LGRB
    copy(20221130,1)	#	PART
    copy(20221202,1)	#	PART
    copy(20221204,1)	#	PART
    copy(20221209,1)	#	LGRB
    copy(20221214,1)	#	PART
    copy(20221217,1)	#	PART
    copy(20221218,1)	#	PART
    copy(20221219,1)	#	PART
    copy(20221220,1)	#	PART
    copy(20221221,1)	#	PART
    copy(20221222,1)	#	PART
    copy(20221223,1)	#	SFLR
    copy(20221224,1)	#	PART
    copy(20221229,1)	#	SFLR
    copy(20221230,1)	#	SFLR
    copy(20221231,1)	#	PART

# Job duration
    duration = d2.split("_")                    # Job start time
    duration1 = int(duration[1])
    hours = trunc(duration1/10000) # 12
    minutes = trunc((duration1 - hours*10000) /100) #10
    seconds = trunc((duration1 - hours*10000 - minutes*100) ) # 05
    job_start = hours * 3600 + minutes * 60 + seconds
    now = datetime.datetime.now()               # Job ending time
    d3 = now.strftime("%Y%m%d_%H%M%S")
    duration = d3.split("_")          
    duration1 = int(duration[1])
    hours = trunc(duration1/10000) # 12
    minutes = trunc((duration1 - hours*10000) /100) #10
    seconds = trunc((duration1 - hours*10000 - minutes*100) ) # 05
    job_duration = hours * 3600 + minutes * 60 + seconds -job_start
    job_duration_mins = trunc((job_duration/60) *100) / 100
    print ('\n',job_duration, ' secs = ', job_duration_mins, ' mins.')
    summary_file.write ('\nJob duration = '+str(job_duration_mins) + ' mins = ' + str(job_duration) + ' secs.')
    summary_file.close()
    summary_file = open(summary_file_name,"r")
    content = summary_file.read()
    print ("\n Files copied -- ", summary_file_name, ":\n",content)

    exit()

#********************************************************************************************************
def copy(file, psummary):

    fileid = str(file)
    trig_datetime = datetime.datetime.strptime(fileid, '%Y%m%d')
    load_dotenv()

    remote_dir = os.getenv("REMOTE_DIR")
    user = os.getenv("REMOTE_USER")
    host = os.getenv("REMOTE_HOST")
    data_dir = os.getenv("LOCAL_DATA_DIR")
    pbub_key = os.getenv("RSA_KEY_PATH")
    transport = paramiko.Transport((host, 22))
    transport.connect(username=user, pkey=paramiko.RSAKey(filename=pbub_key))

    print ('\ntrig_datetime=', type(trig_datetime), trig_datetime, '\n',
'  remote_dir =', type(remote_dir ),remote_dir,'\n',
'  user=', type(user),user,'\n',
'  host=', type(host),host,'\n',
'  data_dir =', type(data_dir ),data_dir ,'\n',
'  pbub_key =', type(pbub_key ),pbub_key )
   
    sftp = paramiko.SFTPClient.from_transport(transport)

    if not os.path.exists(data_dir):
        os.makedirs(data_dir)
    
    trig_index = []
    trig_time = []

    treps = 1    # Number of days to process (always set to 1 for this version)

    for tdelta in range (0, treps):
        readerror = 0
        temp_datetime = trig_datetime + datetime.timedelta(days = tdelta)

        temp_year = int(temp_datetime.strftime('%y'))
        temp_month = int(temp_datetime.strftime('%m'))
        temp_day = int(temp_datetime.strftime('%d'))
        temp_hour = int(temp_datetime.strftime('%H'))
        data_type = 'auxil'
        data_kinds =[
            'att',
            'dt',
            'gti',
            'hk',
            'iat',
            'orb', 
            'tim',
        ]
        fileid = '20'+str(temp_year * 10000 + temp_month * 100 + temp_day)

        for tail in data_kinds:
            filename = 'cgbm_%s.%s' % (fileid, tail)
            remote_file_path ='%s/20%02d/20%02d%02d%02d/%s/%s' % (        
                remote_dir,
                temp_year,
                temp_year,
                temp_month,
                temp_day,
                data_type,
                filename,
            )

            local_file_path = '%s/%s' % (data_dir, filename)

            if tail == data_kinds[0] and psummary == 0:

                print(' remote_file_path = ', remote_file_path,'\n',
'local_file_path = ', local_file_path, filename)
                summary_file.write(summary_file_name + '\n')
                summary_file.write('remote_file_path = '+ remote_file_path+'\n'+
'local_file_path = '+ local_file_path+ filename+'  '+ '\n')
                psummary = 1

            try:
                sftp.get(remote_file_path, local_file_path)            
                print(f'Successfully copied: {filename}')
#            except paramiko.SSHException as e:
#                print(f'Error copying {filename}: {str(e)}')
            except:
                print(f'Error copying {filename}')
                readerror = 1
                break
        
        if readerror == 1:
            readerror = 0
            continue

        data_type = 'monitor'
        data_kinds = [
            'hx1_hth.fits',
            'hx1_lth.fits',
            'hx1_hph.fits',
            'hx1_lph.fits',
            'hx2_hth.fits',
            'hx2_lth.fits',
            'hx2_hph.fits',
            'hx2_lph.fits',
            'sgm_hth.fits',
            'sgm_lth.fits',
            'sgm_hph.fits',
            'sgm_lph.fits',
        ]
        for tail in data_kinds:
            filename = 'cgbm_%s_%s' % (fileid, tail)
            remote_file_path ='%s/20%02d/20%02d%02d%02d/%s/%s' % (
                remote_dir,
                temp_year,
                temp_year,
                temp_month,
                temp_day,
                data_type,
                filename,
            )

            local_file_path = '%s/%s' % (data_dir, filename)
   
            try:
                sftp.get(remote_file_path, local_file_path)
                print(f'Successfully copied: {filename}')
            except paramiko.SSHException as e:
                print(f'Error copying {filename}: {str(e)}') 

        trig_index,trig_time, readerror = find_triggers_Mike6 (fileid)
#        if (readerror == 1) or (len(trig_index) == 0) :
        if len(trig_index) == 0 :
            continue
 
        summary_file.write(str(tdelta+1)+"  "+fileid+"  "+str(trig_index)+str(trig_time)+"\n")
 
    sftp.close()
    transport.close()
    return

#***************************************************************
"""
Script Name: find_triggers_Mike6.py
Description: This script finds onboard triggers using TRIG_STATUS in the hk data.
Author: Yuta Kawakubo
Date: 2023-11-14
Version: 1.0
Modified by Mike Cherry to allow for multiple files
3/4/2024
"""

def find_triggers_Mike6 (fileid):
    """
        This function finds onboard triggers as transitions of TRIG_STATUS.

        Args:
            fileid (str): date of the data (yyyymmdd)
    """

    load_dotenv()
    data_dir = os.getenv("LOCAL_DATA_DIR")
    print (' ')
    if not os.path.exists(data_dir):
        print('No data directory')
        sys.exit()

    trig_index = []
    trig_time = []

    # Loading hk data
    readerror = 0
    hk_file = '%s/cgbm_%s.hk' % (data_dir, fileid)
#    try:
    hk_table = Table.read(hk_file, hdu=1)
#    except:
#        readerror == 1
#        print('readerror == 1')
#        return trig_index, trig_time, readerror

    # Timestamps of CGBM data are end of the bins.
    # Timestamps were changed to center of the bins
    hk_table['TIME'] = hk_table['TIME'] - cg.PERIODIC_HBIN
    hk_table['MDCTIME'] = hk_table['MDCTIME'] - cg.PERIODIC_HBIN

    # The transition of TRIG_STATUS (from 0 to 1) means onboard trigger happened.
    # Calculating
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
    print('%d triggers in %s ' % (len(triggers), fileid))
    print()

    trig_index = []
    trig_time = []

    # Printing the TriggerID and Trigger time
    indexp = 0
    for row in triggers:
        indexp +=1
        print('TriggerID = %d' % (int(row['MDCTIME'])))
        print('Trigger Time (MET) = %.6f' % (float(row['TIME'])))
        print('Trigger Time (MDC) = %.6f' % (float(row['MDCTIME'])))
        trig_index.append (int(row['MDCTIME']))
        trig_time.append (float(row['TIME']))

    print (trig_index)
    print (trig_time,'\n')

    return trig_index,trig_time, readerror
   
if __name__ == "__main__":
    main()