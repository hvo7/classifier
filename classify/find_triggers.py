#!/usr/bin/env python
# coding: utf-8

"""
Script Name: find_trigger.py
Description: This script finds onboard triggers using TRIG_STATUS in the hk data.
Author: Yuta Kawakubo
Date: 2023-11-14
Version: 1.0
"""

from astropy.table import Table
import cgutil as cg
from dotenv import load_dotenv
import os
import sys

args = sys.argv

if len(args) < 2:
    print("Usage: python ./find_triggers.py [fileid yyyymmdd]")
    sys.exit(1)

fileid = str(args[1]) #yyyymmdd

def main(fileid):
    """
        This function finds onboard triggers as transitions of TRIG_STATUS.

        Args:
            fileid (str): date of the data (yyyymmdd)
    """

    load_dotenv()
    data_dir = os.getenv("LOCAL_DATA_DIR")

    if not os.path.exists(data_dir):
        print('No data directory')
        sys.exit()

    # Loading hk data
    hk_file = '%s/cgbm_%s.hk' % (data_dir, fileid)
    hk_table = Table.read(hk_file, hdu=1)

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

    # Printing the TriggerID and Trigger time
    for row in triggers:
        print('TriggerID = %d' % (int(row['MDCTIME'])))
        print('Trigger Time (MET) = %.6f' % (float(row['TIME'])))
        print('Trigger Time (MDC) = %.6f' % (float(row['MDCTIME'])))
        print()

if __name__ == "__main__":
    main(fileid)
