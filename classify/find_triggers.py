#!/usr/bin/env python
# coding: utf-8

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
    load_dotenv()
    data_dir = os.getenv("LOCAL_DATA_DIR")

    if not os.path.exists(data_dir):
        print('No data directory')
        sys.exit()

    hk_file = '%s/cgbm_%s.hk' % (data_dir, fileid)
    hk_table = Table.read(hk_file, hdu=1)
    hk_table['TIME'] = hk_table['TIME'] - cg.PERIODIC_HBIN
    hk_table['MDCTIME'] = hk_table['MDCTIME'] - cg.PERIODIC_HBIN
    hk_table['TRIG_TRANSITION'] = 0
    hk_table['TRIG_TRANSITION'][0] = 0
    hk_table['TRIG_TRANSITION'][1:] = hk_table['TRIG_STATUS'].data[1:] - hk_table['TRIG_STATUS'].data[:-1]

    # I don't know why this is needed.
    bad_mask = (hk_table['TRIG_TRANSITION'] == 0) | (hk_table['TRIG_TRANSITION'] == 1)
    hk_table['TRIG_TRANSITION'][~bad_mask] = 0

    trig_index = hk_table['TRIG_TRANSITION'] == 1

    triggers = hk_table[trig_index]
    print('%d triggers in %s ' % (len(triggers), fileid))
    print()

    for row in triggers:
        print('TriggerID = %d' % (int(row['MDCTIME'])))
        print('Trigger Time (MET) = %.6f' % (float(row['TIME'])))
        print('Trigger Time (MDC) = %.6f' % (float(row['MDCTIME'])))
        print()

if __name__ == "__main__":
    main(fileid)
