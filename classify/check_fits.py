#!/usr/bin/env python
# coding: utf-8

"""
Script Name: check_fits.py
Description: This script checks the header and table in FITS data.
Since rows and columns of some data are too much to display,
this script cannot get all data in FITS data.

Author: Yuta Kawakubo
Date: 2023-12-14
Version: 1.0
"""

from astropy.table import Table
from dotenv import load_dotenv
import os
import sys
import traceback
from astropy.io import fits

#Input line
args = sys.argv

if len(args) != 1:
    print("Usage: python ./check_fits.py")
    sys.exit(1)

def main():
    load_dotenv()
    data_dir = os.getenv("LOCAL_DATA_DIR")

    auxil = {
        'att':'.att',
        'dt':'.dt',
        'gti':'.gti',
        'hk':'.hk',
        'iat':'.iat',
        'orb':'.orb',
        'tim':'.tim'
        }

    monitor = {
        'hx1_hth':'_hx1_hth.fits',
        'hx1_lth': '_hx1_lth.fits',
        'hx1_hph': '_hx1_hph.fits',
        'hx1_lph': '_hx1_lph.fits',
        'hx2_hth': '_hx2_hth.fits',
        'hx2_lth': '_hx2_lth.fits',
        'hx2_hph': '_hx2_hph.fits',
        'hx2_lph': '_hx2_lph.fits',
        'sgm_hth': '_sgm_hth.fits',
        'sgm_lth': '_sgm_lth.fits',
        'sgm_hph': '_sgm_hph.fits',
        'sgm_lph': '_sgm_lph.fits',
    }

    continue_flag = 1
    while continue_flag == 1:
        fileid = input('Please input date of files (yyyymmdd):')
        print('fileid = %s' % fileid)
        print('auxil data:')
        print('att, dt, gti, hk, iat, orb, tim')
        print('monitor data:')
        print('hx1_hth, hx1_lth, hx1_hph, hx1_lph')
        print('hx2_hth, hx2_lth, hx2_hph, hx2_lph')
        print('sgm_hth, sgm_lth, sgm_hph, sgm_lph')
        print('')
        data_type = input('Which data do you want to check:')
        if data_type in auxil.keys():
            filename = 'cgbm_%s%s' % (fileid, auxil[data_type])

        elif data_type in monitor.keys():
            filename = 'cgbm_%s%s' % (fileid, monitor[data_type])

        else:
            print('data_type is not correct')
            continue

        fullpath = '%s/%s' % (data_dir, filename)

        if not os.path.exists(fullpath):
            print('No data are found')

        try:
            hdul = fits.open(fullpath)

        except Exception as e:
            print(f"An error occurred: {e}")
            traceback.print_exc()
            continue

        print('There are %d extensions in the FITS data.' % (len(hdul)))
        for i, header in enumerate(hdul):
            ext_name = header.header.get("EXTNAME", 'PRIMARY')
            print('%d:  %s' % (i, ext_name))

        print('')
        ext_num = int(input('Which extension do you want to see?:'))

        try:
            header = fits.getheader(fullpath, ext=ext_num)
            table = Table.read(fullpath, hdu=ext_num)

        except Exception as e:
            print(f"An error occurred: {e}")
            continue

        print('')
        print('#Header keywords:')
        for key, value in header.items():
            print(f"{key}: {value}")

        print('')
        print('#Columns:')
        print(table.colnames)

        print('')
        print('#Table (the first five rows)')
        print(table[:5])
        continue_flag = int(input('Do you want to continue? (Yes:1, No:0):'))

if __name__ == "__main__":
    main()