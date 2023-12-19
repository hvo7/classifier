#!/usr/bin/env python
# coding: utf-8

"""
Script Name: make_spectrum.py
Description: This script makes spectra of each detector from PH data and save images of spectra.

Author: Yuta Kawakubo
Date: 2023-12-19
Version: 1.0
"""

from astropy.table import Table, vstack
import numpy as np
import cgutil as cg
import matplotlib.pyplot as plt
from dotenv import load_dotenv
import os
import sys

#Input line
args = sys.argv

if len(args) != 4:
    print("Usage: python ./make_spectrum.py fileid(yyyymmdd) chan_type(PHA/PI) ref_met(offset MET)")
    sys.exit(1)

else:
    fileid = str(args[1])
    chan_type = str(args[2])
    ref_met = float(args[3])


def get_phdata(fileid, data_dir, chan_type):
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


def main(fileid, chan_type, ref_met):
    # Input example
    #fileid = '20231110'
    #chan_type = 'PHA'
    #ref_met = 752966385.210472

    load_dotenv()
    data_dir = os.getenv("LOCAL_DATA_DIR")
    ph_table = get_phdata(fileid, data_dir, chan_type)
    ph_table['ZTIME'] = ph_table['TIME'] - cg.PH_HBIN - ref_met


    fig, cv = plt.subplots(3, 1, sharex=True, figsize=(8, 6))
    cv[0].set_xlim(-1800, 1800)
    cv[0].plot(ph_table['ZTIME'], np.sum(ph_table['HXM1L_COUNTS'], axis=1), drawstyle='steps-mid', label='HXM1')
    cv[1].plot(ph_table['ZTIME'], np.sum(ph_table['HXM2L_COUNTS'], axis=1), drawstyle='steps-mid', label='HXM2')
    cv[2].plot(ph_table['ZTIME'], np.sum(ph_table['SGML_COUNTS'], axis=1), drawstyle='steps-mid', label='SGM')

    cv[2].set_xlabel('Time [s] since %.6f' % (ref_met), fontsize=14)
    plt.subplots_adjust(hspace=0, top=0.95, bottom=0.1, left=0.15, right=0.95)

    for ax in cv:
        ax.axvline(x=0, color='red', linestyle='--')
        ax.tick_params(labelsize=14)
        ax.set_ylabel('Counts/s', fontsize=14)

    plt.legend(loc='best')
    plt.show()

    tmin = float(input("Please input the start time:"))
    tmax = float(input("Please input the end time:"))

    tmask = (tmin <= ph_table['ZTIME']) & (ph_table['ZTIME'] <= tmax)
    roi = ph_table[tmask]
    row_num = float(len(roi))

    # HXM1H
    fig = plt.figure(figsize=(8, 6))
    plt.plot(
        np.arange(cg.PH_BINNUM_HIGH),
        np.sum(roi['HXM1H_RATE'], axis=0) / cg.high_ph_width() / row_num,
        drawstyle='steps-mid')
    plt.yscale('log')
    plt.xlabel('PH channel')
    plt.ylabel('Counts /s /ADC channel')
    plt.title('HXM1H spectrum')
    fig.savefig('%s_hxm1h_spec.png' % (fileid))
    plt.show()

    # HXM2H
    fig = plt.figure(figsize=(8, 6))
    plt.plot(
        np.arange(cg.PH_BINNUM_HIGH),
        np.sum(roi['HXM2H_RATE'], axis=0) / cg.high_ph_width() / row_num,
        drawstyle='steps-mid')
    plt.yscale('log')
    plt.xlabel('PH channel')
    plt.ylabel('Counts /s /ADC channel')
    plt.title('HXM2H spectrum')
    fig.savefig('%s_hxm2h_spec.png' % (fileid))
    plt.show()

    # SGMH
    fig = plt.figure(figsize=(8, 6))
    plt.plot(
        np.arange(cg.PH_BINNUM_HIGH),
        np.sum(roi['SGMH_RATE'], axis=0) / cg.high_ph_width() / row_num,
        drawstyle='steps-mid')
    plt.yscale('log')
    plt.xlabel('PH channel')
    plt.ylabel('Counts /s /ADC channel')
    plt.title('SGMH spectrum')
    fig.savefig('%s_sgmh_spec.png' % (fileid))
    plt.show()

    # HXM1L
    fig = plt.figure(figsize=(8, 6))
    plt.plot(
        np.arange(cg.PH_BINNUM_LOW),
        np.sum(roi['HXM1L_RATE'], axis=0) / cg.low_ph_width() / row_num,
        drawstyle='steps-mid')
    plt.yscale('log')
    plt.xlabel('PH channel')
    plt.ylabel('Counts /s /ADC channel')
    plt.title('HXM1L spectrum')
    fig.savefig('%s_hxm1l_spec.png' % (fileid))
    plt.show()

    # HXM2H
    fig = plt.figure(figsize=(8, 6))
    plt.plot(
        np.arange(cg.PH_BINNUM_LOW),
        np.sum(roi['HXM2L_RATE'], axis=0) / cg.low_ph_width() / row_num,
        drawstyle='steps-mid')
    plt.yscale('log')
    plt.xlabel('PH channel')
    plt.ylabel('Counts / ADC channel')
    plt.title('HXM2L spectrum')
    fig.savefig('%s_hxm2l_spec.png' % (fileid))
    plt.show()

    # SGMH
    fig = plt.figure(figsize=(8, 6))
    plt.plot(
        np.arange(cg.PH_BINNUM_LOW),
        np.sum(roi['SGML_RATE'], axis=0) / cg.low_ph_width() / row_num,
        drawstyle='steps-mid')
    plt.yscale('log')
    plt.xlabel('PH channel')
    plt.ylabel('Counts /s /ADC channel')
    plt.title('SGML spectrum')
    fig.savefig('%s_sgml_spec.png' % (fileid))
    plt.show()

if __name__ == "__main__":
    main(fileid, chan_type, ref_met)
