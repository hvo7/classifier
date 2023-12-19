#!/usr/bin/env python
# coding: utf-8

"""
Script Name: cgutil.py
Description: This script include some constants fpr CGBM analysis.
Author: Yuta Kawakubo
Date: 2023-11-14
Version: 1.0
"""


import numpy as np
from astropy.time import Time
from scipy.spatial.transform import Rotation as R

# Offset of the mission elapsed (MET)
T0_STR =  "2000-01-01T00:00:00"
T0 = Time(T0_STR, format="isot", scale="utc")

# Half of the bin width
PERIODIC_HBIN = 0.5
TH_HBIN = 1. / 16.
PH_HBIN = 2.0

PH_BINNUM_HIGH = 102
PH_BINNUM_LOW = 410

# Coordinate conversion quaternion
CGBM2CAL = dict.fromkeys(('HXM1', 'HXM2', 'HXM'),
                         R.from_quat([0.9961947, 0., 0.08715574, 0.]))
CGBM2CAL['SGM'] = R.from_quat([1., 0., 0., 0.])
CAL2CGBM = dict.fromkeys(('HXM1', 'HXM2', 'HXM'),
                         CGBM2CAL['HXM'].inv())
CAL2CGBM['SGM'] = CGBM2CAL['SGM'].inv()
CAL2ISS = R.from_quat([-0.003635, -0.001678, -0.706655, 0.707194])
ISS2CAL = CAL2ISS.inv()

# Parameters

# Pedestal
HXM1H_PED = -2.2
HXM1H_PED_ERR = 0.1
HXM2H_PED = -5.25
HXM2H_PED_ERR = 0.04
SGMH_PED = 10.41
SGMH_PED_ERR = 0.03
HXM1L_PED = -1.24
HXM1L_PED_ERR = 0.04
HXM2L_PED = -0.475
HXM2L_PED_ERR = 0.002
SGML_PED = -0.171
SGML_PED_ERR = 0.003

HXM1_AMP_GAIN=30.458575
HXM2_AMP_GAIN=29.935845
SGM_AMP_GAIN=30.381555


LINR_HXML = 0.75
LINR_SGML = 7.0
LINR_HXMH = 0.024
LINR_SGMH = 0.23

GC_LINE = 511.0
GC_LINE_HXM = 1470.0
GC_LINE_SGM = 2200.0


PH_LLIM_HXMH = 57 #>30.2keV
PH_ULIM_HXMH = 89 #<80.9keV
PH_LLIM_HXML = 51 #>81.0keV
PH_ULIM_HXML = 320 #<2004.0keV

PH_LLIM_SGMH = 33 #>101.9keV
PH_ULIM_SGMH = 81 #<657.6keV
PH_LLIM_SGML = 47 #>658.0keV
PH_ULIM_SGML = 332 #<20048.0keV

def high_th_def():
    """
    This function makes ADC edge of TH High gain.
    """
    th_def = np.array(
        [[155, 410],
         [411, 1002],
         [1003, 2026],
         [2027, 4031]])
    return th_def

def low_th_def():
    """
    This function makes ADC edge of TH Low gain.
    """
    th_def = np.array(
        [[80, 119],
          [120, 219],
          [220, 383],
         [384, 4079]])
    return th_def

def high_th_edge():
    """
    This function makes ADC edge of TH High gain.
    """
    th_def = np.array(
        [[155, 411],
         [411, 1003],
         [1003, 2027],
         [2027, 4032]])
    return th_def

def low_th_edge():
    """
    This function makes ADC edge of TH Low gain.
    """
    th_def = np.array(
        [[80, 120],
          [120, 220],
          [220, 384],
         [384, 4080]])
    return th_def



def high_ph():
    """
    This function makes ADC edge of PH High gain.
    """

    bin_1 = np.array([0])
    bin_2 = np.arange(1, 2, 2)
    bin_4 = np.arange(3, 10, 4)
    bin_8 = np.arange(11, 42, 8)
    bin_16 = np.arange(43, 682, 16)
    bin_64 = np.arange(683, 4075, 64)
    bin_21 = np.array([4075, 4096])

    bin_all = np.concatenate(
            [bin_1, bin_2, bin_4, bin_8, bin_16, bin_64, bin_21],
            axis=0)
    return bin_all


def low_ph():
    """
    This function makes ADC edge of PH High gain.
    """

    bin_2 = np.arange(0, 95, 2)
    bin_4 = np.arange(96, 351, 4)
    bin_8 = np.arange(352, 1375, 8)
    bin_16 = np.arange(1376, 4095, 16)
    last = np.array([4096])

    bin_all = np.concatenate([bin_2, bin_4, bin_8, bin_16, last], axis=0)
    return bin_all


def high_ph_width():
    ph_edge = high_ph()
    width = ph_edge[1:] - ph_edge[:-1]
    return width


def low_ph_width():
    ph_edge = low_ph()
    width = ph_edge[1:] - ph_edge[:-1]
    return width


def high_ph_center():
    ph_edge = high_ph()
    width = high_ph_width()
    center = ph_edge[:-1] + width/2.
    return center


def low_ph_center():
    ph_edge = low_ph()
    width = low_ph_width()
    center = ph_edge[:-1] + width/2.
    return center


def func0(x, a):
    """
    func0 is definition 0-order function of x.
    """
    return a + x*0.


def func1(x, a, b):
    """
    func1 is definition 1-order function of x.
    """
    return b*x + a


def func2(x, a, b, c):
    """
    func2 is definition 2-order function of x.
    """
    return c*(x**2.0) + b*x + a


def func3(x, a, b, c, d):
    """
    func3 is definition 3-order function of x.
    """
    return d*(x**3.0) + c*(x**2.0) + b*x + a

def gauss(x,p0, p1, p2, p3, p4):
    y = p0*(np.exp(-((x-p1)**2)/(2*p2))) + p3*x + p4
    return y


def gauss_alt(x, p0, p1, p2, p3, p4):
    y = (
        (2.355*p0)
        / (p1*p2 * np.sqrt(2.0 * np.pi))
        * np.exp((-1.0/2.0) * ((2.355*(x-p1))/(p1 * p2)) ** 2.0)
        + p3*x
        + p4
        )
    return y


def adc2ene_ngc(det, adc_chan):
    if det == 'HXM1H':
        ene = (adc_chan-HXM1H_PED) * LINR_HXMH

    if det == 'HXM2H':
        ene = (adc_chan-HXM2H_PED) * LINR_HXMH

    elif det == 'SGMH':
        ene = (adc_chan-SGMH_PED) * LINR_SGMH

    elif det == 'HXM1L':
        ene = (adc_chan-HXM1L_PED) * LINR_HXML

    elif det == 'HXM2L':
        ene = (adc_chan-HXM2L_PED) * LINR_HXML

    elif det == 'SGML':
        ene = (adc_chan-SGML_PED) * LINR_SGML

    else:
        raise ValueError
        sys.exit()
    return ene

def adc2ene_gc(det, adc_chan):
    if det == 'HXM1H':
        ene = adc_chan * LINR_HXMH

    elif det == 'HXM2H':
        ene = adc_chan * LINR_HXMH

    elif det == 'SGMH':
        ene = adc_chan * LINR_SGMH

    elif det == 'HXM1L':
        ene = adc_chan * LINR_HXML

    elif det == 'HXM2L':
        ene = adc_chan * LINR_HXML

    elif det == 'SGML':
        ene = adc_chan * LINR_SGML

    else:
        raise ValueError
        sys.exit()
    return ene
