#!/usr/bin/env python
# coding: utf-8

"""
Script Name: cgutil.py
Description: This script include some constants fpr CGBM analysis.
Author: Yuta Kawakubo
Date: 2023-11-14
Version: 1.0
"""

from astropy.time import Time
from scipy.spatial.transform import Rotation as R

# Offset of the mission elapsed (MET)
T0_STR =  "2000-01-01T00:00:00"
T0 = Time(T0_STR, format="isot", scale="utc")

# Half of the bin width
PERIODIC_HBIN = 0.5
TH_HBIN = 1. / 16.

# Coordinate conversion quaternion
CGBM2CAL = dict.fromkeys(('HXM1', 'HXM2', 'HXM'),
                         R.from_quat([0.9961947, 0., 0.08715574, 0.]))
CGBM2CAL['SGM'] = R.from_quat([1., 0., 0., 0.])
CAL2CGBM = dict.fromkeys(('HXM1', 'HXM2', 'HXM'),
                         CGBM2CAL['HXM'].inv())
CAL2CGBM['SGM'] = CGBM2CAL['SGM'].inv()
CAL2ISS = R.from_quat([-0.003635, -0.001678, -0.706655, 0.707194])
ISS2CAL = CAL2ISS.inv()