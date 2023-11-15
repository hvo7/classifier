#!/usr/bin/env python
# coding: utf-8

"""
Script Name: process_trigger.py
Description: This script make some plots and get information to classify triggers.
Author: Yuta Kawakubo
Date: 2023-11-14
Version: 1.0
"""

from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
from astropy.time import Time
from astropy.time import TimeDelta
import astropy.units as u
import cgutil as cg
from scipy.spatial.transform import Rotation as R
from astropy.coordinates import SkyCoord, UnitSphericalRepresentation
from astropy.coordinates import get_sun
import cartopy.crs as ccrs
from dotenv import load_dotenv
import os
import sys

#Input line
args = sys.argv

if len(args) < 4:
    print("Usage: python ./process_trigger.py [fileid yyyymmdd] [trigtime_met] [trigtime_mdc]")
    sys.exit(1)

# Input example
#fileid = '20231004'
#trigtime_met = 749723025.8627577
#trigtime_mdc = 1380443026.7523124

fileid = str(args[1]) #yyyymmdd
trigtime_met = float(args[2])
trigtime_mdc = float(args[3])
trigid = int(trigtime_mdc)

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


def plt_lc_single_det(th_data, trigid, trigtime_met, det):
    """
        This function makes plots of TH light curves and save them as png files.

        Args:
            th_data (table): the table which is an output of get_th_data.
            trigid (int): integer of the trigger time in MDC time
            trigtime_met (float): trigger time in MET
            det (str): detector name (HXM1/HXM2/SGM)

        Returns:
            None
    """

    if det not in ['HXM1', 'HXM2', 'SGM']:
        print('%s is not expected value' % det)
        return None
    fig, cv = plt.subplots(8, 1, sharex=True, figsize=(8, 8))
    cv[0].plot(
        th_data['TIME'] - trigtime_met,
        th_data['%sH_RATE' % (det)].data[:, 0],
        color='k',
        drawstyle='steps-mid',
        label='High ch0'
    )

    cv[1].plot(
        th_data['TIME'] - trigtime_met,
        th_data['%sH_RATE' % (det)].data[:, 1],
        color='k',
        drawstyle='steps-mid',
        label='High ch1'
    )

    cv[2].plot(
        th_data['TIME'] - trigtime_met,
        th_data['%sH_RATE' % (det)].data[:, 2],
        color='k',
        drawstyle='steps-mid',
        label='High ch2'
    )

    cv[3].plot(
        th_data['TIME'] - trigtime_met,
        th_data['%sH_RATE' % (det)].data[:, 3],
        color='k',
        drawstyle='steps-mid',
        label='High ch3'
    )

    cv[4].plot(
        th_data['TIME'] - trigtime_met,
        th_data['%sL_RATE' % (det)].data[:, 0],
        color='k',
        drawstyle='steps-mid',
        label='Low ch0'
    )

    cv[5].plot(
        th_data['TIME'] - trigtime_met,
        th_data['%sL_RATE' % (det)].data[:, 1],
        color='k',
        drawstyle='steps-mid',
        label='Low ch1'
    )

    cv[6].plot(
        th_data['TIME'] - trigtime_met,
        th_data['%sL_RATE' % (det)].data[:, 2],
        color='k',
        drawstyle='steps-mid',
        label='Low ch2'
    )

    cv[7].plot(
        th_data['TIME'] - trigtime_met,
        th_data['%sL_RATE' % (det)].data[:, 3],
        color='k',
        drawstyle='steps-mid',
        label='Low ch3'
    )

    for i in range(8):
        cv[i].tick_params(labelsize=16)
        if i == 8:
            break
        cv[i].set_ylabel('[cts/s]', fontsize=16)
        cv[i].legend(loc='upper right', frameon=False, fontsize=16)
    tmin=-100.
    tmax=300.
    cv[0].set_xlim(tmin, tmax)
    cv[0].set_title('Trigger: %d  %s light curves' % (trigid, det), fontsize=16)
    plt.xlabel('Time since %f [s]' % (trigtime_met), fontsize=16)
    plt.subplots_adjust(
        hspace=0,
        top=0.9,
        bottom=0.1,
        left=0.15,
        right=0.95)
    plt.show()
    return None

def plot_orb(data_dir, fileid, trigtime):
    """
        This function makes an orbit plot and save them as png files.

        Args:
            data_dir (str): path to the data directory
            fileid (str): date of the data (yyyymmdd)
            trigtime (float): trigger time in MET
        Returns:
            None
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
    plt.show()

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
    print('Input MET = %.6f' % (trigtime_met))
    print('Nearest MET = %.6f' % (iat_table['TIME'][nearest_index]))
    print('Time Diff. = %.6f' % (iat_table['TIME'][nearest_index]-trigtime_met))
    print()
    return iss2sky


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
    print('%s zenith (R.A., Dec.) = (%.3f, %.3f)' % (det , sky_pos.ra.deg, sky_pos.dec.deg))
    return sky_pos


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
    """

    trigtime_utc = met2utc(trigtime_met)

    # Get solar position
    sun = get_sun(trigtime_utc)
    solar_angle = get_cgbm_theta_phi(iss2sky, sun.ra.deg, sun.dec.deg, det)
    print('%s solar angle (theta, phi) = (%.3f, %.3f)' % (det , solar_angle[0], solar_angle[1]))
    return solar_angle


def get_angles(data_dir, fileid, trigtime_met):
    """
        This function calculates pointing directions and solar angles and print them.

        Args:
            data_dir (str): path to the data directory
            fileid (str): date of the data (yyyymmdd)
            trigtime_met (float): trigger time in MET
    """

    iat_file = '%s/cgbm_%s.iat' % (data_dir, fileid)
    iat_table = Table.read(iat_file, hdu=1)
    iat_table['TIME'] = iat_table['TIME'] - cg.PERIODIC_HBIN

    iss2sky = get_iss_att(trigtime_met, iat_table)
    zenith_hxm = get_cgbm_pointing(iss2sky, 0., 0., 'HXM')
    zenith_sgm = get_cgbm_pointing(iss2sky, 0., 0., 'SGM')

    solar_hxm = get_solar_incident_angle(iss2sky, trigtime_met, 'HXM')
    solar_sgm = get_solar_incident_angle(iss2sky, trigtime_met, 'SGM')

def process_trigger(data_dir, fileid, trigid, trigtime_met, trigtime_mdc):
    """
        This function make plots and get pointing directions & solar angles.

        Args:
            data_dir (str): path to the data directory
            fileid (str): date of the data (yyyymmdd)
            trigid (int): integer of the trigger time in MDC time
            trigtime_met (float): trigger time in MET
            trigtime_mdc (float): trigger time in MDC time
    """

    th_data = get_th_data(data_dir, fileid)

    # Create CGBM light curves
    plt_lc_single_det(th_data, trigid, trigtime_met, 'HXM1')
    plt_lc_single_det(th_data, trigid, trigtime_met, 'HXM2')
    plt_lc_single_det(th_data, trigid, trigtime_met, 'SGM')

    # Get the ISS orbit
    plot_orb(data_dir, fileid, trigtime_met)
    trigtime_utc = met2utc(trigtime_met)

    # Printing the triggerID and trigger time
    print('TriggerID: %d' % (trigid))
    print('Trigger Time (MDC): %.6f' % (trigtime_mdc))
    print('Trigger Time (MET): %.6f' % (trigtime_met))
    print('Trigger Time (UT): %s' % (trigtime_utc.isot))

    # Printing the pointing directions and solar angles
    get_angles(data_dir, fileid, trigtime_met)


def main(fileid, trigid, trigtime_met, trigtime_mdc):
    load_dotenv()
    data_dir = os.getenv("LOCAL_DATA_DIR")
    if not os.path.exists(data_dir):
        print('No data directory')
        sys.exit()

    process_trigger(data_dir, fileid, trigid, trigtime_met, trigtime_mdc)
    confirm = False
    while not confirm:
        close_flag = input('Do you want to exit? (Yes=1, No=0): ')
        try:
            close_flag = int(close_flag)
        except:
            confirm = False
            close_flag = 0
            continue
        if close_flag == 1:
            confirm = True
        else:
            confirm = False
    sys.exit()

if __name__ == "__main__":
    main(fileid, trigid, trigtime_met, trigtime_mdc)
