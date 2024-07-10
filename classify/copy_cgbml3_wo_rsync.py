#!/usr/bin/env python
# coding: utf-8

"""
Script Name: copy_cgbml3_wo_rsync.py
Description: This script get CGBM data from the data server in Aoyama Gakuin University. Compatible with Windows
Author: Yuta Kawakubo
Date: 2023-12-15
Version: 1.0
"""

import sys
import datetime
# import subprocess
from dotenv import load_dotenv
import os
import paramiko

# Input line
args = sys.argv

if len(args) < 2:
    print("Usage: python ./copy_cgbml3_wo_rsync.py [fileid yyyymmdd]")
    sys.exit(1)

fileid = str(args[1]) # yyyymmdd


def main(fileid):
    """
        This function executes rsync commands to get CGBM data.

        Args:
            fileid (str): date of the data (yyyymmdd)
    """

    # .env is required to set up the environments
    load_dotenv()
    trig_datetime = datetime.datetime.strptime(fileid, '%Y%m%d')

    remote_dir = os.getenv("REMOTE_DIR")
    user = os.getenv("REMOTE_USER")
    host = os.getenv("REMOTE_HOST")
    data_dir = os.getenv("LOCAL_DATA_DIR")
    pbub_key = os.getenv("RSA_KEY_PATH")

    transport = paramiko.Transport((host, 22))
    transport.connect(username=user, pkey=paramiko.RSAKey(filename=pbub_key))

    sftp = paramiko.SFTPClient.from_transport(transport)

    if not os.path.exists(data_dir):
        os.makedirs(data_dir)

    # In the case of we want data of multiple days
    for tdelta in [0]:
        temp_datetime = trig_datetime + datetime.timedelta(days=tdelta)
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
        for tail in data_kinds:
            filename = 'cgbm_%s.%s' % (fileid, tail)
            remote_file_path = '%s/20%02d/20%02d%02d%02d/%s/%s' % (
                remote_dir,
                temp_year,
                temp_year,
                temp_month,
                temp_day,
                data_type,
                filename,
            )

            local_file_path = '%s/%s' % (data_dir, filename)

            print('Copying %s...' % filename)
            try:
                sftp.get(remote_file_path, local_file_path)
                print(f'Successfully copied: {filename}')
            except paramiko.SSHException as e:
                print(f'Error copying {filename}: {str(e)}')

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
            remote_file_path = '%s/20%02d/20%02d%02d%02d/%s/%s' % (
                remote_dir,
                temp_year,
                temp_year,
                temp_month,
                temp_day,
                data_type,
                filename,
            )

            local_file_path = '%s/%s' % (data_dir, filename)
            print('Copying %s...' % filename)

            try:
                sftp.get(remote_file_path, local_file_path)
                print(f'Successfully copied: {filename}')
            except paramiko.SSHException as e:
                print(f'Error copying {filename}: {str(e)}')

    sftp.close()
    transport.close()

if __name__ == "__main__":
    main(fileid)