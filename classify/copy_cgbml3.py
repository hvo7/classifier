#!/usr/bin/env python
# coding: utf-8

"""
Script Name: copy_cgbml3.py
Description: This script get CGBM data from the data server in Aoyama Gakuin University. Compatible with Linux 
Author: Yuta Kawakubo
Date: 2023-11-14
Version: 1.0
"""

import sys
import datetime
import subprocess
from dotenv import load_dotenv
import os

#Input line
args = sys.argv

if len(args) < 2:
    print("Usage: python ./copy_cgbml3.py [fileid yyyymmdd]")
    sys.exit(1)

fileid = str(args[1]) #yyyymmdd

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

    if not os.path.exists(data_dir):
        os.makedirs(data_dir)

    # In the case of we want data of multiple days
    for tdelta in [0] :
        temp_datetime = trig_datetime + datetime.timedelta(days=tdelta)
        temp_year = int(temp_datetime.strftime('%y'))
        temp_month = int(temp_datetime.strftime('%m'))
        temp_day = int(temp_datetime.strftime('%d'))
        temp_hour = int(temp_datetime.strftime('%H'))
        data_type = 'auxil'
        cmd_args = [
                'rsync',
                '-av',
                '-u',
                 '--progress',
                '--timeout=10',
                '--update',
                '%s@%s:%s/20%02d/20%02d%02d%02d/%s/*' % (
                    user,
                    host,
                    remote_dir,
                    temp_year,
                    temp_year,
                    temp_month,
                    temp_day,
                    data_type),
                data_dir]
        subprocess.check_call(cmd_args)
        result = subprocess.run(cmd_args, text=True, capture_output=True)
        if result.returncode != 0:
            print("rsync for %s failed" % data_type)
            print(f"Error output:\n{result.stderr}")

        data_type = 'monitor'
        cmd_args = [
            'rsync',
            '-av',
            '-u',
            '--progress',
            '--timeout=10',
            '--update',
            '%s@%s:%s/20%02d/20%02d%02d%02d/%s/*' % (
                user,
                host,
                remote_dir,
                temp_year,
                temp_year,
                temp_month,
                temp_day,
                data_type),
            data_dir]
        # print(cmd_args)
        subprocess.check_call(cmd_args)
        result = subprocess.run(cmd_args, text=True, capture_output=True)
        if result.returncode != 0:
            print("rsync for %s failed" % data_type)
            print(f"Error output:\n{result.stderr}")

    print("rsync has finished")

if __name__ == "__main__":
    main(fileid)