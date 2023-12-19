Description: This plaintext shows the usage of Python scripts.
Author: Yuta Kawakubo
Date: 2023-11-14
Last Update: 2023-12-19

# How to get cgbm data
You need to create .env file to define the environments (REMOTE_HOST, REMOTE_USER, REMOTE_DIR,
LOCAL_DATA_DIR, and RSA_KEY_PATH) to access the data server.

(venv) ykawakubo@procyon:classify % python ./copy_cgbml3.py
Usage: python ./copy_cgbml3.py [fileid yyyymmdd]
(venv) ykawakubo@procyon:classify % python ./copy_cgbml3.py 20231110 #yyyymmdd

Since copy_cgbml3.py needs rsync, and it does not work on Windows, I prepared copy_cgbml3_wo_rsync.py
to do the same thing without rsync on Windows.

(venv) ykawakubo@procyon:classify % python ./copy_cgbml3_wo_rsync.py 20231110
Copying cgbm_20231110.att...
Successfully copied: cgbm_20231110.att
Copying cgbm_20231110.dt...
Successfully copied: cgbm_20231110.dt
Copying cgbm_20231110.gti...
Successfully copied: cgbm_20231110.gti
Copying cgbm_20231110.hk...
Successfully copied: cgbm_20231110.hk
Copying cgbm_20231110.iat...
Successfully copied: cgbm_20231110.iat
Copying cgbm_20231110.orb...
Successfully copied: cgbm_20231110.orb
Copying cgbm_20231110.tim...
Successfully copied: cgbm_20231110.tim
Copying cgbm_20231110_hx1_hth.fits...
Successfully copied: cgbm_20231110_hx1_hth.fits
Copying cgbm_20231110_hx1_lth.fits...
Successfully copied: cgbm_20231110_hx1_lth.fits
Copying cgbm_20231110_hx1_hph.fits...
Successfully copied: cgbm_20231110_hx1_hph.fits
Copying cgbm_20231110_hx1_lph.fits...
Successfully copied: cgbm_20231110_hx1_lph.fits
Copying cgbm_20231110_hx2_hth.fits...
Successfully copied: cgbm_20231110_hx2_hth.fits
Copying cgbm_20231110_hx2_lth.fits...
Successfully copied: cgbm_20231110_hx2_lth.fits
Copying cgbm_20231110_hx2_hph.fits...
Successfully copied: cgbm_20231110_hx2_hph.fits
Copying cgbm_20231110_hx2_lph.fits...
Successfully copied: cgbm_20231110_hx2_lph.fits
Copying cgbm_20231110_sgm_hth.fits...
Successfully copied: cgbm_20231110_sgm_hth.fits
Copying cgbm_20231110_sgm_lth.fits...
Successfully copied: cgbm_20231110_sgm_lth.fits
Copying cgbm_20231110_sgm_hph.fits...
Successfully copied: cgbm_20231110_sgm_hph.fits
Copying cgbm_20231110_sgm_lph.fits...
Successfully copied: cgbm_20231110_sgm_lph.fits

# How to check CGBM FITS data
You may want to see the contents of the CGBM FITS data.
You can load FITS data with astropy.Table, and you can access the table in Python easily.
But, it may need some knowledge of Python and astropy. So, I prepared a simple script (check_fits.py)
to check the structures of CGBM FITS data.
This script shows a list of extensions in the FITS data, header keywords, and a list of columns
in the extension. Although this script shows the first five rows of the data in the extension,
I suggest not to use this script to get values in FITS data for analysis.
If you have HEASOFT, you can see the structure of the FITS data with 'fstruct'.
Also, you can check data with 'fv' or 'fdump'.

(venv) ykawakubo@procyon:classify % python ./check_fits.py
Please input date of files (yyyymmdd):20231110
fileid = 20231110
auxil data:
att, dt, gti, hk, iat, orb, tim
monitor data:
hx1_hth, hx1_lth, hx1_hph, hx1_lph
hx2_hth, hx2_lth, hx2_hph, hx2_lph
sgm_hth, sgm_lth, sgm_hph, sgm_lph
Which data do you want to check:att
There are 2 extensions in the FITS data.
0:  PRIMARY
1:  ATTITUDE
Which extension do you want to see?:1
Header keywords:
XTENSION: BINTABLE
BITPIX: 8
NAXIS: 2
NAXIS1: 50
NAXIS2: 114957
PCOUNT: 0
GCOUNT: 1
TFIELDS: 5
EXTNAME: ATTITUDE
TTYPE1: TIME
TFORM1: D
TUNIT1: s
TTYPE2: MDCTIME
TFORM2: D
TUNIT2: s
TTYPE3: QPARAM
TFORM3: 4D
TTYPE4: ATT_FLAG
TFORM4: B
TTYPE5: ATT_RESID
TFORM5: B
TELESCOP: CALET
INSTRUME: CGBM
OBS_ID: 20231110
MJDREFI: 51544
MJDREFF: 0.00074287037037037
TIMEREF: LOCAL
TASSIGN: SATELLITE
TIMESYS: TT
TIMEUNIT: s
CLOCKAPP: True
DATE-OBS: 2023-11-10T00:01:00
DATE-END: 2023-11-11T00:01:01
TSTART: 752889665.075952
TSTOP: 752976066.41187
HDUCLASS: OGIP
HDUCLAS1: TEMPORALDATA
HDUCLAS2: ASPECT
ORIGIN: JAXA, WCOC
CREATOR: cgbmgenatt version 0.2.0 (2018-02-19)
PROCVER: 4.0.3.2
CALDBVER: hxm151005_sgm151005
SEQPNUM: 1
DATE: 2023-11-10T00:00:00.000
L3CORR: True
CHECKSUM: WY8VZX7VWX7VWX7V
DATASUM: 3607934091
HISTORY: START cgbmgenatt at 2023-11-12T09:05:05
HISTORY:
HISTORY: - previous:
HISTORY: /mnt/cgbmdata2/cgbm_heasarc_archive/cgbmarchL3HEA/obs/2023/20231109/auxi
HISTORY: l/cgbm20231109.att
HISTORY: - infile: /tmp/cgbmarch_mkatt.py0vz36Z/asc_att_20231110.fits
HISTORY: - procnum: 13
HISTORY: - timlistfile: 25684048
HISTORY:
HISTORY: END cgbmgenatt
HISTORY: Corrected by correct_l3.py version 3.2.22 (2023-07-24)
HISTORY: at 2023-11-12T09:07:15
Columns:
['TIME', 'MDCTIME', 'QPARAM', 'ATT_FLAG', 'ATT_RESID']
Table (the first five rows)
       TIME            MDCTIME                          QPARAM                   ATT_FLAG ATT_RESID
        s                 s
----------------- ------------------ ------------------------------------------- -------- ---------
752889665.0759516 1383609600.4341407 -0.13415424525737762 .. -0.5696675181388855      233         2
752889665.5879164 1383609600.9460938 -0.13435865938663483 .. -0.5695774555206299      233         5
 752889666.075005  1383609601.433172 -0.13449132442474365 .. -0.5695376396179199      233         4
752889666.5860319 1383609601.9441874 -0.13448673486709595 .. -0.5695946216583252      233         3
752889667.0774179 1383609602.4355626  -0.13461360335350037 .. -0.569556474685669      233         4
Do you want to continue? (Yes:1, No:0):0

# How to find triggers
(venv) ykawakubo@procyon:classify % python ./find_triggers.py
Usage: python ./find_triggers.py [fileid yyyymmdd]
(venv) ykawakubo@procyon:classify % python ./find_triggers.py 20231110

1 triggers in 20231110

TriggerID = 1383686318
Trigger Time (MET) = 752966385.210472
Trigger Time (MDC) = 1383686318.874000

# How to check the trigger
(venv) ykawakubo@procyon:classify % python ./process_trigger.py
Usage: python ./process_trigger.py [fileid yyyymmdd] [trigtime_met] [trigtime_mdc]

(venv) ykawakubo@procyon:classify % python ./process_trigger.py 20231110 752966385.210472 1383686318.874000

TriggerID: 1383686318
Trigger Time (MDC): 1383686318.874000
Trigger Time (MET): 752966385.210472
Trigger Time (UT): 2023-11-10T21:19:40.210
Input MET = 752966385.210472
Nearest MET = 752966385.210535
Time Diff. = 0.000063

HXM zenith (R.A., Dec.) = (68.508, 60.695)
SGM zenith (R.A., Dec.) = (65.945, 50.797)
HXM solar angle (theta, phi) = (133.430, 338.304)
SGM solar angle (theta, phi) = (142.579, 333.780)
Do you want to exit? (Yes=1, No=0): 1

# How to make spectra
(venv) ykawakubo@procyon:classify % python ./make_spectrum.py
Usage: python ./make_spectrum.py fileid(yyyymmdd) chan_type(PHA/PI) ref_met(offset MET)

(venv) ykawakubo@procyon:classify % python ./make_spectrum.py 20231110 PHA 752966385.210472
Warning: chan_type is not PI.
Warning: Energy calibration was not applied.
Please input the start time:-1200
Please input the end time:1200