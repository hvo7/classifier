Description: This plaintext shows the usage of Python scripts.
Author: Yuta Kawakubo
Date: 2023-11-14

# How to get cgbm data
You need to create .env file to define the environments (REMOTE_HOST, REMOTE_USER, REMOTE_DIR, LOCAL_DATA_DIR)
for the data server.

(venv) ykawakubo@procyon:classify % python ./copy_cgbml3.py
Usage: python ./copy_cgbml3.py [fileid yyyymmdd]
(venv) ykawakubo@procyon:classify % python ./copy_cgbml3.py 20231110 #yyyymmdd


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


