import os


full_file_path = os.path.realpath(__file__)
drive, filename = os.path.split(full_file_path)
file_prefix, extention = filename.rsplit('.', 1)
