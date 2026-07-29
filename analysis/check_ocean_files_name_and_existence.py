from check_utils import check_files

for reg in ['A', 'W']:

   status = check_files(institute='IGE-CNRS-UGA',model='NEMO3.6',dir_model='./',region=reg)
