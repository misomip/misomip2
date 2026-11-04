from check_utils import check_dims_vars_attrs

for reg in ['A', 'W']:
  status = check_dims_vars_attrs(institute='IGE-CNRS-UGA',model='NEMO3.6',dir_model='.',region=reg)

