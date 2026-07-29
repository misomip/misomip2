from download_utils import download_MIPkit_from_zenodo

for reg in ['A', 'W']:

  status = download_MIPkit_from_zenodo(region=reg,output_dir='MIPkit/MIPkit-'+reg)
