def define_models():
     """ Used to define the input directories, the list of models, 
         institute and associated colors in all scripts

     Usage: 
         dir_MIPkit,dir_models,mod,inst,mcolor = define_models()
     """

     dir_MIPkit = './MIPkit' # contains MIPkit-A and MIPkit-W
     dir_models = './DATA' # contains all the model netcdf files
         
     mod    = ['NEMO4.0'      , 'MITgcm'      , 'ROMS' , 'FESOM2'      ]
     inst   = ['IGE-CNRS-UGA' , 'UCLA-UMD'    , 'UTAS' , 'AWI'         ]
     mcolor = ['deepskyblue'  , 'darkmagenta' , 'gold' , 'forestgreen' ]  # 'silver', 'darkorange', 'mediumblue', 'gray'

     return dir_MIPkit,dir_models,mod,inst,mcolor
