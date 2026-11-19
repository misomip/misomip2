def define_models(region,scenario='hind'):
     """ Used to define the input directories, the list of models, 
         institute and associated colors in all scripts
     
         region [string]: 'A' (Amundsen) or 'W' (Weddell)
         scenario [string]: 'hind' (default), 'warm', ctrl, ... 

         Usage: 
             dir_MIPkit,dir_models,mod,inst,mcolor = define_models()
     """

     dir_MIPkit = './MIPkit' # contains MIPkit-A and MIPkit-W
     dir_models = './DATA' # contains all the model netcdf files
        
     # list all models, institutes, and corresponding mcolor
     mod    = ['NEMO4.0'      , 'MITgcm'      , 'ROMS' , 'FESOM2'      , 'MOM6'   , 'MITgcm66j'        , 'MITgcm-67d' , 'MITgcm-67s' ]
     inst   = ['IGE-CNRS-UGA' , 'UCLA-UMD'    , 'UTAS' , 'AWI'         , 'ANU'    , 'DartmouthCollege' , 'BAS'        , 'BAS'        ]
     mcolor = ['deepskyblue'  , 'darkmagenta' , 'gold' , 'forestgreen' , 'silver' , 'darkorange'       , 'mediumblue' , 'mediumblue' ]  # 'gray' 

     # Removing non relevant models (GO BY DECREASING INDEX !!)
     if   ( ( region == 'A' ) & ( scenario == 'hind' ) ): 
        print('WARNING: Not including ',mod[6], inst[6],' for region ',region,' and scenario ',scenario) # BAS_MITgcm-67j
        del mod[6], inst[6], mcolor[6]
        print('WARNING: Not including ',mod[4], inst[4],' for region ',region,' and scenario ',scenario) # ANU_MOM6
        del mod[4], inst[4], mcolor[4]
     elif ( ( region == 'W' ) & ( scenario == 'hind' ) ):
        print('WARNING: Not including ',mod[7], inst[7],' for region ',region,' and scenario ',scenario) # BAS_MITgcm-67s
        del mod[7], inst[7], mcolor[7]
        print('WARNING: Not including ',mod[5], inst[5],' for region ',region,' and scenario ',scenario) # DartmouthCollege_MITgcm66j
        del mod[5], inst[5], mcolor[5]
        print('WARNING: Not including ',mod[2], inst[2],' for region ',region,' and scenario ',scenario) # UTAS_ROMS
        del mod[2], inst[2], mcolor[2]
     elif ( ( scenario == 'warm' ) | ( scenario == 'ctrl' ) ):  
        del mod[2:], inst[2:], mcolor[2:]

     return dir_MIPkit,dir_models,mod,inst,mcolor
