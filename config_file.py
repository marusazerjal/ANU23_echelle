

def set_config():
    config = {
        "folder" : "/media/Onion/Data/ANU23echelle/20181118/",
        "delete_prev" : True,        ### delete all previous reductions
        "binning" : 1,                 ### Y binning 
        "bias" : "Bias",              ### bias header
        "flat" : "Flat",              ### flat header
        "arc" : "ThAr",               ### ThAr header
        "norder" : 24,                ### Number of orders to extract ### was 17
        "jdheader" : "MJD-OBS",       ### original JD header that's written directly into each obs without modification -- not the BJD!!!
        "init_wave_obj" : "HD1388",  ### Use this object's spectrum for initial wavelength solution
        "run_analysis" : True,         ### Run the RV and spectral typing analysis
        "spectral_library" : "/media/Onion/Data/spectral_library/"   ### path for spectral library
    }


    return config

