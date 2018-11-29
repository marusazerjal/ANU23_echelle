

def set_config():
    config = {
        "folder" : "/media/Onion/Data/ANU23echelle/bin2test/",
        "delete_prev" : True,        ### delete all previous reductions
        "binning" : 2,                 ### Y binning 
        "bias" : "bias",              ### bias header
        "flat" : "Flat",              ### flat header
        "arc" : "ThAr",               ### ThAr header
        "norder" : 24,                ### Number of orders to extract ### Default is 24
        "jdheader" : "MJD-OBS",       ### original JD header that's written directly into each obs without modification -- not the BJD!!!
        "init_wave_obj" : "HD96700",  ### Use this object's spectrum for initial wavelength solution
        "run_analysis" : False,         ### Run the RV and spectral typing analysis
        "spectral_library" : "/media/Onion/Data/spectral_library/"   ### path for spectral library
    }


    return config

