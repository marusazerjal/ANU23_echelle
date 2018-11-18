import os,sys,string,glob
from numpy import *
import matplotlib.pyplot as plt
import pandas
from astropy.time import Time

basedir = "/home/yjzhou/echelle_pipeline/"



def create_database():
    columns = ["objectname"]

    database = pandas.DataFrame(["TEST"],columns=columns)
    database['filepath'] = pandas.Series(zeros(len(database)))
    database['hjd'] = pandas.Series(zeros(len(database)))
    database['ra'] = pandas.Series(zeros(len(database)))
    database['dec'] = pandas.Series(zeros(len(database)))
    database['gmag'] = pandas.Series(zeros(len(database)))
    database['exptime'] = pandas.Series(zeros(len(database)))
    database['teff'] = pandas.Series(zeros(len(database)))
    database['logg'] = pandas.Series(zeros(len(database)))
    database['feh'] = pandas.Series(zeros(len(database)))
    database['snr'] = pandas.Series(zeros(len(database)))
    database['bcorr'] = pandas.Series(zeros(len(database)))
    database['telluricRV'] = pandas.Series(zeros(len(database)))
    database['vsini'] = pandas.Series(zeros(len(database)))
    database['lsdRV'] = pandas.Series(zeros(len(database)))
    database['ccfRV'] = pandas.Series(zeros(len(database)))
    for order in range(25):
        database['order'+str(order)] = pandas.Series(zeros(len(database)))

    database.to_csv(basedir+"database.csv",index=False)

def add_to_database(rvobs):
    rvobs = pandas.read_csv(rvobs,header=None)
    
    header = ["objectname","filepath","hjd","ra","dec","gmag","exptime","teff","logg","feh","snr","bcorr","telluricRV","vsini","lsdRV","ccfRV"]
    headerlen = len(header)
    for i in range(len(header),len(rvobs.iloc[0])):
        header.append("order"+str(i-headerlen))

    rvobs.columns = header
    print rvobs

    database = pandas.read_csv(basedir+"database.csv")

    ### check if observation exists
    if min(abs(database["hjd"]-float(rvobs["hjd"]))) < 1./(60*60*24): ### less than 1s:
        print "Observation exists"
        mask = abs(database["hjd"]-float(rvobs["hjd"])) < 1./(60*60*24)
        indx = arange(len(database))
        database = database.drop(indx[mask])

    print "appending observation"
    database = database.append(rvobs,sort=False,ignore_index=True)
    database.to_csv(basedir+"database.csv",index=False)


def add_folder(folder):
    rvlist = sort(glob.glob(folder+"*.rv"))
    for rvobs in rvlist:
        add_to_database(rvobs)
        
if __name__ == "__main__":


    if sys.argv[1] == "--create_new":
        create_database()
    
    if sys.argv[1] == "--add_rv":
        rvobs = sys.argv[2]
        add_to_database(rvobs)

    if sys.argv[1] == "--add_folder":
        folder = sys.argv[2]
        add_folder(folder)
