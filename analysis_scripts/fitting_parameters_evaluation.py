#!/usr/bin/python
#
# Batch multifit analysis.
# 
# Colenso 05_30_18

import glob
import os, time, datetime

import multiprocessing
import signal
import subprocess
import sys
import threading
from tqdm import tqdm 
import shutil 

import numpy
import math
from subprocess import *
import datareader
import re
import arraytoimage
import i3togrid
from PIL import Image
import storm_analysis.sa_library.sa_h5py as saH5Py

import storm_analysis.sa_library.batch_run as batchRun
import storm_analysis.sa_library.datareader as datareader

def fitting_params(expfolder, max_processes, selected_channels, movie_idx, debug=False):
    
    if debug: print('Entered fitting params') 
    # set path to raw data files
    #expfolder = "Z:\\Chenghang\\chenghaz006_1.31.2019\\"
    expfolder = os.path.normpath(expfolder) 
    XML_folder = expfolder + "\\XMLs\\"

    # set directory for raw STORM .DAX files
    input_directory = expfolder + "\\acquisition\\"
    # make directory for bead registration output

    # create new directory for bead alignment files
    if not os.path.exists(input_directory + "bead_registration\\"):
        os.mkdir(input_directory + "bead_registration\\")
    bead_registration = input_directory + "bead_registration\\"

    # create new directory for bin files
    if not os.path.exists(input_directory + "bins\\"):
        os.mkdir (input_directory + "bins\\")
    bin_folder = input_directory + "bins\\"

    if debug: print('Created directories') 
    
    print(len(glob.glob(bin_folder + "647storm_*_mlist.bin")))

    # If there are no previous STORM molecule lists found, start the STORM analysis
    if len(glob.glob(bin_folder + "647storm_*_mlist.bin")) == 0: 
        # Which color channels are you analyzing?   
        #channels = ["750","647","561","488","IRbead","Visbead"]
        #channels = ["750", "647","561","488"]
        channels = selected_channels 
        cmd_lines = []
        # setup automated STORM analysis
        for channel in (channels):
                    
         #for  i in range (iterations):
            base = str(channel)
                       
            # set output directory for analyzed molecule lists
            output_directory = input_directory + "bins\\"
            # set path to analysis parameters
            
            multi_xml = XML_folder + base + ".xml"
            # set minimum length of .DAX file for analysis
            
            
            minimum_length = 100
            # create list of all .DAX files in STORM directory (all if -1 otherwise select specific frames) 
            if movie_idx >= 0:
                dax_files = glob.glob(input_directory + base + "*{}.dax".format(movie_idx))
            else:
                dax_files = glob.glob(input_directory + base + "*{}.dax")
            # dax_files =glob.glob(input_directory + base + "*{}.dax".format(movie_idx)) if movie_idx >= 0 else glob.glob(input_directory + base + "*{}.dax")
           
            
            # set path to multi-fit analysis code
            mufit_exe = "C:/Program Files/Python36/Lib/site-packages/storm_analysis/sCMOS/scmos_analysis.py"

            # start processes
            # for each STORM movie
            for movie_file in dax_files:
                # read the file information and determine if it is long enough
                movie_obj = datareader.inferReader(movie_file)
                if(movie_obj.filmSize()[2] > minimum_length):
                    # set file name parameters for jobs
                    print("Analyzing:", movie_file,"with xml file",multi_xml)
                    basename = os.path.basename(movie_file)
                    mlistname = output_directory + "/" + basename[:-4] + "_mlist.hdf5"
                    # start fitting
                    cmd_lines.append(['python', mufit_exe,
                                      "--movie", movie_file,
                                      "--bin", mlistname,
                                      "--xml", multi_xml])
            # run the fitting in a batch with the specified number of concurrent processes
        print(cmd_lines)
        if debug:
            import pdb; pdb.set_trace() 
        
        batchRun.batchRun(cmd_lines, max_processes = max_processes)
                    
    if debug: print('Ran Fitting!')
    
    # where are your data?
    storm_folder = expfolder + "\\acquisition\\"
    input_directory = storm_folder + "bins\\"
    output_directory = storm_folder + "real_bins\\"

    # Mkdir for the output folder.
    if not os.path.exists(output_directory):
        os.mkdir(output_directory)

    # which imaging channels are you saving?
    channels = ["750storm", "647storm", "561storm","488storm"]

    cmd_lines = []
    
    if debug: print('Automated Analysis') 
    
    # setup automated analysis
    for channel in tqdm(channels):
     #for  i in range (iterations):
        base = str(channel)
     # create list of all .hdf5 files in BIN directory
        h5_files =glob.glob(input_directory + base + "*.hdf5")
     # set path to multi-fit analysis code
        hdf5_to_bin = "C:/Program Files/Python36/Lib/site-packages/storm_analysis/sa_utilities/hdf5_to_bin.py"
        
        if debug:
            import pdb
            pdb.set_trace()

     # start processes
        for h5_file in tqdm(h5_files):
                
            print ("Found:", h5_file)
            h5_filename = os.path.basename(h5_file)
            print("Converting:",h5_filename)
            bin_out_name = output_directory + h5_filename[:-5] + ".bin"
            cmd_lines.append(['python', hdf5_to_bin,
                                "--hdf5", h5_file,
                                "--bin", bin_out_name
                                ])
    if debug:
        import pdb
        pdb.set_trace()
    
    print(cmd_lines)
    batchRun.batchRun(cmd_lines, max_processes = max_processes)
    
    if debug: print('Done!') 


if __name__ == "__main__":
    exp_folder = "C:/Users/Vatsal/QT_Projects/GUI_XY"
    # exp_folder = "C:/Users/Jerry Yang/Desktop/Test_bead_data" 
    selected_channels = ['750']
    max_processes = 50
    fitting_params(exp_folder, max_processes, selected_channels, debug=True)
