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
import storm_analysis.Name_Correction as N_C

from .gen_bead_warp_hcam import * 

def batch_fitting(expfolder, max_processes, selected_channels, debug=False):
    # set path to raw data files
    #expfolder = "Z:\\Chenghang\\chenghaz006_1.31.2019\\"
    expfolder = os.path.normpath(expfolder)
    XML_folder = expfolder + "\\XMLs\\"

    # set directory for raw STORM .DAX files
    input_directory = expfolder + "\\acquisition\\"

    #Correct the file names if there is less than 11 images taken. (might need modification if there are more than 100 images in the future)
    # h5_files =glob.glob(input_directory + "647storm_" + "*.dax")
    # test = os.path.basename(h5_files[0])
    # if len(test) == 14:
        # N_C.CorrectIt(input_directory)


    # make directory for bead registration output
    # create new directory for bead alignment files
    if not os.path.exists(input_directory + "bead_registration\\"):
        os.mkdir(input_directory + "bead_registration\\")
    bead_registration = input_directory + "bead_registration\\"

    # create new directory for bin files
    if not os.path.exists(input_directory + "bins\\"):
        os.mkdir (input_directory + "bins\\")
    bin_folder = input_directory + "bins\\"

    # If there are no previous STORM molecule lists found, start the STORM analysis
    cmd_lines = []
    channels = selected_channels 
    
    for channel in channels:
        
        if len(glob.glob(bin_folder + "{}storm_*_mlist.bin".format(channel))) == 0: 
        
            base = str(channel)
            # set output directory for analyzed molecule lists
            output_directory = input_directory + "bins\\"
            # set path to analysis parameters
            multi_xml = XML_folder + base + ".xml"
            # set minimum length of .DAX file for analysis
            minimum_length = 0
            # create list of all .DAX files in STORM directory
            dax_files =glob.glob(input_directory + base + "*.dax")
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
        batchRun.batchRun(cmd_lines, max_processes = max_processes)

    # If there are no previous bead molecule lists found, start the bead analysis
    if not os.path.isfile(bin_folder + "Regbead_1_00_647mlist.bin"):

        # set bead fitting channels
        channel = "Regbead"

        # setup automated bead analysis
        base = str(channel)
        # set output directory for analyzed molecule lists
        output_directory = input_directory + "bins\\"
        # set path to analysis parameters
        xml = ["Visbead_647.xml","Visbead_561.xml","Visbead_488.xml","IRbead_750.xml","IRbead_647.xml"]
       
        # set minimum length of .DAX file for analysis
        minimum_length = 1
     
        # create list of all .DAX files in bead directory
        dax_files =glob.glob(input_directory + base + "*.dax")
        #print (dax_files)
        
        # set path to multi-fit analysis code
        mufit_exe = "C:/Program Files/Python36/Lib/site-packages/storm_analysis/sCMOS/scmos_analysis.py"

        # start processes
        cmd_lines = []
        for movie_file in dax_files:
            movie_name = os.path.basename(movie_file[:-4])
            # run analysis on Visbeads
            movie_obj = datareader.inferReader(movie_file)
            if(movie_obj.filmSize()[2] > minimum_length):

                basename = os.path.basename(movie_file)
                mlistname = output_directory + "/" + basename[:-4]
                # fit 647 portion of Visbead dax file
                cmd_lines.append(['python', mufit_exe,
                                  "--movie", movie_file,
                                  "--bin", mlistname + "_Vis647mlist.hdf5",
                                  "--xml", XML_folder + xml[0]])
                print("Analyzing:", movie_file,"with xml file",xml[0])
                # fit 561 portion of Visbead dax file
                cmd_lines.append(['python', mufit_exe,
                                  "--movie", movie_file,
                                  "--bin", mlistname + "_561mlist.hdf5",
                                  "--xml", XML_folder + xml[1]])
                print("Analyzing:", movie_file,"with xml file",xml[1])
                # fit 488 portion of Visbead dax file
                cmd_lines.append(['python', mufit_exe,
                                  "--movie", movie_file,
                                  "--bin", mlistname + "_488mlist.hdf5",
                                  "--xml", XML_folder + xml[2]])
                print("Analyzing:", movie_file,"with xml file",xml[2])

                mlistname = output_directory + "/" + basename[:-4]
                # fit 750 portion of IRbead dax file
                cmd_lines.append(['python', mufit_exe,
                                  "--movie", movie_file,
                                  "--bin", mlistname + "_750mlist.hdf5",
                                  "--xml", XML_folder + xml[3]])
                print("Analyzing:", movie_file,"with xml file",xml[3])

                mlistname = output_directory + "/" + basename[:-4]
                # fit 750 portion of IRbead dax file
                cmd_lines.append(['python', mufit_exe,
                                  "--movie", movie_file,
                                  "--bin", mlistname + "_IR647mlist.hdf5",
                                  "--xml", XML_folder + xml[4]])
                print("Analyzing:", movie_file,"with xml file",xml[4])
                
        batchRun.batchRun(cmd_lines, max_processes = max_processes)

    # where are your data?
    storm_folder = expfolder + "\\acquisition\\"
    input_directory = storm_folder + "bins\\"
    output_directory = storm_folder + "real_bins\\"

    # Mkdir for the output folder.
    if not os.path.exists(output_directory):
        os.mkdir(output_directory)

    # which imaging channels are you saving?
    channels = ["750storm", "647storm", "561storm","488storm"]

    # setup automated analysis
    for channel in channels:
     #for  i in range (iterations):
        base = str(channel)
     # create list of all .hdf5 files in BIN directory
        h5_files =glob.glob(input_directory + base + "*.hdf5")
     # set path to multi-fit analysis code
        hdf5_to_bin = "C:/Program Files/Python36/Lib/site-packages/storm_analysis/sa_utilities/hdf5_to_bin.py"

     # start processes
        cmd_lines = []
        for h5_file in h5_files:
                
            print ("Found:", h5_file)
            h5_filename = os.path.basename(h5_file)
            print("Converting:",h5_filename)
            bin_out_name = output_directory + h5_filename[:-5] + ".bin"
            cmd_lines.append(['python', hdf5_to_bin,
                                "--hdf5", h5_file,
                                "--bin", bin_out_name
                                ])
        batchRun.batchRun(cmd_lines, max_processes = max_processes)
            
            
    # Convert the .hdf5 files to .tiff images. 
    storm_image_scale = int(10)
    # where do you want to save the data?
    output_directory = expfolder + "\\stormtiffs\\"
    # Mkdir for the output folder. --- 7.30.2018. 
    if not os.path.exists(output_directory):
        os.mkdir(output_directory)
    # which imaging channels are you saving?
    channels = ["750storm", "647storm", "561storm","488storm"]
    # what sigma to use for rendering images in [channels]?
    sigma_channels = [1,0.5,0.5,1]
    # setup automated analysis
    for channel in channels:
     #for  i in range (iterations):
        base = str(channel)
     # create list of all .hdf5 files in BIN directory
        h5_files =glob.glob(input_directory + base + "*.hdf5")
     # set path to multi-fit analysis code
        hdf5_to_image = "C:/Program Files/Python36/Lib/site-packages/storm_analysis/sa_utilities/hdf5_to_image.py"

     # start processes
        cmd_lines = []
        for h5_file in h5_files:
            if base == "750storm":
                sigma = sigma_channels[0]
            elif base == "647storm":
                sigma = sigma_channels[1]
            elif base == "561storm":
                sigma = sigma_channels[2]
            elif base == "488storm":
                sigma = sigma_channels[3]
                
            print ("Found:", h5_file)
            h5_filename = os.path.basename(h5_file)
            print("Analyzing:",h5_filename,"with sigma ",str(sigma))
            tiff_out_name = output_directory + h5_filename[:-5] + ".tiff"
            cmd_lines.append(['python', hdf5_to_image,
                                "--image", tiff_out_name,
                                "--bin", h5_file,
                                "--scale", str(storm_image_scale),
                                "--sigma", str(sigma)])
        batchRun.batchRun(cmd_lines, max_processes = max_processes)
            

    # perform bead registration
        
    # MATLAB Code 
    # print ("starting matlab gen_bead_warp analysis")
    # if not os.path.isfile(bead_registration + 'Cumulative_distribution_for_registration_self.tif'):
    #    matsub = ("""matlab -nosplash -nodisplay -r "arg1='""" + expfolder + """'; gen_bead_warp_hcam" """)
    #    subprocess.Popen(matsub ,shell=True)
        
    # while not os.path.isfile(bead_registration + 'Cumulative_distribution_for_registration_self.tif'):
    #     time.sleep(10)
    # print ("Bead alignment is finished")  
    # print("The test part is skipped" + "\n" + "Done! ")


    #The following codes are for testing purpose. No need if everything works well in previous codes. 
    #Note that the following code won't end as expected for a customized because the test images are saved under 'Z:\Colenso\05_25_18_sample5\acquisition\bead_registration'

        #test warping transform on beads (only if two bead passes)
    #if not os.path.isfile(bead_registration + 'Cumulative_distribution_for_registration_1_based_on_2.fig'):
    #    print ("starting matlab_test_bead_warp")
    #    matsub = ("""matlab -nosplash -nodisplay -r "arg1='""" + expfolder + """'; test_bead_warp_hcam" """)
    #    subprocess.Popen(matsub ,shell=True)

    #while not os.path.isfile(bead_registration + 'Cumulative_distribution_for_registration_1_based_on_2.fig'):
    #    time.sleep(10)
    #print ("Bead testing is finished")  
if __name__ == "__main__":
    exp_folder = "C:/Users/Vatsal/QT_Projects/New_Pipeline/4color_test"
    # exp_folder = "C:/Users/Jerry Yang/Desktop/Test_bead_data" 
    selected_channels = ['750storm_']
    max_processes = 50
    batch_fitting(exp_folder, max_processes, selected_channels, debug=True)
