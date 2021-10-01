import os
import time
import glob
import sys
import shutil
import datetime
import subprocess
from cproc import *
from reorder import *

def z_align(expfolder):
    exp_folder = os.path.normpath(expfolder) + "\\analysis\\"
    print(exp_folder)
    exp_folder_fixed = exp_folder.replace(os.sep, "/")
    print(exp_folder_fixed)
    ISanalysisfolder = exp_folder + "individual_sections\\"  
    slicenum = len(os.listdir(ISanalysisfolder))

    if not os.path.exists(exp_folder + "unaligned_original\\"):
        print ('copying original data')
        shutil.copytree(exp_folder + "unaligned\\",exp_folder + "unaligned_original\\")

    #for both conv and storm images, normalize intensity across sections
    #and threshold images
    if not os.path.exists(exp_folder + "unaligned\\for_align\\"):
        print ('normalizing alignment channel images')
        matsub = ("""matlab -nosplash -nodisplay -r "arg1='""" +
              exp_folder + """'; wga_norm_and_thresh_488wga" """)
        subprocess.check_call(matsub ,shell=True)
        while not os.path.exists(exp_folder + "unaligned\\for_align\\" +
                                 "%03d" % (slicenum-1) + ".png"):
            print ("...WGA normalization processing...")
            time.sleep(60) 

    #perform rigid align
    if not os.path.exists(exp_folder + "rigid_align"):
        print ('running rigid alignment')
        fijisub = ('C:/Users/Vatsal/Fiji.app/ImageJ-win64 ' +
                   '-Xms50g -Xmx50g -Xincgc -XX:MaxPermSize=256m ' + 
                   '-XX:PermSize=256m -XX:NewRatio=5 -XX:CMSTriggerRatio=50 ' +
                   ' -- --no-splash -macro fiji_z_align1_488wga.py "' +
                   exp_folder_fixed + '"')
        subprocess.check_call(fijisub ,shell=True)
        while not os.path.exists(exp_folder + "rigid_align\\for_align\\" +
                                 "%03d" % (slicenum-1) + ".png"):
            print ("...rigid alignment processing...")
            time.sleep(60) 

    #manually determine angle of rotation and cropping
    if not os.path.exists(exp_folder + "cropped"):
        print ('initializing image crop')
        matsub2 = ("""matlab -nosplash -r "arg1='""" + exp_folder
                   + """'; crop_datastack_488wga" """)
        subprocess.check_call(matsub2 ,shell=True)
        while not os.path.exists(exp_folder + "cropped\\for_align\\" +
                                 "%03d" % (slicenum-1) + ".tif"):
            print ("...cropping images now...")
            time.sleep(60) 

        
    #perform elastic align
    if not os.path.exists(exp_folder + "elastic_align"):
        print ('running elastic alignment')
        fijisub = ('C:/Users/Chenghang/Fiji.app/ImageJ-win64 ' +
                   '-Xms50g -Xmx50g -Xincgc -XX:MaxPermSize=256m ' + 
                   '-XX:PermSize=256m -XX:NewRatio=5 -XX:CMSTriggerRatio=50 ' +
                   ' -- --no-splash -macro fiji_z_align2_488wga.py "' +
                   exp_folder_fixed + '"')
        subprocess.check_call(fijisub ,shell=True)
        while not os.path.exists(exp_folder + "elastic_align\\for_align\\" +
                                 "%03d" % (slicenum-1) + ".png"):
            print ("...elastic alignment processing...")
            time.sleep(120) 

        
    #if elastic align ran, but images didn't finish saving out, finish saving images
    if os.path.exists(exp_folder + "elastic_align"):
        print ('checking for missing images and resaving')
        cropfiles = os.listdir(exp_folder + "cropped\\storm_merged_ds\\")
        elasticfiles = os.listdir(exp_folder + "elastic_align\\storm_merged_ds\\")
        if len(cropfiles)>len(elasticfiles):
            fijisub = ('C:/Users/Vatsal/Fiji.app/ImageJ-win64 ' +
                       '-Xms50g -Xmx50g -Xincgc -XX:MaxPermSize=256m ' +
                       '-XX:PermSize=256m -XX:NewRatio=5 -XX:CMSTriggerRatio=50 ' +
                       ' -- --no-splash -macro fiji_z_align3.py "' + exp_folder_fixed + '"')
            subprocess.check_call(fijisub ,shell=True)
            while not os.path.exists(exp_folder + "elastic_align\\for_align\\" +
                                 "%03d" % (slicenum-1) + ".png"):
                print ("...elastic images saving out...")
                time.sleep(60)
                
if __name__ == "__main__":
    exp_folder = "C:/Users/Vatsal/QT_Projects/New_Pipeline/4color_test"
    z_align(exp_folder)
    #parameters to use for alignment
        #SIFT
            #Sigma-1.2
            #Steps-6
            #minSize-20
            #maxSize-400
            #Bins+Size = 8
            #ratio = 0.98
            #maxEpsilon = 30
            #InlierRatio = 0.03
            #NumInlier = 5
        #Elastic
            #Downsample = 0.1
            #Search Size = 20
            #Block Size = 400
            #Mesh = 36
            #geometric constraints 0.1 10 0.9
            #Sigma = 200
            #MaxDistance = 10
            #MaxRelative = 3
        #spring const try 0.5
            #Everything else stays the same
