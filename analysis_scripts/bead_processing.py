
import glob
import os
import sys
from PIL import Image, ImageFilter, ImageOps
import sa_library.datareader as daxspereader
import numpy
import subprocess
from scipy import misc, ndimage
import re
import sa_library.arraytoimage as arraytoimage
import sa_library.i3togrid as i3togrid
import math
import time

import yaml 

# where are your data?


def run_bead_processing(expfolder): 

    #expfolder = "Z:\\Chenghang\\chenghaz007_1.31.2019\\"
    fijifolder = "C:\\Users\\Chenghang\\Fiji.app\\macros\\"

    acq_folder = os.path.normpath(expfolder) + "\\acquisition\\"

    # Read in parameters from yaml file 
    with open('./configs/bead_analysis_params.yml') as f:
            config = yaml.load(f)

    shape = (config['shape_h'], config['shape_w'])
    scaled_shape = (config['shape_h'] * config['scale'], 
                    config['shape_w'] * config['scale'])

    frame_750 = config['frame_750']
    frame_647 = config['frame_647']
    frame_561 = config['frame_561']
    frame_488 = config['frame_488']


    print(acq_folder, acq_folder + "distcorr\\")
    if not os.path.exists(acq_folder + "avg561ffc.tif"):

        if not os.path.isdir(acq_folder + "distcorr\\"):
            os.mkdir(acq_folder + "distcorr\\")

        # determine the 99th percentile of intensity values for all pixels in the conventional images

        # find all the Visconv movies
        Visconv_files = glob.glob(acq_folder + 'Conv_' + '*.dax')
        
        if len(Visconv_files)>0:
            cnt = 0
        # pad matrices
            aperc_v488 = [0]*len(Visconv_files)
            aperc_v561 = [0]*len(Visconv_files)
            aperc_v647 = [0]*len(Visconv_files)
            for file in Visconv_files:
        # read 647 image intensities and load to matrix
                dax_file = daxspereader.inferReader(file)
                image = dax_file.loadAFrame(frame_647).astype(numpy.uint16)
                aperc_v647[cnt] = numpy.percentile(image, 99.999)
        # read 561 image intensities and load to matrix
                image = dax_file.loadAFrame(frame_561).astype(numpy.uint16)
                aperc_v561[cnt] = numpy.percentile(image, 99.999)
        # read 488 image intensities and load to matrix
                image = dax_file.loadAFrame(frame_488).astype(numpy.uint16)
                aperc_v488[cnt] = numpy.percentile(image, 99.999)
               
                cnt = cnt+1

        # find all the IRconv movies
        IRconv_files = glob.glob(acq_folder + 'Conv_' + '*.dax')

        assert len(Visconv_files) + len(IRconv_files) > 0, "Files are missing! (Visconv or IRConv)" 

        if len(IRconv_files)>0:
            cnt = 0
        # pad matrices
            aperc_IR750 = [0]*len(IRconv_files)
            aperc_IR647 = [0]*len(IRconv_files)
            for file in IRconv_files:
        # read 750 image intensities and load to matrix
                dax_file = daxspereader.inferReader(file)
                image = dax_file.loadAFrame(frame_750).astype(numpy.uint16)
                aperc_IR750[cnt] = numpy.percentile(image, 99.999)
        # read 647 image intensities and load to matrix
                image = dax_file.loadAFrame(frame_647).astype(numpy.uint16)
                aperc_IR647[cnt] = numpy.percentile(image, 99.999)
                
                cnt = cnt+1

        # compute the mean 99th percentile for all images and save in a list          
        rel_conv_ints = [0]*5
        rel_conv_ints[0] =  numpy.mean(aperc_v488)/256
        rel_conv_ints[1] =  numpy.mean(aperc_v561)/256
        rel_conv_ints[2] =  numpy.mean(aperc_v647)/256
        rel_conv_ints[3] =  numpy.mean(aperc_IR647)/256
        rel_conv_ints[4] =  numpy.mean(aperc_IR750)/256
        print (rel_conv_ints,"are the relative conventional intensities")

        # determine the 99th percentile of intensity values for all pixels in the ffc images

        # find all the VisFFC movies
        VisFFC_files = glob.glob(acq_folder + 'VisFFC_' + '*.dax')
        print(VisFFC_files)
        if len(VisFFC_files)>0:
            cnt = 0
        # pad matrices
            aperc_v488 = [0]*len(VisFFC_files)
            aperc_v561 = [0]*len(VisFFC_files)
            aperc_v647 = [0]*len(VisFFC_files)
            for file in VisFFC_files:
        # read 647 image intensities and load to matrix
                dax_file = daxspereader.DaxReader(file)
                image = dax_file.loadAFrame(frame_647).astype(numpy.uint16)
                aperc_v647[cnt] = numpy.percentile(image, 99.95)
        # read 561 image intensities and load to matrix
                image = dax_file.loadAFrame(frame_561).astype(numpy.uint16)
                aperc_v561[cnt] = numpy.percentile(image, 99.95)
        # read 488 image intensities and load to matrix
                image = dax_file.loadAFrame(frame_488).astype(numpy.uint16)
                aperc_v488[cnt] = numpy.percentile(image, 99.95)

                cnt = cnt+1

        # find all the IRFFC movies
        IRFFC_files = glob.glob(acq_folder + 'FFC_' + '*.dax')
        if len(IRFFC_files)>0:
            cnt = 0
        # pad matrices
            aperc_IR750 = [0]*len(IRFFC_files)
            aperc_IR647 = [0]*len(IRFFC_files)
            for file in IRFFC_files:
        # read 750 image intensities and load to matrix (adjusted to exclude small fraction of saturated pixels)
                dax_file = daxspereader.DaxReader(file)
                image = dax_file.loadAFrame(frame_750).astype(numpy.uint16)
                aperc_IR750[cnt] = numpy.percentile(image, 99.95)
        # read 647 image intensities and load to matrix (adjusted to exclude small fraction of saturated pixels)
                image = dax_file.loadAFrame(frame_647).astype(numpy.uint16) # FRAME 39 
                aperc_IR647[cnt] = numpy.percentile(image, 99.95)

                cnt = cnt+1

        # compute the mean 99th percentile for all images and save in a list 
        rel_ffc_ints = [0]*5
        rel_ffc_ints[0] =  numpy.mean(aperc_v488)/256
        rel_ffc_ints[1] =  numpy.mean(aperc_v561)/256
        rel_ffc_ints[2] =  numpy.mean(aperc_v647)/256
        rel_ffc_ints[3] =  numpy.mean(aperc_IR647)/256
        rel_ffc_ints[4] =  numpy.mean(aperc_IR750)/256
        print (rel_ffc_ints,"are the relative FFC intensities")

        # locate, correct the intensity, and save the flat field correction images

        # find all the VisFFC movies
        VisFFC_files = glob.glob(acq_folder + "FFC*.dax")
        if len(VisFFC_files)>0:
            print ("saving out VisFFC files")
            for file in VisFFC_files:
                print ("File:", os.path.basename(file))
        # load 488 FFC images 
                dax_file = daxspereader.inferReader(file)
                image = dax_file.loadAFrame(frame_488).astype(numpy.uint16) # FRAME 79
        # normalize histogram
                image = numpy.floor_divide(image,int(rel_ffc_ints[0]))
                #print (int(rel_ffc_ints[0]), " is the 488 ffc mean ")
        # generate image and convert to grayscale
                pilimage = Image.fromarray(image,'I;16')
                pilimage = pilimage.convert('L')
               # pilimage = pilimage.rotate(-90)
               # pilimage = pilimage.transpose(Image.FLIP_LEFT_RIGHT)
        # save the result
                name = os.path.basename(file)
                pilimage.save(acq_folder + "488" + name[:-4] + ".tif")
                
        # load 561 FFC images 
                image = dax_file.loadAFrame(9).astype(numpy.uint16) # FRAME 59
        # normalize histogram
                image = numpy.floor_divide(image,int(rel_ffc_ints[1]))
        # generate image and convert to grayscale
                pilimage = Image.fromarray(image,'I;16')
                pilimage = pilimage.convert('L')
               # pilimage = pilimage.rotate(-90)
                #pilimage = pilimage.transpose(Image.FLIP_LEFT_RIGHT)
        # save the result
                name = os.path.basename(file)
                pilimage.save(acq_folder + "561" + name[:-4] + ".tif")

        # load 647 FFC images 
                image = dax_file.loadAFrame(4).astype(numpy.uint16) # FRAME 39
        # normalize histogram
                image = numpy.floor_divide(image,int(rel_ffc_ints[2]))
        # generate image and convert to grayscale
                pilimage = Image.fromarray(image,'I;16')
                pilimage = pilimage.convert('L')
                #pilimage = pilimage.rotate(-90)
                #pilimage = pilimage.transpose(Image.FLIP_LEFT_RIGHT)
        # save the result
                name = os.path.basename(file)
                pilimage.save(acq_folder + "647" + name[:-4] + ".tif")

            # generate average FFC images from stack
            shape = (896, 896, 1) 
            
        # set FFC image range
            for i in range(9):
        # load image data as array and resize
                im = Image.open(acq_folder + '488FFC_' + str(i) + '.tif')
                imnp = numpy.array(im)
                print(imnp.shape)
                imnp = numpy.reshape(imnp,shape)
                if i == 0:
                    imstack = imnp
                else:  
                    imstack = numpy.concatenate((imstack, imnp), axis=2)
        # average images 
            avgim = numpy.average(imstack,axis=2)
            pilimage = Image.fromarray(avgim)
        # blur image 
            ffc488np = ndimage.gaussian_filter(pilimage,60)
            ffc488np[ffc488np == 0] = 1
            ffc488mean = numpy.mean(ffc488np)
            print (ffc488mean)
            print (" is the 488 ffc mean")
        # convert to 8bit grayscale and save
            pilimage = Image.fromarray(ffc488np)
            # pilimage = pilimage.convert('L')
            pilimage.save(acq_folder + "avg488ffc.tif")

        # set FFC image range     
            for i in range(9):
                # load image data as array and resize
                im = Image.open(acq_folder + '561FFC_' + str(i) + '.tif')
                imnp = numpy.array(im)
                imnp = numpy.reshape(imnp,shape)
                if i == 0:
                    imstack = imnp
                else:  
                    imstack = numpy.concatenate((imstack, imnp), axis=2)
        # average images 
            avgim = numpy.average(imstack,axis=2)
            pilimage = Image.fromarray(avgim)
        # blur image 
            ffc561np = ndimage.gaussian_filter(pilimage,60)
            ffc561np[ffc561np == 0] = 1
            ffc561mean = numpy.mean(ffc561np)
            print (ffc561mean)
            print (" is the 561 ffc mean")
        # convert to 8bit grayscale and save
            pilimage = Image.fromarray(ffc561np)
            #pilimage = pilimage.convert('L')
            pilimage.save(acq_folder + "avg561ffc.tif")

        # set FFC image range     
            for i in range(9):
        # load image data as array and resize
                im = Image.open(acq_folder + '647FFC_' + str(i) + '.tif')
                imnp = numpy.array(im)
                imnp = numpy.reshape(imnp,shape)
                if i == 0:
                    imstack = imnp
                else:  
                    imstack = numpy.concatenate((imstack, imnp), axis=2)
        # average images 
            avgim = numpy.average(imstack,axis=2)
            pilimage = Image.fromarray(avgim)
        # blur image 
            ffcVis647np = ndimage.gaussian_filter(pilimage,60)
            ffcVis647np[ffcVis647np == 0] = 1
            ffcVis647mean = numpy.mean(ffcVis647np)
            print (ffcVis647mean)
            print (" is the 647vis ffc mean")
        # convert to 8bit grayscale and save
            pilimage = Image.fromarray(ffcVis647np)
            pilimage = pilimage.convert('L')
            pilimage.save(acq_folder + "avgVis647ffc.tif")

            
        # find all the IRFFC movies
        IRFFC_files = glob.glob(acq_folder + "FFC*.dax")
        if len(IRFFC_files)>0:
            for file in IRFFC_files:
        # load 750 FFC images 
                dax_file = daxspereader.inferReader(file)
                image = dax_file.loadAFrame(1).astype(numpy.uint16)
        # normalize histogram
                image = numpy.floor_divide(image,int(rel_ffc_ints[4]))
        # generate image and convert to 8bit grayscale
                pilimage = Image.fromarray(image,'I;16')
                pilimage = pilimage.convert('L')
              #  pilimage = pilimage.rotate(-90)
               # pilimage = pilimage.transpose(Image.FLIP_LEFT_RIGHT)
        # save the result
                name = os.path.basename(file)
                pilimage.save(acq_folder + "750" + name[:-4] + ".tif")

        # set FFC image range    
            for i in range(9):
        # load image data as array and resize
                im = Image.open(acq_folder + '647FFC_' + str(i) + '.tif')
                imnp = numpy.array(im)
                imnp = numpy.reshape(imnp,shape)
                if i == 0:
                    imstack = imnp
                else:  
                    imstack = numpy.concatenate((imstack, imnp), axis=2)
        # average images 
            avgim = numpy.average(imstack,axis=2)
            pilimage = Image.fromarray(avgim)
        # blur image 
            ffcIR647np = ndimage.gaussian_filter(pilimage,60)
            ffcIR647np[ffcIR647np == 0] = 1
            ffcIR647mean = numpy.mean(ffcIR647np)
            print (ffcIR647mean)
            print (" is the 647IR ffc mean")
        # save
            pilimage = Image.fromarray(ffcIR647np)
            pilimage = pilimage.convert('L')
            pilimage.save(acq_folder + "avgIR647ffc.tif")

        # set FFC image range    
            for i in range(9):
        # load image data as array and resize
                im = Image.open(acq_folder + '750FFC_' + str(i) + '.tif')
                imnp = numpy.array(im)
                imnp = numpy.reshape(imnp,shape)
                if i == 0:
                    imstack = imnp
                else:  
                    imstack = numpy.concatenate((imstack, imnp), axis=2)
        # average images 
            avgim = numpy.average(imstack,axis=2)
            pilimage = Image.fromarray(avgim)
        # blur image 
            ffc750np = ndimage.gaussian_filter(pilimage, 60)
            ffc750np[ffc750np == 0] = 1
            ffc750mean = numpy.mean(ffc750np)
            print (ffc750mean)
            print (" is the 750 ffc mean")
        # save
            pilimage = Image.fromarray(ffc750np)
            pilimage = pilimage.convert('L')
            pilimage.save(acq_folder + "avg750ffc.tif")

        # save out images for distortion correction generation

        # find the IRFFC dax movies
        files = glob.glob(acq_folder + "IRffc*.dax")
        for i in range(9):
            if len(files)>0:
                im = Image.open(acq_folder + '647ffc_' + str(i) + '.tif')
                imnp = numpy.array(im)*ffcIR647mean
                corr = numpy.array(imnp/ffcIR647np)
            else:
                im = Image.open(acq_folder + '647ffc_' + str(i) + '.tif')            
                imnp = numpy.array(im)*ffcVis647mean
                corr = numpy.array(imnp/ffcVis647np)
            pilimage = Image.fromarray(corr)
            pilimage = pilimage.resize(scaled_shape,Image.BILINEAR)
            pilimage = pilimage.convert('L')
            pilimage = ImageOps.autocontrast(pilimage,cutoff=0)
            pilimage.save(acq_folder + 'distcorr/647ffc_' + "%02d" % i + '.tif')

    # perform lens distortion correction

    # pass argument to FIJI (inelegant, but functional)
    #arg2 = fijifolder

    # initialize lens distortion correction
    #print ("starting distortion correction generation")
    #if not os.path.isfile(acq_folder + "dist_corr/distCorr.txt"):
        # set cmd line arguments and open subprocess
        #fijisub = ('C:\\Users\\Colenso\\Fiji.app\\ImageJ-win64' + ' -macro fiji_distcorr_make.py ' + arg2)
        #print (fijisub)
        #subprocess.Popen(fijisub, shell=True)
        #print('popen is running')

    # look for distortion correction file to finish process
    #while not os.path.isfile(acq_folder + "dist_corr/distCorr.txt"):
        #time.sleep(10)
    print ("bead images saved")  
    return "success" 


