import sys 
import glob
import os
import numpy as np
import multiprocessing as mp
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

from scipy.interpolate import make_interp_spline,BSpline

from scipy.stats import zscore 
from scipy.interpolate import UnivariateSpline

from PIL import Image
from imageio import imwrite

import yaml                 


def cal_hist(wgafile,num_images):
    print(wgafile)
    A = mpimg.imread(wgafile)
    hist,bins = np.histogram(A.ravel(),255,[1,255])
    return hist
    
def new_hist(varuse, wgafile, path4, path4a):
    print(wgafile)
    
    import pdb; pdb.set_trace() 
    
    hgram4 = varuse / sum(varuse)
    A = mpimg.imread(wgafile)
    hist,bins = np.histogram(A.ravel(),256,[0,255])
    hist_cum = np.cumsum(hist)
    a = np.array(hist[0])
    b = hgram4*(sum(hist)-hist[0])
    hgram4a = np.concatenate((a,b),axis=None)

    import pdb; pdb.set_trace() 

    constructed = np.array([])
    for i in range(256):
        new = np.array(i*np.ones(int(np.round(hgram4a[i]))))
        constructed = np.concatenate((constructed, new),axis=None)
    InCaseShort = 255 * np.ones(300)
    constructed = np.concatenate((constructed,InCaseShort),axis=None)
    
    new = np.array(range(256))
    for i in range(256):
        if i == 0:
            first = 0
        else:
            first = hist_cum[i-1]
        second = hist_cum[i]
        result = np.average(constructed[first:second])
        new[i] = np.round(result)

    B = A.ravel()
    C = new[B]
    C = C.reshape(6400,6400)
    C = np.uint8(C)
    #D = imresize(C,(640,640))
    D = np.array(Image.fromarray(C).resize((640,640)))
    name = wgafile.split('\\')
    name = name[-1]
    #!!!!!!
    #Important: need to double check whether the file name should start with 0 or 1. 
    imwrite(path4 + name , C)
    imwrite(path4a + name , D)

    return 0

    
def wga_norm_and_thresh(exp_folder, alignment_channel): 

    # Assert alignment_channel is correct 
    assert alignment_channel in [561, 488]
    
    print('Loading in analysis config') 
    
    # Read in parameters from yaml file 
    with open('./configs/bead_analysis_params.yml') as f:
            config = yaml.load(f)

    shape = (config['shape_h'], config['shape_w']) 
    
    exp_folder = os.path.normpath(exp_folder) + "\\"
    storm_merged_path = exp_folder + 'unaligned\\storm_merged\\'
    conv_align_path = exp_folder + 'unaligned\\conv_{}\\'.format(str(alignment_channel))

    storm_merged_files = glob.glob(storm_merged_path + '*.tif')
    num_merged_images = len(storm_merged_files) 
    
    wga_files = glob.glob(conv_align_path + '*.tif')
    num_wga_images = len(wga_files) 
    
    assert num_merged_images == num_wga_images, "Number of images must match!"
    num_images = num_merged_images
    
    hy3c = np.zeros((num_images, 255))
    hy4c = np.zeros((num_images, 255))

    hy3cb = np.zeros((num_images, 255))
    hy4cb = np.zeros((num_images, 255))
    
    import pdb; pdb.set_trace() 
    
    print('Calculating histograms!') 
    
    print(num_images)
    for i in range(num_images): 
        hy3c[i] = cal_hist(storm_merged_files[i], num_images) # storm_merged 
        hy4c[i] = cal_hist(wga_files[i], num_images) # conv_561 
    
    # Normalizing counts to 0-1 range 
    hy3cb = hy3c / hy3c.sum(axis=1, keepdims=True)        
    hy4cb = hy4c / hy4c.sum(axis=1, keepdims=True)   

    chan = hy4cb 
    varuse4 = np.zeros([num_images, 255])
    
    x_hist = np.arange(1,255) 
    x_sections = np.arange(0, num_images)
    
    print('Thresholding!!') 

    import pdb; pdb.set_trace() 

    for i in range(255): 
        zthresh = 3 
        curr_param = chan[:, i] # Distribution of channel i values across all images 
        
        mean = np.mean(curr_param, axis=0)
        sd = np.std(curr_param, axis=0)        
        distance_from_mean = abs(chan[:, i] - mean)
        mask = distance_from_mean < zthresh * sd
        
        # Select which sections can be used for smooth interpolation 
        currfitx = x_sections[mask]
        currfity = curr_param[mask] 
                
        # currfitx = (currfitx - np.mean(currfitx)) / (np.std(currfitx) + 0.00001)
        # currfity = (currfity - np.mean(currfity)) / (np.std(currfity) + 0.00001)
      
        spl = UnivariateSpline(currfitx, currfity)
        spl.set_smoothing_factor(0.9) 
    
        varuse4[:, i] = spl(np.arange(0,num_images))
            
    path4 = exp_folder + 'unaligned\\for_align\\'
    path4a = exp_folder + 'unaligned\\for_align_ds\\'
    
    
    print('Saving out new images!') 
    
    if not os.path.exists(path4):
        os.mkdir(path4)
    if not os.path.exists(path4a): 
        os.mkdir(path4a)
    
    hgram4 = varuse4 / varuse4.sum(axis=0, keepdims=True) # Normalize over the channels for each image 
    
    for i in range(num_images): 
        new_hist(varuse4[i,:],wga_files[i],path4,path4a)    
        i = i+1

    print('Done!') 
    return True
    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        