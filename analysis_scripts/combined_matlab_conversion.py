import numpy as np
import sys
import os
import pathlib
import sys
import h5py
from PIL import Image

from .corr_mols_v8 import corr_mols

import skimage
import skimage.transform
from skimage.transform import (estimate_transform)

import matplotlib.pyplot as plt
import yaml 

import tifffile as timg;
from tifffile import imsave;
import skimage as skimg;
import scipy
import cv2
from PIL import Image

from .py_dftregistration import dftregistration;

from .imresize import imresize;
from .imresize import imadjust;

def ecdf(data):
    #compute ECDF
    x = np.sort(data)
    n = x.size
    y = np.arange(1, n+1) / n
    return(x,y)

'''
%function f = gen_bead_warp(arg1)
% this script is used to generate a bead map and corresponding warping 
% transforms based on bead images collected through the emission filter
% wheel. IRbeads are visible 750/647 channels (100nm 715/755 beads) and 
% Visbeads are visible in 647/561/488 channels (200nm beads).
% Transforms map 750/561/488 to the 647 channel.
% Updated 06_10_18 for analysis of .hdf5 files
% Colenso Speer
% make sure to change frames in hdf5read to properly grab beads images
'''

##local_exp = 'Y:\\Chenghang\\06_Testing\\GUI_test\\'

def gen_bead_warp_and_chrom_align(local_exp, alignment_channel):


    acq_path = local_exp + "acquisition\\"
    base_pth = acq_path + "bins\\"
    save_pth = acq_path + "bead_registration\\"

     
    # Read in parameters from yaml file 
    with open('./configs/bead_analysis_params.yml') as f:
        config = yaml.load(f)

    frame_750 = config['frame_750']
    frame_647 = config['frame_647']
    frame_561 = config['frame_561']
    frame_488 = config['frame_488']

    #% set scale for image expansion
    scale = config['scale']
     
    #% size of field in pixels
    xdiv = config['shape_w']*scale
    ydiv = config['shape_h']*scale

    #% pixel size in final images
    nm_per_pix = 155/scale

    #% axis/label/title fontsize for plots
    ax_fs = 10 
    la_fs = 10
    ti_fs = 10


    if not os.path.exists(save_pth):
        os.mkdir(save_pth)


    fileExt = r"*.dax"
    numbeadmov=len(list(pathlib.Path(acq_path).glob(fileExt)))

    if numbeadmov <=1:
        sys.exit()

    track_file_nms1 = list(pathlib.Path(base_pth).glob(r"Regbead_1_*750mlist.hdf5"))
    track_file_nms2 = list(pathlib.Path(base_pth).glob(r"Regbead_1_*IR647mlist.hdf5"))
    track_file_nms3 = list(pathlib.Path(base_pth).glob(r"Regbead_1_*Vis647mlist.hdf5"))
    track_file_nms4 = list(pathlib.Path(base_pth).glob(r"Regbead_1_*561mlist.hdf5"))
    track_file_nms5 = list(pathlib.Path(base_pth).glob(r"Regbead_1_*488mlist.hdf5"))

        
    for p in range(1, 5):
        if p == 1:
            length = 1;
            #% radius tolerance for matching left/right channels
            match_radius = 1*scale; 
        elif p==2:
            length = 6;
            #% radius tolerance for matching left/right channels
            match_radius = 0.5*scale; 
        else:
            length = len(track_file_nms1);
            #% radius tolerance for matching left/right channels
            match_radius = 0.1*scale; 

        list_S750_1_x = []
        list_S750_1_y = []
            
        list_S750_2_x = []
        list_S750_2_y = []
            
        list_S561_1_x = []
        list_S561_1_y = []
            
        list_S561_2_x = []
        list_S561_2_y = []
            
        list_S488_1_x = []
        list_S488_1_y = []
            
        list_S488_2_x = []
        list_S488_2_y = []
            
        for k in range(0, length):
            #% load in 750 mlist 
            track_file_nm = track_file_nms1[k]
            print(track_file_nm);
        
            #% find molecules identified in frame 1
            hf = h5py.File(track_file_nm, 'r')
            X_value = hf.get('/fr_1/x')
            Y_value = hf.get('/fr_1/y')
            
            #% resize x and y values for larger canvas
            M2_1_x = (np.array(X_value))*scale
            M2_1_y = (np.array(Y_value))*scale     
        
            track_file_nm = track_file_nms2[k]
            print(track_file_nm);
        
            #% find molecules identified in frame 43
            hf = h5py.File(track_file_nm, 'r')
            X_value = hf.get('/fr_7/x')
            Y_value = hf.get('/fr_7/y')
        
            #% resize x and y values for larger canvas
            M2_2_x = (np.array(X_value))*scale
            M2_2_y = (np.array(Y_value))*scale      
        
            track_file_nm = track_file_nms3[k]
            print(track_file_nm);
        
            #% find molecules identified in frame 61
            hf = h5py.File(track_file_nm, 'r')
            X_value = hf.get('/fr_10/x')
            Y_value = hf.get('/fr_10/y')
        
            #% resize x and y values for larger canvas
            M2_3_x = (np.array(X_value))*scale
            M2_3_y = (np.array(Y_value))*scale      
             
            track_file_nm = track_file_nms4[k]
            print(track_file_nm);
        
            #% find molecules identified in frame 83
            hf = h5py.File(track_file_nm, 'r')
            X_value = hf.get('/fr_14/x')
            Y_value = hf.get('/fr_14/y')
        
            #% resize x and y values for larger canvas
            M2_4_x = (np.array(X_value))*scale
            M2_4_y = (np.array(Y_value))*scale        
        
            track_file_nm = track_file_nms5[k]
            print(track_file_nm);
        
            #% find molecules identified in frame 123
            hf = h5py.File(track_file_nm, 'r')
            X_value = hf.get('/fr_21/x')
            Y_value = hf.get('/fr_21/y')
        
            #% resize x and y values for larger canvas
            M2_5_x = (np.array(X_value))*scale
            M2_5_y = (np.array(Y_value))*scale       
        
            set488_inds, = np.where((M2_5_x>0) & (M2_5_x<=xdiv) & (M2_5_y<=ydiv));
            set488_pos_x = (M2_5_x)[set488_inds];
            set488_pos_y = (M2_5_y)[set488_inds];
            
            set561_inds, = np.where((M2_4_x>0) & (M2_4_x<=xdiv) & (M2_4_y<=ydiv));
            set561_pos_x = (M2_4_x)[set561_inds];
            set561_pos_y = (M2_4_y)[set561_inds];
        
            setV647_inds, = np.where((M2_3_x>0) & (M2_3_x<=xdiv) & (M2_3_y<=ydiv));
            setV647_pos_x = (M2_3_x)[setV647_inds];
            setV647_pos_y = (M2_3_y)[setV647_inds];
        
            set647_inds, = np.where((M2_2_x>0) & (M2_2_x<=xdiv) & (M2_2_y<=ydiv));
            set647_pos_x = (M2_2_x)[set647_inds];
            set647_pos_y = (M2_2_y)[set647_inds];
        
            set750_inds, = np.where((M2_1_x>0) & (M2_1_x<=xdiv) & (M2_1_y<=ydiv));
            set750_pos_x = (M2_1_x)[set750_inds];
            set750_pos_y = (M2_1_y)[set750_inds];
        
            #% make new transform variables
            if p == 1:
                tform_750start = estimate_transform('similarity', np.array([[1,0],[0,1]]), np.array([[1,0],[0,1]]))
                tform_561start = estimate_transform('similarity', np.array([[1,0],[0,1]]), np.array([[1,0],[0,1]]))
                tform_488start = estimate_transform('similarity', np.array([[1,0],[0,1]]), np.array([[1,0],[0,1]]))
                 
                tform_750_2_647 = tform_750start
                tform_561_2_647 = tform_561start
                tform_488_2_647 = tform_488start           
            else:
                tform_750start=tform_750_2_647
                tform_561start=tform_561_2_647
                tform_488start=tform_488_2_647
         
            # Call corr_mols function
            matched750_set1_inds,matched750_set2_inds,unmatched750_set1_inds,unmatched750_set2_inds = corr_mols(set750_pos_x,set750_pos_y,set647_pos_x,set647_pos_y,tform_750start, match_radius);              
            matched561_set1_inds,matched561_set2_inds,unmatched561_set1_inds,unmatched561_set2_inds = corr_mols(set561_pos_x,set561_pos_y,setV647_pos_x,setV647_pos_y,tform_561start, match_radius);      
            matched488_set1_inds,matched488_set2_inds,unmatched488_set1_inds,unmatched488_set2_inds = corr_mols(set488_pos_x,set488_pos_y,setV647_pos_x,setV647_pos_y,tform_488start, match_radius);  
            
            # % display the bead matching efficacy
            print(str(len(matched750_set1_inds)) + "/" + str((2*len(matched750_set1_inds)+len(unmatched750_set1_inds)+len(unmatched750_set2_inds))/2) + " 750 molecules matched")
            print(str(len(matched561_set1_inds)) + "/" + str((2*len(matched561_set1_inds)+len(unmatched561_set1_inds)+len(unmatched561_set2_inds))/2) + " 561 molecules matched")
            print(str(len(matched488_set1_inds)) + "/" + str((2*len(matched488_set1_inds)+len(unmatched488_set1_inds)+len(unmatched488_set2_inds))/2) + " 488 molecules matched")
            print("bead image number: " + str(k+1))
              
            #% return matched molecules, set1
            S750_1_x = (M2_1_x)[(set750_inds)[matched750_set1_inds]]; 
            S750_1_y = (M2_1_y)[(set750_inds)[matched750_set1_inds]];             
            #% return matched molecules, set2
            S750_2_x = (M2_2_x)[(set647_inds)[matched750_set2_inds]]; 
            S750_2_y = (M2_2_y)[(set647_inds)[matched750_set2_inds]];
            
            #% return matched molecules, set1
            S561_1_x = (M2_4_x)[(set561_inds)[matched561_set1_inds]]; 
            S561_1_y = (M2_4_y)[(set561_inds)[matched561_set1_inds]];             
            #% return matched molecules, set2
            S561_2_x = (M2_3_x)[(setV647_inds)[matched561_set2_inds]]; 
            S561_2_y = (M2_3_y)[(setV647_inds)[matched561_set2_inds]]; 
        
            #% return matched molecules, set1
            S488_1_x = (M2_5_x)[(set488_inds)[matched488_set1_inds]]; 
            S488_1_y = (M2_5_y)[(set488_inds)[matched488_set1_inds]]; 
            #% return matched molecules, set2
            S488_2_x = (M2_3_x)[(setV647_inds)[matched488_set2_inds]]; 
            S488_2_y = (M2_3_y)[(setV647_inds)[matched488_set2_inds]]; 
            
            #adding matched molecules to list
            list_S750_1_x.append(S750_1_x)
            list_S750_1_y.append(S750_1_y)
            
            list_S750_2_x.append(S750_2_x)
            list_S750_2_y.append(S750_2_y)
            
            list_S561_1_x.append(S561_1_x)
            list_S561_1_y.append(S561_1_y)
            
            list_S561_2_x.append(S561_2_x)
            list_S561_2_y.append(S561_2_y)
            
            list_S488_1_x.append(S488_1_x)
            list_S488_1_y.append(S488_1_y)
            
            list_S488_2_x.append(S488_2_x)
            list_S488_2_y.append(S488_2_y)
           
                
        comb_set750_1_pos_x = np.array([]); #% set1 = left channel (pixel units)
        comb_set750_1_pos_y = np.array([]); 
        comb_set750_2_pos_x = np.array([]); #% set2 = right channel (pixel units)
        comb_set750_2_pos_y = np.array([]);
        
        comb_set561_1_pos_x = np.array([]); #% set1 = left channel (pixel units)
        comb_set561_1_pos_y = np.array([]); 
        comb_set561_2_pos_x = np.array([]); #% set2 = right channel (pixel units)
        comb_set561_2_pos_y = np.array([]);
        
        comb_set488_1_pos_x = np.array([]); #% set1 = left channel (pixel units)
        comb_set488_1_pos_y = np.array([]); 
        comb_set488_2_pos_x = np.array([]); #% set2 = right channel (pixel units)
        comb_set488_2_pos_y = np.array([]); 
           
        for k in range(0, length):
            try:
                comb_set750_1_pos_x = np.concatenate([comb_set750_1_pos_x, list_S750_1_x[k]]);
                comb_set750_1_pos_y = np.concatenate([comb_set750_1_pos_y, list_S750_1_y[k]]);
                comb_set750_2_pos_x = np.concatenate([comb_set750_2_pos_x, list_S750_2_x[k]]);
                comb_set750_2_pos_y = np.concatenate([comb_set750_2_pos_y, list_S750_2_y[k]]);
            except:
                print("Something went wrong")
                
        for k in range(0, length):
            try:
                comb_set561_1_pos_x = np.concatenate([comb_set561_1_pos_x, list_S561_1_x[k]]);
                comb_set561_1_pos_y = np.concatenate([comb_set561_1_pos_y, list_S561_1_y[k]]);
                comb_set561_2_pos_x = np.concatenate([comb_set561_2_pos_x, list_S561_2_x[k]]);
                comb_set561_2_pos_y = np.concatenate([comb_set561_2_pos_y, list_S561_2_y[k]]);
            except:
                print("Something went wrong")
          
        for k in range(0, length):
            try:
                comb_set488_1_pos_x = np.concatenate([comb_set488_1_pos_x, list_S488_1_x[k]]);
                comb_set488_1_pos_y = np.concatenate([comb_set488_1_pos_y, list_S488_1_y[k]]);
                comb_set488_2_pos_x = np.concatenate([comb_set488_2_pos_x, list_S488_2_x[k]]);
                comb_set488_2_pos_y = np.concatenate([comb_set488_2_pos_y, list_S488_2_y[k]]);
            except:
                print("Something went wrong")
          
        #% find mapping from set2 (right) onto set1 (left), save tform_right2left
        set750_1_points = np.column_stack([comb_set750_1_pos_x, comb_set750_1_pos_y]);
        set750_2_points = np.column_stack([comb_set750_2_pos_x, comb_set750_2_pos_y]);
        
        set561_1_points = np.column_stack([comb_set561_1_pos_x, comb_set561_1_pos_y]);
        set561_2_points = np.column_stack([comb_set561_2_pos_x, comb_set561_2_pos_y]);
        
        set488_1_points = np.column_stack([comb_set488_1_pos_x, comb_set488_1_pos_y]);
        set488_2_points = np.column_stack([comb_set488_2_pos_x, comb_set488_2_pos_y]);
        
        # Only do "affine" transformations for now
        tform_750_2_647 = estimate_transform('similarity', set750_1_points, set750_2_points);     
        tform_561_2_647 = estimate_transform('similarity', set561_1_points, set561_2_points);   
        tform_488_2_647 = estimate_transform('similarity', set488_1_points, set488_2_points); 
         
        warped_set750_2_pos_x = (tform_750_2_647.inverse(set750_2_points))[:,0];
        warped_set750_2_pos_y = (tform_750_2_647.inverse(set750_2_points))[:,1];
        
        warped_set561_2_pos_x = (tform_561_2_647.inverse(set561_2_points))[:,0];
        warped_set561_2_pos_y = (tform_561_2_647.inverse(set561_2_points))[:,1];         
            
        warped_set488_2_pos_x = (tform_488_2_647.inverse(set488_2_points))[:,0];
        warped_set488_2_pos_y = (tform_488_2_647.inverse(set488_2_points))[:,1];  
        
        #get errors
        set750_2_warp_error_x = comb_set750_1_pos_x - warped_set750_2_pos_x;
        set750_2_warp_error_y = comb_set750_1_pos_y - warped_set750_2_pos_y;
        std_set750_2_warp_error_x = np.std(set750_2_warp_error_x);
        std_set750_2_warp_error_y = np.std(set750_2_warp_error_y);
        std_total750_warp_error = np.std(np.sqrt(np.square(set750_2_warp_error_x) + np.square(set750_2_warp_error_y)));
        
        set750_2_orig_error_x = comb_set750_1_pos_x - comb_set750_2_pos_x ;
        set750_2_orig_error_y = comb_set750_1_pos_y - comb_set750_2_pos_y;
        std_set750_2_orig_error_x = np.std(set750_2_orig_error_x);
        std_set750_2_orig_error_y = np.std(set750_2_orig_error_y);
        
        set561_2_warp_error_x = comb_set561_1_pos_x - warped_set561_2_pos_x;
        set561_2_warp_error_y = comb_set561_1_pos_y - warped_set561_2_pos_y;
        std_set561_2_warp_error_x = np.std(set561_2_warp_error_x);
        std_set561_2_warp_error_y = np.std(set561_2_warp_error_y);
        std_total561_warp_error = np.std(np.sqrt(np.square(set561_2_warp_error_x) + np.square(set561_2_warp_error_y)));
        
        set561_2_orig_error_x = comb_set561_1_pos_x - comb_set561_2_pos_x;
        set561_2_orig_error_y = comb_set561_1_pos_y - comb_set561_2_pos_y;
        std_set561_2_orig_error_x = np.std(set561_2_orig_error_x);
        std_set561_2_orig_error_y = np.std(set561_2_orig_error_y);
        
        set488_2_warp_error_x = comb_set488_1_pos_x - warped_set488_2_pos_x;
        set488_2_warp_error_y = comb_set488_1_pos_y - warped_set488_2_pos_y;
        std_set488_2_warp_error_x = np.std(set488_2_warp_error_x);
        std_set488_2_warp_error_y = np.std(set488_2_warp_error_y);
        std_total488_warp_error = np.std(np.sqrt(np.square(set488_2_warp_error_x) + np.square(set488_2_warp_error_y)));
        
        set488_2_orig_error_x = comb_set488_1_pos_x - comb_set488_2_pos_x;
        set488_2_orig_error_y = comb_set488_1_pos_y - comb_set488_2_pos_y;
        std_set488_2_orig_error_x = np.std(set488_2_orig_error_x);
        std_set488_2_orig_error_y = np.std(set488_2_orig_error_y);
        
        #ecdf
        cdf750_x, cdf750_y = ecdf(nm_per_pix*np.sqrt(np.square(set750_2_warp_error_x) + np.square(set750_2_warp_error_y)));
        cdf90_750 = (cdf750_x[np.where(cdf750_y>0.9)])[0];
        print("90% of 750 beads aligned to ", str(cdf90_750), "nm, using ", str(len(set750_1_points[:,1])), " beads")

        cdf561_x, cdf561_y = ecdf(nm_per_pix*np.sqrt(np.square(set561_2_warp_error_x) + np.square(set561_2_warp_error_y)));
        cdf90_561 = (cdf561_x[np.where(cdf561_y>0.9)])[0];
        print("90% of 561 beads aligned to ", str(cdf90_561), "nm, using ", str(len(set561_1_points[:,1])), " beads")

        cdf488_x, cdf488_y = ecdf(nm_per_pix*np.sqrt(np.square(set488_2_warp_error_x) + np.square(set488_2_warp_error_y)));
        cdf90_488 = (cdf488_x[np.where(cdf488_y>0.9)])[0];
        print("90% of 488 beads aligned to ", str(cdf90_488), "nm, using ", str(len(set488_1_points[:,1])), " beads")
        
        print("End of p loop: " + str(p))

    #plot deviations due to chromatic abberation
    fac = 20

    fig, axs = plt.subplots(2, 3)

    axs[0, 0].set_title('561 to 647, no warping 20x mag. residuals', fontsize = 5, fontweight='bold')
    axs[0, 0].scatter(comb_set561_2_pos_x,comb_set561_2_pos_y, color='lime', marker = ',', s = .005)
    q1 = axs[0, 0].quiver(comb_set561_2_pos_x, comb_set561_2_pos_y, fac*set561_2_orig_error_x, fac*set561_2_orig_error_y, color = 'm', width=0.002,angles='xy', scale_units='xy', scale=1)

    axs[0, 1].set_title('488 to 647, no warping 20x mag. residuals', fontsize = 5, fontweight='bold')
    axs[0, 1].scatter(comb_set488_2_pos_x,comb_set488_2_pos_y, color='lime', marker = ',', s = .005)
    q2 = axs[0, 1].quiver(comb_set488_2_pos_x, comb_set488_2_pos_y, fac*set488_2_orig_error_x, fac*set488_2_orig_error_y, color = 'm', width=0.002,angles='xy', scale_units='xy', scale=1)
      
    axs[1, 0].set_title('561 to 647, affine warp, 20x mag. residuals', fontsize = 5, fontweight='bold')
    axs[1, 0].scatter(comb_set561_2_pos_x,comb_set561_2_pos_y, color='lime', marker = ',', s = .005)
    q3 = axs[1, 0].quiver(comb_set561_2_pos_x, comb_set561_2_pos_y, fac*set561_2_warp_error_x, fac*set561_2_warp_error_y, color = 'm', width=0.002, angles='xy', scale_units='xy', scale=1)

    axs[1, 1].set_title('488 to 647, affine warp, 20x mag. residuals', fontsize = 5, fontweight='bold')
    axs[1, 1].scatter(comb_set488_2_pos_x,comb_set488_2_pos_y, color='lime', marker = ',', s = .005)
    q4 = axs[1, 1].quiver(comb_set488_2_pos_x, comb_set488_2_pos_y, fac*set488_2_warp_error_x, fac*set488_2_warp_error_y, color = 'm', width=0.002,angles='xy', scale_units='xy', scale=1)

    fac=40

    axs[0, 2].set_title('750 to 647, no warping 40x mag. residuals', fontsize = 5, fontweight='bold')
    axs[0, 2].scatter(comb_set750_2_pos_x,comb_set750_2_pos_y, color='lime', marker = ',', s = .005)
    q5 = axs[0, 2].quiver(comb_set750_2_pos_x, comb_set750_2_pos_y, fac*set750_2_orig_error_x, fac*set750_2_orig_error_y, color = 'm', width=0.002,angles='xy', scale_units='xy', scale=1)
        
    axs[1, 2].set_title('750 to 647, affine warp, 40x mag. residuals', fontsize = 5, fontweight='bold')
    axs[1, 2].scatter(comb_set750_2_pos_x,comb_set750_2_pos_y, color='lime', marker = ',', s = .005)
    q6 = axs[1, 2].quiver(comb_set750_2_pos_x, comb_set750_2_pos_y, fac*set750_2_warp_error_x, fac*set750_2_warp_error_y, color = 'm', width=0.002,angles='xy', scale_units='xy', scale=1)


    for ax in axs.flat:
        ax.set_xlabel('X axis [pixel]', fontsize = 5)
        ax.set_ylabel('Y axis [pixel]', fontsize = 5)
        ax.xaxis.set_tick_params(labelsize=5)
        ax.yaxis.set_tick_params(labelsize=5)
        ax.margins(0.02)
        ax.invert_yaxis()
        ax.set_aspect('equal', adjustable='box')
        

    plt.tight_layout()
    #save out beadmap figures
    plt.savefig(save_pth+'BeadMap_for_registration_self.png', bbox_inches='tight', dpi=1200)


    #plot total offsets for warping
    fig, (ax1, ax2, ax3) = plt.subplots(1,3)

    cdf488_x, cdf488_y = ecdf(nm_per_pix*np.sqrt(np.square(set488_2_warp_error_x) + np.square(set488_2_warp_error_y)));     
    ax1.scatter(cdf488_x,cdf488_y, color='g', s = 0.1)   

    cdf561_x, cdf561_y = ecdf(nm_per_pix*np.sqrt(np.square(set561_2_warp_error_x) + np.square(set561_2_warp_error_y)));    
    ax2.scatter(cdf561_x,cdf561_y, color='g', s = 0.1) 


    cdf750_x, cdf750_y = ecdf(nm_per_pix*np.sqrt(np.square(set750_2_warp_error_x) + np.square(set750_2_warp_error_y)));    
    ax3.scatter(cdf750_x,cdf750_y, color='g', s = 0.1)

    ax1_title = 'Cumulative distribution self-warping 488 to 647 \n STD =' + str(round(std_total488_warp_error*nm_per_pix,1)) + 'nm \n ' + '90% of 488 beads aligned to ' + str(round(cdf90_488,3)) + 'nm'
    ax2_title = 'Cumulative distribution self-warping 561 to 647 \n STD =' + str(round(std_total561_warp_error*nm_per_pix,1)) + 'nm \n ' + '90% of 561 beads aligned to ' + str(round(cdf90_561,3)) + 'nm'
    ax3_title = 'Cumulative distribution self-warping 750 to 647 \n STD =' + str(round(std_total750_warp_error*nm_per_pix,1)) + 'nm \n ' + '90% of 750 beads aligned to ' + str(round(cdf90_750,3)) + 'nm'

    ax1.set_title(ax1_title, fontsize = 4, fontweight='bold')
    ax1.set_xlabel('distance error [nm]', fontsize = 5)
    ax1.set_ylabel('cumulative probability', fontsize = 5)
    ax2.set_title(ax2_title, fontsize = 4, fontweight='bold')
    ax2.set_xlabel('distance error [nm]', fontsize = 5)
    ax2.set_ylabel('cumulative probability', fontsize = 5)
    ax3.set_title(ax3_title, fontsize = 4, fontweight='bold')
    ax3.set_xlabel('distance error [nm]]', fontsize = 5)
    ax3.set_ylabel('cumulative probability', fontsize = 5)

    ax1.xaxis.set_tick_params(labelsize=5)
    ax1.yaxis.set_tick_params(labelsize=5)
    ax2.xaxis.set_tick_params(labelsize=5)
    ax2.yaxis.set_tick_params(labelsize=5)
    ax3.xaxis.set_tick_params(labelsize=5)
    ax3.yaxis.set_tick_params(labelsize=5)


    ax1.grid(linewidth = 0.15)
    ax2.grid(linewidth = 0.15)
    ax3.grid(linewidth = 0.15)

    ax1.margins(0)
    ax2.margins(0)
    ax3.margins(0)

    fig.tight_layout()
    #save out probability distribution alignment figures
    plt.savefig(save_pth+'Cumulative_distribution_for_registration_self.png', bbox_inches='tight', dpi=1500)

    '''
    image_chrom_align_hcam starts here.
    '''

    # arg1 = "C:\\Users\\Colenso\\Desktop\\GUI_XY_test\\";
    # local_exp =  arg1;
    rel_conv_ints = '1111'; #%rel_conv_ints will be used to rescale the pixel intensities for different channels. 
    analysisfolder = local_exp + "analysis/";
    ISanalysisfolder = analysisfolder + "individual_sections/";

    '''
    %folders for merged image ouput
    mergeconv = cat(2,analysisfolder, 'unaligned/conv_merged/');
    mergestorm = cat(2, analysisfolder, 'unaligned/storm_merged/');
    conv561only = cat(2, analysisfolder, 'unaligned/conv_561/');
    mergeconv488561 = cat(2, analysisfolder, 'unaligned/conv_561_488/');
    '''

    mergeconv = analysisfolder + "unaligned/conv_merged/";
    mergestorm = analysisfolder + "unaligned/storm_merged/";

    convAlignonly = analysisfolder + "unaligned/conv_{}/".format(str(alignment_channel))

    mergeconv488561 = analysisfolder + "unaligned/conv_561_488/";
    aligned_storm = analysisfolder + "unaligned/storm_{}/".format(str(alignment_channel))

    '''
    % determine number of sections
    slices = (numel(dir(fullfile(ISanalysisfolder, '0*')))-1);
    '''

    fileExt = r"0*"
    slices=len(list(pathlib.Path(ISanalysisfolder).glob(fileExt)));
    ##slices = 1

    '''
    for slice = 0:0
        
    filename.storm488 = strcat(ISanalysisfolder,sprintf('%04d',slice),'/rawimages/for_matlab/488storm_',sprintf('%03d',slice),'.tiff');
    filename.storm561 = strcat(ISanalysisfolder,sprintf('%04d',slice),'/rawimages/for_matlab/561storm_',sprintf('%03d',slice),'.tiff');
    filename.storm647 = strcat(ISanalysisfolder,sprintf('%04d',slice),'/rawimages/for_matlab/647storm_',sprintf('%03d',slice),'.tiff');
    filename.storm750 = strcat(ISanalysisfolder,sprintf('%04d',slice),'/rawimages/for_matlab/750storm_',sprintf('%03d',slice),'.tiff');
    filename.convVis488 = strcat(ISanalysisfolder,sprintf('%04d',slice),'/rawimages/for_matlab/488Visconv_',sprintf('%03d',slice),'.tif');
    filename.convVis561 = strcat(ISanalysisfolder,sprintf('%04d',slice),'/rawimages/for_matlab/561Visconv_',sprintf('%03d',slice),'.tif');
    filename.convVis647 = strcat(ISanalysisfolder,sprintf('%04d',slice),'/rawimages/for_matlab/647Visconv_',sprintf('%03d',slice),'.tif');
    filename.convIR647 = strcat(ISanalysisfolder,sprintf('%04d',slice),'/rawimages/for_matlab/647IRconv_',sprintf('%03d',slice),'.tif');
    filename.convIR750 = strcat(ISanalysisfolder,sprintf('%04d',slice),'/rawimages/for_matlab/750IRconv_',sprintf('%03d',slice),'.tif');
    filenameout.storm488 = strcat(ISanalysisfolder,sprintf('%04d',slice),'/aligned/488storm_',sprintf('%03d',slice),'.tif');
    filenameout.storm647 = strcat(ISanalysisfolder,sprintf('%04d',slice),'/aligned/647storm_',sprintf('%03d',slice),'.tif');
    filenameout.storm561 = strcat(ISanalysisfolder,sprintf('%04d',slice),'/aligned/561storm_',sprintf('%03d',slice),'.tif');
    filenameout.storm750 = strcat(ISanalysisfolder,sprintf('%04d',slice),'/aligned/750storm_',sprintf('%03d',slice),'.tif');
    filenameout.convVis488 = strcat(ISanalysisfolder,sprintf('%04d',slice),'/aligned/488Visconv_',sprintf('%03d',slice),'.tif');
    filenameout.convVis561 = strcat(ISanalysisfolder,sprintf('%04d',slice),'/aligned/561Visconv_',sprintf('%03d',slice),'.tif');
    filenameout.convVis647 = strcat(ISanalysisfolder,sprintf('%04d',slice),'/aligned/647Visconv_',sprintf('%03d',slice),'.tif');
    filenameout.convIR647 = strcat(ISanalysisfolder,sprintf('%04d',slice),'/aligned/647IRconv_',sprintf('%03d',slice),'.tif');
    filenameout.convIR750 = strcat(ISanalysisfolder,sprintf('%04d',slice),'/aligned/750IRconv_',sprintf('%03d',slice),'.tif');
    % filenames for saving merged STORM, conventional, and alignment images
    filenameout.stormmerge      = strcat(mergestorm,sprintf('%03d', slice),'.tif');
    filenameout.mergeconv       = strcat(mergeconv,sprintf('%03d', slice),'.tif');
    filenameout.conv561only     = strcat(conv561only,sprintf('%03d', slice),'.tif');
    filenameout.mergeconv488561 = strcat(mergeconv488561,sprintf('%03d', slice),'.tif');
    '''

    for slice in range(0,slices):
        
        filename_storm488 = ISanalysisfolder + str(slice).zfill(4)  + "/rawimages/for_matlab/" + "488storm_" + str(slice).zfill(3) + ".tiff";
        filename_storm561 = ISanalysisfolder + str(slice).zfill(4)  + "/rawimages/for_matlab/" + "561storm_" + str(slice).zfill(3) + ".tiff";
        filename_storm647 = ISanalysisfolder + str(slice).zfill(4)  + "/rawimages/for_matlab/" + "647storm_" + str(slice).zfill(3) + ".tiff";
        filename_storm750 = ISanalysisfolder + str(slice).zfill(4)  + "/rawimages/for_matlab/" + "750storm_" + str(slice).zfill(3) + ".tiff";

        filename_convVis488 = ISanalysisfolder + str(slice).zfill(4)  + "/rawimages/for_matlab/" + "488Visconv_" + str(slice).zfill(3) + ".tif";
        filename_convVis561 = ISanalysisfolder + str(slice).zfill(4)  + "/rawimages/for_matlab/" + "561Visconv_" + str(slice).zfill(3) + ".tif";
        filename_convVis647 = ISanalysisfolder + str(slice).zfill(4)  + "/rawimages/for_matlab/" + "647Visconv_" + str(slice).zfill(3) + ".tif";   
        filename_convIR647  = ISanalysisfolder + str(slice).zfill(4)  + "/rawimages/for_matlab/" + "647IRconv_" + str(slice).zfill(3) + ".tif";
        filename_convIR750  = ISanalysisfolder + str(slice).zfill(4)  + "/rawimages/for_matlab/" + "750IRconv_" + str(slice).zfill(3) + ".tif";   

        filenameout_storm488  = ISanalysisfolder + str(slice).zfill(4)  + "/aligned/" + "488storm_" + str(slice).zfill(3) + ".tif";
        filenameout_storm647  = ISanalysisfolder + str(slice).zfill(4)  + "/aligned/" + "647storm_" + str(slice).zfill(3) + ".tif";
        filenameout_storm561  = ISanalysisfolder + str(slice).zfill(4)  + "/aligned/" + "561storm_" + str(slice).zfill(3) + ".tif";
        filenameout_storm750  = ISanalysisfolder + str(slice).zfill(4)  + "/aligned/" + "750storm_" + str(slice).zfill(3) + ".tif";

        filenameout_convVis488  = ISanalysisfolder + str(slice).zfill(4)  + "/aligned/" + "488Visconv_" + str(slice).zfill(3) + ".tif";
        filenameout_convVis561  = ISanalysisfolder + str(slice).zfill(4)  + "/aligned/" + "561Visconv_" + str(slice).zfill(3) + ".tif";
        filenameout_convVis647  = ISanalysisfolder + str(slice).zfill(4)  + "/aligned/" + "647Visconv_" + str(slice).zfill(3) + ".tif";
        filenameout_convIR647   = ISanalysisfolder + str(slice).zfill(4)  + "/aligned/" + "647IRconv_" + str(slice).zfill(3) + ".tif";
        filenameout_convIR750   = ISanalysisfolder + str(slice).zfill(4)  + "/aligned/" + "750IRconv_" + str(slice).zfill(3) + ".tif";

        filenameout_stormmerge        = mergestorm + str(slice).zfill(3) + ".tif";
        filenameout_mergeconv         = mergeconv + str(slice).zfill(3) + ".tif";
        filenameout_convAlignonly       = convAlignonly + str(slice).zfill(3) + ".tif";
        filenameout_cmergeconv488561  = mergeconv488561 + str(slice).zfill(3) + ".tif";
        filenameout_aligned_storm = aligned_storm + str(slice).zfill(3) + ".tif";

        '''
        % load image data
        im.conv488 = im2double(imread(filename.convVis488));
        im.conv488 = imresize(im.conv488,10)./ str2double(rel_conv_ints(4));
        im.conv488adj = imadjust(im.conv488,stretchlim(im.conv488,[0.8 1]),[0 1]);
        im.conv561 = im2double(imread(filename.convVis561));
        im.conv561 = imresize(im.conv561,10)./str2double(rel_conv_ints(3));
        im.conv561adj = imadjust(im.conv561,stretchlim(im.conv561,[0 1]),[0 1]);
        im.conv647 = im2double(imread(filename.convVis647));
        im.conv647 = imresize(im.conv647,10)./str2double(rel_conv_ints(2));
        im.conv647adj = imadjust(im.conv647,stretchlim(im.conv647,[0 1]),[0 1]);
        im.convIR647 = im2double(imread(filename.convIR647));
        im.convIR647 = imresize(im.convIR647,10)./str2double(rel_conv_ints(2));
        im.convIR647adj = imadjust(im.convIR647,stretchlim(im.convIR647,[0 1]),[0 1]);
        im.conv750 = im2double(imread(filename.convIR750));
        im.conv750 = imresize(im.conv750,10)./str2double(rel_conv_ints(1));
        im.conv750adj = imadjust(im.conv750,stretchlim(im.conv750,[0 1]),[0 1]);
        disp(strcat('images loaded for slice # ',sprintf('%04d',slice)));
        '''
        
        im_conv488    = skimg.img_as_float64(timg.imread(filename_convVis488));
        im_conv488    = imresize(im_conv488,10);
        p2, p98 = np.percentile(im_conv488, (80, 100))
        if p2 < 0:
            p2 = 0.0
        if p98 > 1:
            p98 = 1.0
        im_conv488adj = imadjust(im_conv488, p2, p98,0, 1);
     
        im_conv561    = skimg.img_as_float64(timg.imread(filename_convVis561));
        im_conv561    = imresize(im_conv561,10);
        p2, p98 = np.percentile(im_conv561, (0, 100))
        if p2 < 0:
            p2 = 0.0
        if p98 > 1:
            p98 = 1.0
        im_conv561adj = imadjust(im_conv561, p2, p98,0, 1);
        
        im_conv647    = skimg.img_as_float64(timg.imread(filename_convVis647));
        im_conv647    = imresize(im_conv647,10);
        p2, p98 = np.percentile(im_conv647, (0, 100))
        if p2 < 0:
            p2 = 0.0
        if p98 > 1:
            p98 = 1.0
        im_conv647adj = imadjust(im_conv647, p2, p98 ,0, 1);
        
        im_convIR647    = skimg.img_as_float64(timg.imread(filename_convIR647));
        im_convIR647    = imresize(im_convIR647,10);
        p2, p98 = np.percentile(im_convIR647, (0, 100))
        if p2 < 0:
            p2 = 0.0
        if p98 > 1:
            p98 = 1.0
        im_convIR647adj = imadjust(im_convIR647, p2, p98, 0, 1);
        
        im_conv750     = skimg.img_as_float64(timg.imread(filename_convIR750));
        im_conv750     = imresize(im_conv750,10);
        p2, p98 = np.percentile(im_conv750, (0, 100))
        if p2 < 0:
            p2 = 0.0
        if p98 > 1:
            p98 = 1.0
        im_conv750adj  = imadjust(im_conv750, p2, p98, 0, 1); 
        
        print('images loaded for slice # ' + str(slice).zfill(4));    
        
        '''
        if exist(filename.storm488)==2
        im.storm488 = im2double(imread(filename.storm488));
        else
        im.storm488 = im.conv488adj;
        disp('storm488 is not used, loading cov488 instead');
        end
        if exist(filename.storm561)==2
        im.storm561 = im2double(imread(filename.storm561));
        else
        im.storm561 = im.conv561adj;
        disp('storm561 is not used, loading cov561 instead');
        end
        if exist(filename.storm647)==2 %#ok<*EXIST>
        im.storm647 = im2double(imread(filename.storm647));
        else
        im.storm647 = im.conv647adj;
        disp('storm647 is not used, loading cov647 instead');
        end
        if exist(filename.storm750)==2
        im.storm750 = im2double(imread(filename.storm750));
        else
        im.storm750 = im.conv750adj;
        disp('storm750 is not used, loading cov750 instead');
        end
        '''
        
        """ Run through all channels, even if not all of it is used """ 
        
        if os.path.isfile(filename_storm488):
            im_storm488 = skimg.img_as_float64(timg.imread(filename_storm488));
        else:
            im_storm488 = im_conv488adj;
            print('storm488 is not used, loading cov488 instead');
            
        if os.path.isfile(filename_storm561):
            im_storm561 = skimg.img_as_float64(timg.imread(filename_storm561));
        else:
            im_storm561 = im_conv561adj;
            print('storm561 is not used, loading cov561 instead');     
            
        if os.path.isfile(filename_storm647):
            im_storm647 = skimg.img_as_float64(timg.imread(filename_storm647));
        else:
            im_storm647 = im_conv647adj;
            print('storm647 is not used, loading cov647 instead');         
            
        if os.path.isfile(filename_storm750):
            im_storm750 = skimg.img_as_float64(timg.imread(filename_storm750));
        else:
            im_storm750 = im_conv750adj;
            print('storm750 is not used, loading cov750 instead');         
            
        '''
        %correct for storm to conv drift
        [output.storm488] = dftregistration(fft2(im.conv488adj),fft2(im.storm488),1);
        
        [output.storm561] = dftregistration(fft2(im.conv561adj),fft2(im.storm561),1);
        [output.storm647] = dftregistration(fft2(im.conv647adj),fft2(im.storm647),1);
        [output.storm750] = dftregistration(fft2(im.conv750adj),fft2(im.storm750),1);
        '''        

        print("Start dftregistration")
        output_storm488 = dftregistration(np.fft.fft2(im_conv488adj), np.fft.fft2(im_storm488), 1);
        
        output_storm561 = dftregistration(np.fft.fft2(im_conv561adj), np.fft.fft2(im_storm561), 1);
        
        output_storm647 = dftregistration(np.fft.fft2(im_conv647adj), np.fft.fft2(im_storm647), 1);
        
        output_storm750 = dftregistration(np.fft.fft2(im_conv750adj), np.fft.fft2(im_storm750), 1);
        print("Finished dftregistration")

        '''
        xform647 = [ 1  0  0
              0  1  0
             (output.storm647(4)) (output.storm647(3))  1 ];
        tform_translate647 = maketform('affine',xform647); %#ok<*MTFA1>
        imagesize = size(im.conv647);
        xdata = [1 imagesize(2)];
        ydata = [1 imagesize(1)];
        imreg.storm647 = imtransform(im.storm647, tform_translate647, 'XData',xdata,'YData',ydata);%#ok<*DIMTRNS>
        % 
        xform561 = [ 1  0  0
               0  1  0
              (output.storm561(4)) (output.storm561(3))  1 ];
        tform_translate561 = maketform('affine',xform561); %#ok<*MTFA1>
        imreg.storm561 = imtransform(im.storm561, tform_translate561, 'XData',xdata,'YData',ydata); 
        % 
        xform488 = [  1  0  0
               0  1  0
              (output.storm488(4)) (output.storm488(3))  1 ];
        tform_translate488 = maketform('affine',xform488);
        [imreg.storm488] = imtransform(im.storm488, tform_translate488, 'XData',xdata,'YData',ydata);
        xform750 = [  1  0  0
              0  1  0
            (output.storm750(4)) (output.storm750(3))  1 ];
        tform_translate750 = maketform('affine',xform750);
        [imreg.storm750] = imtransform(im.storm750, tform_translate750, 'XData',xdata,'YData',ydata);    
        '''
        
        xform647 = np.array([[1, 0, (-1)*(output_storm647[2]) ], 
                             [0, 1, (-1)*(output_storm647[1])], 
                             [0, 0,  1 ]]);  
        transposeim_storm647 = np.transpose(im_storm647)
        imreg_storm647 = scipy.ndimage.affine_transform(transposeim_storm647,xform647);
        
        xform561 = np.array([[1, 0, (-1)*(output_storm561[2]) ], 
                             [0, 1, (-1)*(output_storm561[1])], 
                             [0, 0,  1 ]]);  
        transposeim_storm561 = np.transpose(im_storm561)
        imreg_storm561 = scipy.ndimage.affine_transform(transposeim_storm561,xform561);
        
        xform488 = np.array([[1, 0, (-1)*(output_storm488[2]) ], 
                             [0, 1, (-1)*(output_storm488[1])], 
                             [0, 0,  1 ]]);  
        transposeim_storm488 = np.transpose(im_storm488)
        imreg_storm488 = scipy.ndimage.affine_transform(transposeim_storm488,xform488);    
        
        xform750 = np.array([[1, 0, (-1)*(output_storm750[2]) ], 
                             [0, 1, (-1)*(output_storm750[1])], 
                             [0, 0,  1 ]]);  
        transposeim_storm750 = np.transpose(im_storm750)
        imreg_storm750 = scipy.ndimage.affine_transform(transposeim_storm750,xform750);    
        
        '''
        [im.conv488_warp] = imtransform(im.conv488, tform_488_2_647,'XData',[1 imagesize(2)],'YData',[1 imagesize(1)]);
        [im.conv561_warp] = imtransform(im.conv561, tform_561_2_647,'XData',[1 imagesize(2)],'YData',[1 imagesize(1)]);
        [im.conv750_warp] = imtransform(im.conv750, tform_750_2_647,'XData',[1 imagesize(2)],'YData',[1 imagesize(1)]);
        '''
        
        rows=len(im_conv488[:,0]);
        cols =len(im_conv488[0,:]);
        tmp_trans = tform_488_2_647.params;
        M1 = np.float32([tmp_trans[0,:],tmp_trans[1,:]]);
        transposeim_conv488 = np.transpose(im_conv488);
        im_conv488_warp= cv2.warpAffine(transposeim_conv488,M1,(cols,rows));

        rows=len(im_conv561[:,0]);
        cols =len(im_conv561[0,:]);
        tmp_trans = tform_561_2_647.params;
        M1 = np.float32([tmp_trans[0,:],tmp_trans[1,:]]);
        transposeim_conv561 = np.transpose(im_conv561);
        im_conv561_warp= cv2.warpAffine(transposeim_conv561,M1,(cols,rows));   
        
        rows=len(im_conv750[:,0]);
        cols =len(im_conv750[0,:]);
        tmp_trans = tform_750_2_647.params;
        M1 = np.float32([tmp_trans[0,:],tmp_trans[1,:]]);
        transposeim_conv750 = np.transpose(im_conv750);
        im_conv750_warp= cv2.warpAffine(transposeim_conv750,M1,(cols,rows));      
        
        '''
        [im.storm488_warp] = imtransform(imreg.storm488, tform_488_2_647,'XData',[1 imagesize(2)],'YData',[1 imagesize(1)]);
        [im.storm561_warp] = imtransform(imreg.storm561, tform_561_2_647,'XData',[1 imagesize(2)],'YData',[1 imagesize(1)]);
        [im.storm750_warp] = imtransform(imreg.storm750, tform_750_2_647,'XData',[1 imagesize(2)],'YData',[1 imagesize(1)]);
        '''
        rows=len(imreg_storm488[:,0]);
        cols =len(imreg_storm488[0,:]);
        tmp_trans = tform_488_2_647.params;
        M1 = np.float32([tmp_trans[0,:],tmp_trans[1,:]]);
        transposeim_storm488 = np.transpose(imreg_storm488);
        im_storm488_warp = cv2.warpAffine(transposeim_storm488,M1,(cols,rows));
        
        rows=len(imreg_storm561[:,0]);
        cols =len(imreg_storm561[0,:]);
        tmp_trans = tform_561_2_647.params;
        M1 = np.float32([tmp_trans[0,:],tmp_trans[1,:]]);
        transposeim_storm561 = np.transpose(imreg_storm561);
        im_storm561_warp = cv2.warpAffine(transposeim_storm561,M1,(cols,rows));
        
        rows=len(imreg_storm750[:,0]);
        cols =len(imreg_storm750[0,:]);
        tmp_trans = tform_750_2_647.params;
        M1 = np.float32([tmp_trans[0,:],tmp_trans[1,:]]);
        transposeim_storm750 = np.transpose(imreg_storm750);
        im_storm750_warp = cv2.warpAffine(transposeim_storm750,M1,(cols,rows));
        
        print('storm images # ' + str(slice).zfill(4) + ' aligned');    
        
        '''
        %save files
        imwrite(im.conv488_warp, filenameout.convVis488);
        imwrite(im.conv561_warp, filenameout.convVis561);
        imwrite(im.conv647,      filenameout.convVis647);
        imwrite(im.conv750_warp, filenameout.convIR750);
        '''
        
        '''
        Last step before write to image files 
        '''

        im_conv488_warp_t = np.transpose(im_conv488_warp);
        #im_conv488_warp_w = im_conv488_warp_t.astype(np.float32);
        im_conv488_warp_w = (im_conv488_warp_t*255.999).astype(np.uint8);
        
        im_conv561_warp_t = np.transpose(im_conv561_warp);
        #im_conv561_warp_w = im_conv561_warp_t.astype(np.float32);  
        im_conv561_warp_w = (im_conv561_warp_t*255.999).astype(np.uint8);
        
        im_conv647_t      = im_conv647;
        #im_conv647_w      = im_conv647_t.astype(np.float32);
        im_conv647_w = (im_conv647_t*255.999).astype(np.uint8);
        
        im_conv750_warp_t = np.transpose(im_conv750_warp);
        #im_conv750_warp_w = im_conv750_warp_t.astype(np.float32);
        im_conv750_warp_w = (im_conv750_warp_t*255.999).astype(np.uint8);

        imsave(filenameout_convVis488,im_conv488_warp_w);
        imsave(filenameout_convVis561,im_conv561_warp_w);
        imsave(filenameout_convVis647,im_conv647_w);
        imsave(filenameout_convIR750, im_conv750_warp_w);
      
        '''
        imwrite(imreg.storm647, filenameout.storm647);
        imwrite(im.storm561_warp, filenameout.storm561);
        imwrite(im.storm488_warp, filenameout.storm488);
        imwrite(im.storm750_warp, filenameout.storm750);
        disp(strcat('chromatic alignment done for slice #',sprintf('%04d',slice)));
        '''


        imreg_storm647_t = np.transpose(imreg_storm647);
        #imreg_storm647_w = imreg_storm647_t.astype(np.float32);
        imreg_storm647_w = (imreg_storm647_t*255.999).astype(np.uint8);
        
        im_storm561_warp_t = im_storm561_warp;
        #im_storm561_warp_w = im_storm561_warp_t.astype(np.float32);
        im_storm561_warp_w = (im_storm561_warp_t*255.999).astype(np.uint8);
        
        im_storm488_warp_t = im_storm488_warp;
        #im_storm488_warp_w = im_storm488_warp_t.astype(np.float32);
        im_storm488_warp_w = (im_storm488_warp_t*255.999).astype(np.uint8);
        
        im_storm750_warp_t = im_storm750_warp;
        #im_storm750_warp_w = im_storm750_warp_t.astype(np.float32);
        im_storm750_warp_w = (im_storm750_warp_t*255.999).astype(np.uint8);

        imsave(filenameout_storm647,imreg_storm647_w);
        imsave(filenameout_storm561,im_storm561_warp_w);
        imsave(filenameout_storm488,im_storm488_warp_w);
        imsave(filenameout_storm750,im_storm750_warp_w);
        
        print('chromatic alignment done for slice #' + str(slice).zfill(4));      
        
        '''
        im.conv_rgb = cat(3, im.conv647,im.conv750_warp,im.conv488_warp,im.conv561_warp);
        imwrite(im.conv_rgb, filenameout.mergeconv);
        %im.storm488_warp = zeros(size(im.storm647));
        im.storm_rgb = cat(3, imreg.storm647,im.storm750_warp,im.storm488_warp,im.storm561_warp);
        imwrite(im.storm_rgb, filenameout.stormmerge);
        imwrite(im.conv561_warp, filenameout.conv561only);    
        '''

    ## this saveout needs to be corrected to generate RGB outputs

        '''
        im_conv_rgb = np.array([im_conv750_warp_w,im_conv647_w, im_conv488_warp_w,im_conv561_warp_w]);
        im_conv_rgb_t = np.transpose(im_conv_rgb);
        im_conv_rgb_w = im_conv_rgb_t.astype(np.uint16);
        gim_conv_rgb = imsave(filenameout_mergeconv, im_conv_rgb_w);
        '''
        
        #im_conv_rgb = (np.dstack((im_conv647_w,im_conv750_warp_w,im_conv488_warp_w,im_conv561_warp_w)) * 255.999) .astype(np.uint8)
        #imsave(filenameout_mergeconv, im_conv_rgb);

        im_conv_rgb = np.stack((im_conv647_w, im_conv750_warp_w, im_conv488_warp_w, im_conv561_warp_w), axis=2)
        Image.fromarray(im_conv_rgb.astype(np.uint8)).save(filenameout_mergeconv)
        
        '''
        im_storm_rgb = np.array([im_storm750_warp_w,imreg_storm647_w,im_storm488_warp_w, im_storm561_warp_w]);
        im_storm_rgb_t = np.transpose(im_storm_rgb);
        im_storm_rgb_w = im_storm_rgb_t.astype(np.uint16);
        gim_storm_rgb = imsave(filenameout_stormmerge, im_storm_rgb_w);
        '''
        
        #im_storm_rgb = (np.dstack((imreg_storm647_w,im_storm750_warp_w,im_storm488_warp_w,im_storm561_warp_w)) * 255.999) .astype(np.uint8)
        #imsave(filenameout_stormmerge, im_storm_rgb);

        im_storm_rgb = np.stack((im_storm647, im_storm750_warp_w, im_storm488_warp_w, im_storm561_warp_w), axis=2)
        Image.fromarray(im_storm_rgb.astype(np.uint8)).save(filenameout_stormmerge)
        
        '''
        %im.align561488 = imfuse(im.conv561_warp,im.conv488_warp);
        %imwrite(im.align561488, filenameout.mergeconv488561);
        '''
        
        print('Finished writing all images for loop #' + str(slice).zfill(4)); 

    print('image_chrom_align complete');  
    
    return "success" 