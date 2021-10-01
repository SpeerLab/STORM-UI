import numpy as np
import os
import pathlib
import sys
import h5py
from PIL import Image

from analysis_scripts.corr_mols_v8 import corr_mols
# from corr_mols_v8 import corr_mols

import skimage
import skimage.transform
from skimage.transform import (estimate_transform)

import matplotlib.pyplot as plt

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

def gen_bead_warp(local_exp):

    local_exp = os.path.normpath(local_exp)

    acq_path = local_exp + "\\acquisition\\"
    base_pth = acq_path + "\\bins\\"
    save_pth = acq_path + "\\bead_registration\\"

    #% set scale for image expansion
    scale = 10
     
    #% size of field in pixels
    xdiv = 896*scale
    ydiv = 896*scale

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

            import pdb; pdb.set_trace() 

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
    q1 = axs[0, 0].quiver(comb_set561_2_pos_x, comb_set561_2_pos_y, fac*set561_2_orig_error_x, fac*set561_2_orig_error_y, color = 'r', width=0.002,angles='xy', scale_units='xy', scale=1)

    axs[0, 1].set_title('488 to 647, no warping 20x mag. residuals', fontsize = 5, fontweight='bold')
    axs[0, 1].scatter(comb_set488_2_pos_x,comb_set488_2_pos_y, color='lime', marker = ',', s = .005)
    q2 = axs[0, 1].quiver(comb_set488_2_pos_x, comb_set488_2_pos_y, fac*set488_2_orig_error_x, fac*set488_2_orig_error_y, color = 'r', width=0.002,angles='xy', scale_units='xy', scale=1)
      
    axs[1, 0].set_title('561 to 647, affine warp, 20x mag. residuals', fontsize = 5, fontweight='bold')
    axs[1, 0].scatter(comb_set561_2_pos_x,comb_set561_2_pos_y, color='lime', marker = ',', s = .005)
    q3 = axs[1, 0].quiver(comb_set561_2_pos_x, comb_set561_2_pos_y, fac*set561_2_warp_error_x, fac*set561_2_warp_error_y, color = 'r', width=0.002, angles='xy', scale_units='xy', scale=1)

    axs[1, 1].set_title('488 to 647, affine warp, 20x mag. residuals', fontsize = 5, fontweight='bold')
    axs[1, 1].scatter(comb_set488_2_pos_x,comb_set488_2_pos_y, color='lime', marker = ',', s = .005)
    q4 = axs[1, 1].quiver(comb_set488_2_pos_x, comb_set488_2_pos_y, fac*set488_2_warp_error_x, fac*set488_2_warp_error_y, color = 'r', width=0.002,angles='xy', scale_units='xy', scale=1)

    fac=40

    axs[0, 2].set_title('750 to 647, no warping 40x mag. residuals', fontsize = 5, fontweight='bold')
    axs[0, 2].scatter(comb_set750_2_pos_x,comb_set750_2_pos_y, color='lime', marker = ',', s = .005)
    q5 = axs[0, 2].quiver(comb_set750_2_pos_x, comb_set750_2_pos_y, fac*set750_2_orig_error_x, fac*set750_2_orig_error_y, color = 'r', width=0.002,angles='xy', scale_units='xy', scale=1)

    axs[1, 2].set_title('750 to 647, affine warp, 40x mag. residuals', fontsize = 5, fontweight='bold')
    axs[1, 2].scatter(comb_set750_2_pos_x,comb_set750_2_pos_y, color='lime', marker = ',', s = .005)
    q6 = axs[1, 2].quiver(comb_set750_2_pos_x, comb_set750_2_pos_y, fac*set750_2_warp_error_x, fac*set750_2_warp_error_y, color = 'r', width=0.002,angles='xy', scale_units='xy', scale=1)


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
    
    return True 
    
if __name__ == "__main__":
    exp_folder = "C:/Users/Vatsal/QT_Projects/New_Pipeline/4color_test"
    gen_bead_warp(exp_folder)
