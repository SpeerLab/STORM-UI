import sys
import glob
import os
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.image as mpimg

import scipy
import imageio

import cv2
import pathlib
import skimage

import PIL
from PIL import Image

from interactive_crop import interactive_crop 

def crop_datastack(exp_folder, alignment_channel): 
    
    # Assert alignment_channel is correct 
    assert alignment_channel in [561, 488]
    
    path = exp_folder + "\\rigid_align\\"
    storm_path = path + "\\storm_merged\\"
    conv_path = path + "\\conv_merged\\"
    conv_path_ds = path + "\\conv_merged_ds\\"
    align_path = path + "\\for_align\\"
    wga_path = path + "\\conv_{}\\".format(str(alignment_channel))
    
    tif_ext = r"*.tif"
    png_ext = r"*.png"
    
    wga_files = glob.glob(wga_path + tif_ext) + glob.glob(wga_path + png_ext) 
    storm_files = glob.glob(storm_path + tif_ext) + glob.glob(storm_path + png_ext) 
    conv_files = glob.glob(conv_path + tif_ext) + glob.glob(conv_path + png_ext) 
    conv_files_ds = glob.glob(conv_path_ds + tif_ext) + glob.glob(conv_path_ds + png_ext) 
    
    num_images = len(storm_files) 
    
    C1 = []
    
    print("normalizing conv-1")

    # normalize intensity of conv_merged images
    for k in range(num_images):
        A = imageio.imread(conv_files_ds[k])
        C1.append(A)
    C1 = np.stack(C1, -1) # h x w x 3 x n_images 
    C2 = C1.max(axis=3)
    
    # determine angle of rotation 
    ang = -8 
    C3 = skimage.transform.rotate(C2,ang)

    # interactive crop 
    crop_reg = interactive_crop(C3)
    x_start, y_start, x_end, y_end = crop_reg 
    C4 = C3[x_start:x_end, y_start:y_end, :] 
    print('Cropping Coordinates: {}, Angle: {}'.format(crop_reg, ang))

    crop_storm = np.ceil(crop_reg*10)
    crop_storm_x_range = list(range(crop_storm[0], crop_storm[2]))
    crop_storm_y_range = list(range(crop_storm[1], crop_storm[3]))
    
    crop_reg = np.ceil(crop_reg)
    
    print('Making directories!') 
    if not os.path.exists(exp_folder + "\\cropped\\"):
        os.mkdir(exp_folder + "\\cropped\\")
        os.mkdir(exp_folder + "\\cropped\\storm_merged\\")
        os.mkdir(exp_folder + "\\cropped\\storm_merged_ds\\")
        os.mkdir(exp_folder + "\\cropped\\conv_merged\\")
        os.mkdir(exp_folder + "\\cropped\\conv_merged_ds\\") 
        os.mkdir(exp_folder + "\\cropped\\conv_{}\\".format(alignment_channel))
        os.mkdir(exp_folder + "\\cropped\\conv_{}_ds\\".format(alignment_channel))
        os.mkdir(exp_folder + "\\cropped\\for_align\\")
        os.mkdir(exp_folder + "\\cropped\\for_align_ds\\")
        
    print('Applying angle and crop to full size images and then saving out!') 
   
    for i in range(num_images): 
        print(k)
        A1 = cv2.imread(storm_path+storm_files[k])
        A2 = skimage.transform.rotate(A1,ang)
        A3 = (A2[crop_storm_x_range, crop_storm_y_range]).astype(np.float32)
        cv2.imwrite(exp_folder+"/cropped/storm_merged/"+str(k).zfill(3) +'.tif', A3)
        A3_small = skimage.transform.imresize(A3, 0.1)
        cv2.imwrite(exp_folder+"/cropped/storm_merged_ds/"+str(k).zfill(3) +'.tif', A3_small)

        A1c = cv2.imread(conv_path+conv_files[k])
        A2c = skimage.transform.rotate(A1c,ang)
        A3c = (A2[crop_storm_x_range, crop_storm_y_range]).astype(np.float32)
        cv2.imwrite(exp_folder+"/cropped/conv_merged/"+str(k).zfill(3) +'.tif', A3c)
        A3c_small = skimage.transform.imresize(A3c, 0.1)
        cv2.imwrite(exp_folder+"/cropped/conv_merged_ds/"+str(k).zfill(3) +'.tif', A3_small)

        A1w1 = cv2.imread(align_path+align_files[k])
        A2w1 = skimage.transform.rotate(A1w1, ang)
        A3w1 = (A2w1[crop_storm_x_range, crop_storm_y_range]).astype(np.float32)        
        cv2.imwrite(exp_folder+"/cropped/for_align/"+str(k).zfill(3) +'.tif', A3w1)
        A3w1_small = skimage.transform.imresize(A3w1, 0.1)
        cv2.imwrite(exp_folder+"/cropped/for_align_ds/"+str(k).zfill(3) +'.tif', A3w1_small)
        
        A1w = cv2.imread(wga_path+wga_files[k])
        A2w = skimage.transform.rotate(A1w,ang)
        A3w = (A2w[crop_storm_x_range, crop_storm_y_range]).astype(np.float32)
        cv2.imwrite(exp_folder+"/cropped/conv_{}/".format(alignment_channel)+str(k).zfill(3) +'.tif', A3w)
        A3w_small = skimage.transform.imresize(A3w, 0.1)
        cv2.imwrite(exp_folder+"/cropped/conv_{}_ds/".format(alignment_channel)+str(k).zfill(3) +'.tif', A3w_small)

    # Load in new ds images and check max project 
    
    path = exp_folder + "\\cropped\\"
    storm_path = exp_folder + "\\cropped\\storm_merged_ds\\"
    conv_path = exp_folder + "\\cropped\\conv_merged_ds\\"
    conv_path_ds = exp_folder + "\\cropped\\conv_{}_ds\\".format(alignment_channel)
    wga_path = exp_folder + "\\cropped\\for_align_ds\\"
       
    wga_files = glob.glob(wga_path + tif_ext) + glob.glob(wga_path + png_ext) 
    storm_files = glob.glob(storm_path + tif_ext) + glob.glob(storm_path + png_ext) 
    conv_files = glob.glob(conv_path + tif_ext) + glob.glob(conv_path + png_ext) 
    conv_files_ds = glob.glob(conv_path_ds + tif_ext) + glob.glob(conv_path_ds + png_ext) 

    num_images = len(storm_files) 
    
    C1 = []
    
    print("normalizing conv-1 2.0")

    # normalize intensity of conv_merged images
    for k in range(num_images):
        A = imageio.imread(conv_files_ds[k])
        C1.append(A)
    C1 = np.stack(C1, -1) # h x w x 3 x n_images 
    
    cv2.imwrite(exp_folder + '/cropped/conv_xy_projection.tif', 
                C1.max(axis=3)) # h x w x 3 
    cv2.imwrite(exp_folder + '/cropped/conv_xz_projection.tif', 
                C1.max(axis=1).squeeze().transpose(0,2,1))  # h x 3 x n_images -> h x n_images x 3 
    cv2.imwrite(exp_folder + '/cropped/conv_yz_projection.tif', 
                C1.max(axis=0).squeeze().transpose(0,2,1))  # w x 3 x n_images -> w x n_images x 3
                
                
    print("normalizing storm-1 2.0")

    # normalize intensity of storm_merged images
    for k in range(num_images):
        A = imageio.imread(storm_files[k])
        C1.append(A)
    C1 = np.stack(C1, -1) # h x w x 3 x n_images 
    
    cv2.imwrite(exp_folder + '/cropped/storm_xy_projection.tif', 
                C1.max(axis=3)) # h x w x 3 
    cv2.imwrite(exp_folder + '/cropped/storm_xz_projection.tif', 
                C1.max(axis=1).squeeze().transpose(0,2,1))  # h x 3 x n_images -> h x n_images x 3 
    cv2.imwrite(exp_folder + '/cropped/storm_yz_projection.tif', 
                C1.max(axis=0).squeeze().transpose(0,2,1))  # w x 3 x n_images -> w x n_images x 3
                
    print("normalizing wga-1 2.0")
    # normalize intensity of aligned images
    for k in range(num_images):
        A = imageio.imread(wga_files[k])
        C1.append(A)
    C1 = np.stack(C1, -1) # h x w x 3 x n_images 
    
    cv2.imwrite(exp_folder + '/cropped/wga_xy_projection.tif', 
                C1.max(axis=3)) # h x w x 3 
    cv2.imwrite(exp_folder + '/cropped/wga_xz_projection.tif', 
                C1.max(axis=1).squeeze().transpose(0,2,1))  # h x 3 x n_images -> h x n_images x 3 
    cv2.imwrite(exp_folder + '/cropped/wga_yz_projection.tif', 
                C1.max(axis=0).squeeze().transpose(0,2,1))  # w x 3 x n_images -> w x n_images x 3
                
    print('Done with cropping!')         
    
    
    return True
    
    
    

