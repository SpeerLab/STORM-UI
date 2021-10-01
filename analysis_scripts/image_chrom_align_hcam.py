import tifffile as timg;
from tifffile import imsave;
import skimage as skimg;
import scipy
import cv2
from PIL import Image

from py_dftregistration import dftregistration;

from imresize import imresize;
from imresize import imadjust;

arg1 = "C:\\Users\\Colenso\\Desktop\\GUI_XY_test\\";
local_exp =  arg1;
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
conv561only = analysisfolder + "unaligned/conv_561/";
mergeconv488561 = analysisfolder + "unaligned/conv_561_488/";
aligned_storm_561 = analysisfolder + "unaligned/storm_561/";

'''
% determine number of sections
slices = (numel(dir(fullfile(ISanalysisfolder, '0*')))-1);
'''

fileExt = r"0*"
slices=len(list(pathlib.Path(ISanalysisfolder).glob(fileExt)));
#slices = 1

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
    filenameout_conv561only       = conv561only + str(slice).zfill(3) + ".tif";
    filenameout_cmergeconv488561  = mergeconv488561 + str(slice).zfill(3) + ".tif";
    filenameout_aligned_storm_561 = aligned_storm_561 + str(slice).zfill(3) + ".tif";

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