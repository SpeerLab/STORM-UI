%cropping images

%read in ds conv images
delete(gcp)
parpool(28)
%%
%base_path = 'Y:/backup/ysigal/';
exp_folder = arg1
%base_path = 'C:/Batch_analysis_colenso/';
%exp_folder = 'Z:/Colenso/10_02_18_SAC/analysis/';
%base_path = '/n/contefs1/backup/ysigal/';
%exp_folder = [base_path 'REIT_3rd_proc/'];


path = [exp_folder  'rigid_align/'];

stormpath = [path 'storm_merged/'];
convpath =  [path 'conv_merged/'];
convpathds =  [path 'conv_merged_ds/'];
alignpath = [path 'for_align/'];
alignfiles = [dir([alignpath '*.tif']); dir([alignpath '*.png'])];
wgapath = [path 'conv_488/'];
wgafiles = [dir([wgapath '*.tif']); dir([wgapath '*.png'])];
stormfiles = [dir([stormpath '*.tif']); dir([stormpath '*.png'])];
convfiles = [dir([convpath '*.tif']); dir([convpath '*.png'])];
convfilesds = [dir([convpathds '*.tif']); dir([convpathds '*.png'])];
numel(convfilesds)
numel(convfiles)
num_images = numel(stormfiles);
info = imfinfo([convpathds convfilesds(1,1).name]);
%%
disp('normalizing conv-1')
C1 = zeros((info(1,1).Height),(info(1,1).Width),3,num_images,'uint8');

%%normalize intensity of conv_merged images
parfor k = 1:num_images
    A = imread([convpathds convfilesds(k,1).name]);
    C1(:,:,:,k) = A;
end
C2 = max(C1,[],4);

%%
%determine angle of rotation

ang = -8;
C3 = imrotate(C2,ang);

figure;
imshow(C3(:,:,:))
%%
%select region containing data to crop
figure;
[C4 crop_reg] = imcrop(C3(:,:,:));
imshow(C4)
%imrect
%adjust crop region for full size image
crop_storm = ceil(crop_reg*10)
crop_reg = ceil(crop_reg)
save([path 'rot_crop.mat'],'crop_reg','ang')
%%

if exist([exp_folder 'cropped/'])~=7
    mkdir([exp_folder 'cropped/']);
	mkdir([exp_folder 'cropped/storm_merged/'])
	mkdir([exp_folder 'cropped/storm_merged_ds/'])
	mkdir([exp_folder 'cropped/conv_merged/'])
	mkdir([exp_folder 'cropped/conv_merged_ds/'])	
    mkdir([exp_folder 'cropped/conv_488/'])
	mkdir([exp_folder 'cropped/conv_488_ds/'])	
    mkdir([exp_folder 'cropped/for_align/'])
	mkdir([exp_folder 'cropped/for_align_ds/'])
end
parfor k = 1:num_images
    disp(k)
    A1 = imread([stormpath stormfiles(k,1).name]);
    A2 = imrotate(A1,ang);
    A3 = imcrop(A2,crop_storm);
    imwrite(A3, [exp_folder 'cropped/storm_merged/' ...
        sprintf('%03d',k) '.tif']);
    A3_small = imresize(A3, 0.1);
    imwrite(A3_small, [exp_folder 'cropped/storm_merged_ds/' ...
        sprintf('%03d',k) '.tif']);
    
    A1c = imread([convpath convfiles(k,1).name]);
    A2c = imrotate(A1c,ang);
    A3c = imcrop(A2c,crop_storm);
    imwrite(A3c, [exp_folder 'cropped/conv_merged/' ...
        sprintf('%03d',k) '.tif']);
    A3c_small = imresize(A3c, 0.1);
    imwrite(A3c_small, [exp_folder 'cropped/conv_merged_ds/' ...
        sprintf('%03d',k) '.tif']);
    
    A1w1 = imread([alignpath alignfiles(k,1).name]);
    A2w1 = imrotate(A1w1,ang);
    A3w1 = imcrop(A2w1,crop_storm);
    imwrite(A3w1, [exp_folder 'cropped/for_align/' ...
        sprintf('%03d',k) '.tif']);
    A3w1_small = imresize(A3w1, 0.1);
    imwrite(A3w1_small, [exp_folder 'cropped/for_align_ds/' ...
        sprintf('%03d',k) '.tif']);
    
        A1w = imread([wgapath wgafiles(k,1).name]);
    A2w = imrotate(A1w,ang);
    A3w = imcrop(A2w,crop_storm);
    imwrite(A3w, [exp_folder 'cropped/conv_488/' ...
        sprintf('%03d',k) '.tif']);
    A3w_small = imresize(A3w, 0.1);
    imwrite(A3w_small, [exp_folder 'cropped/conv_488_ds/' ...
        sprintf('%03d',k) '.tif']);
end

		%apply angle and crop to fullsize images and save out
        %
%        load in new ds images and check max project

path = [exp_folder  'cropped/'];

stormpath = [path 'storm_merged_ds/'];
convpath =  [path 'conv_merged_ds/'];
convpathds =  [path 'conv_merged_ds/'];
wgapath = [path 'for_align_ds/'];
wgafiles = [dir([wgapath '*.tif']) dir([wgapath '*.png'])];
stormfiles = [dir([stormpath '*.tif']) dir([stormpath '*.png'])];
convfiles = [dir([convpath '*.tif']) dir([convpath '*.png'])];
convfilesds = [dir([convpathds '*.tif']) dir([convpathds '*.png'])];
numel(convfilesds)
numel(convfiles)
num_images = numel(stormfiles);
info = imfinfo([convpathds convfilesds(1,1).name]);
%
disp('normalizing conv-1')
C1 = zeros((info(1,1).Height),(info(1,1).Width),3,num_images,'uint8');

%%normalize intensity of conv_merged images
parfor k = 1:num_images
    A = imread([convpath convfiles(k,1).name]);
    C1(:,:,:,k) = A;
end
C2 = max(C1,[],4);
%figure;
imwrite(C2,[path 'conv_xy_projection.tif'])

C2 = permute(squeeze(max(C1,[],2)),[1 3 2]);
%figure;
imwrite(C2,[path 'conv_xz_projection.tif'])

C2 = permute(squeeze(max(C1,[],1)),[1 3 2]);
%figure;
imwrite(C2,[path 'conv_yz_projection.tif'])

C1 = zeros((info(1,1).Height),(info(1,1).Width),3,num_images,'uint8');

%%normalize intensity of conv_merged images
parfor k = 1:num_images
    A = imread([stormpath stormfiles(k,1).name]);
    C1(:,:,:,k) = A;
end
C2 = max(C1,[],4);
%figure;
imwrite(C2,[path 'storm_xy_projection.tif'])

C2 = permute(squeeze(max(C1,[],2)),[1 3 2]);
%figure;
imwrite(C2,[path 'storm_xz_projection.tif'])

C2 = permute(squeeze(max(C1,[],1)),[1 3 2]);
%figure;
imwrite(C2,[path 'storm_yz_projection.tif'])
%
C1 = zeros((info(1,1).Height),(info(1,1).Width),num_images,'uint8');

%%normalize intensity of conv_merged images
parfor k = 1:num_images
    A = imread([wgapath wgafiles(k,1).name]);
    C1(:,:,k) = A;
end
C2 = max(C1,[],3);
%figure;
imwrite(C2,[path 'wga_xy_projection.tif'])

C2 = squeeze(max(C1,[],2));
%figure;
imwrite(C2,[path 'wga_xz_projection.tif'])

C2 = squeeze(max(C1,[],1));
%figure;
imwrite(C2,[path 'wga_yz_projection.tif'])

delete(gcp);
quit