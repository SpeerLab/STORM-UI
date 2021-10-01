%normalize image intensity and threshold images

%equalize histogram across sections
%clear all
%base_path = 'Y:/backup/ysigal/';
%base_path = '/n/contefs1/backup/ysigal/';
%arg1 = 'Z:/Colenso/10_02_18_SAC/analysis/';
exp_folder = arg1;

path2 = [exp_folder 'unaligned/storm_merged/'];
path3 = [exp_folder 'unaligned/conv_488/'];
stormfiles = dir([path2 '*.tif']);

num_images = numel(stormfiles);
%

%hy1c = zeros(num_images,255);
%hy2c = zeros(num_images,255);
hy3c = zeros(num_images,255);
hy4c = zeros(num_images,255);
pc=parcluster('local');
delete(pc.Jobs);
pc.NumWorkers=12;
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool(12)
end
parfor j=1:num_images
    disp(j)
    
    A = double(imread([path2 stormfiles(j,1).name]));

    A3 = A(:,:,3);
    A3a = A3(find(A3));
    [hy hx] = hist(A3a,1:1:255);
    hy3c(j,:) = hy;
    
    A = double(imread([path3 stormfiles(j,1).name]));
    A4 = A(:,:,1);
    A4a = A4(find(A4));
    [hy hx] = hist(A4a,1:1:255);
    hy4c(j,:) = hy;
end

x_hist = 1:255;
x_sec = 1:num_images;
%%
for i=1:numel(stormfiles)
    hy3cb(i,:) = hy3c(i,:)./sum(hy3c(i,:));
    hy4cb(i,:) = hy4c(i,:)./sum(hy4c(i,:));
end
%%
% chan = hy3cb;
% for idx = 1:255
%     disp(idx)
%     zthresh = 3;
%     currparam = chan(:,idx);
% currfitx = x_sec(~(zscore(currparam) <-zthresh | zscore(currparam) >zthresh));
% currfity = currparam(~(zscore(currparam) <-zthresh | zscore(currparam) >zthresh));
% [xData, yData] = prepareCurveData( currfitx', currfity );
% 
% % Set up fittype and options.
% ft = fittype( 'smoothingspline' );
% opts = fitoptions( 'Method', 'SmoothingSpline' );
% opts.Normalize = 'on';
% 
% opts.SmoothingParam = 0.8;
% 
% % Fit model to data.
% [fitresult, gof] = fit( xData, yData, ft, opts );
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult,  xData',  yData' );
% varuse3(idx,:) = feval(fitresult,1:num_images);
% end
%
chan = hy4cb;
for idx = 1:255
    disp(idx)
    zthresh = 3;
    currparam = chan(:,idx);
currfitx = x_sec(~(zscore(currparam) <-zthresh | zscore(currparam) >zthresh));
currfity = currparam(~(zscore(currparam) <-zthresh | zscore(currparam) >zthresh));
%[xData, yData, weights] = prepareCurveData( currfitx', currfity, 5*(717-currfitx)' );
[xData, yData] = prepareCurveData( currfitx', currfity );

% Set up fittype and options.
ft = fittype( 'smoothingspline' );
opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.Normalize = 'on';
%opts.SmoothingParam = 0.9999999646249706;
%opts.Weights = weights;

opts.SmoothingParam = 0.9;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
%figure( 'Name', 'untitled fit 1' );
%h = plot( fitresult,  xData',  yData' );
varuse4(idx,:) = feval(fitresult,1:num_images);
end
%
%figure;
%fplot(hy1cb(:,100))
path4 = [exp_folder  'unaligned/for_align/'];
path4a = [exp_folder  'unaligned/for_align_ds/'];
if exist(path4)~=7
    mkdir(path4);
	mkdir([path4a])
end
num_images = numel(stormfiles);
%
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool(32)
end
parfor i=1:num_images
    disp(i)
% hgram3 = varuse3(:,i)/sum(varuse3(:,i));
hgram4 = varuse4(:,i)/sum(varuse4(:,i));
%
%     A = (imread([path2 stormfiles(i,1).name]));
%     A3 = A(:,:,3);
    %
    A = (imread([path3 stormfiles(i,1).name]));
    A4 = A(:,:,1);

% hgram3a = cat(1,numel(find(A3==0)),hgram3*numel(find(A3)));
hgram4a = cat(1,numel(find(A4==0)),hgram4*numel(find(A4)));
%
% B3 = histeq(A3,hgram3a);
B4 = histeq(A4,hgram4a);

% B3(A3<1)=0;
B4(A4<1)=0;

%     Balign = imadd(B4,B3);
    Balign = B4;
    Balign_small = imresize(Balign, 0.1);

    imwrite(Balign, [path4 sprintf('%03d',i) '.png']);
    imwrite(Balign_small, [path4a sprintf('%03d',i) '.png']); 


end
delete(gcp);
quit
