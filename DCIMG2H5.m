clear all;

metadata.mainFolder='C:\Users\Z\Documents\SLab\20260210\obj16X08W_ASAP6c_M1\roi2';
metadata.outputFolder = 'results';
metadata.fps = 1000; % hz
metadata.exposure = 0.996;  % ms
metadata.bias = 100;
metadata.fullwellcap = 15000;
metadata.bitdepth = 16;
metadata.quantumeff = 0.8;

nPreview = 10;


% Create folders and collect dcimg files
metadata.h5Path = fullfile(metadata.mainFolder, metadata.outputFolder, "dataset.h5");
metadata.savePath = fullfile(metadata.mainFolder, metadata.outputFolder);
if ~exist(metadata.savePath, 'dir')
    mkdir(metadata.savePath);
end
dcimgFileList = dir(fullfile(metadata.mainFolder, '**\*.dcimg'));
a= {dcimgFileList.folder}';
[~,idx]=unique(a,'stable');
dcimgFileList=dcimgFileList(idx);
if isempty(dcimgFileList)
    error('no extract output file detected in any subfolder')
end

% Read a few frames in tiff for preview and roi determination in ImageJ
options.parallel = false;
[movie,~,summary]=loadDCIMG(fullfile(metadata.mainFolder, dcimgFileList(1).name),[1,nPreview], options);
imwrite(uint16(movie(:,:,1)), fullfile(metadata.savePath, 'InitialFrame.tif'));
for i = 2:size(movie,3)
    imwrite(uint16(movie(:,:,i)), fullfile(metadata.savePath, 'InitialFrame.tif'), 'WriteMode', 'append');
end

clear movie options
%% Load & Convert .dcimg data

options.cropROI = [702, 0, 350, 200];
% options.resize = true;
% options.scale_factor = 0.5;
% [movie,~,summary]=loadDCIMG(fullfile(mainFolder, dcimgFileList(1).name),[1,100], options);
options.h5Path = fullfile(metadata.savePath, "dataset.h5");
options.frameRange = [1,200000];
options.binning = 2;
options.parallel = true;
options.imshow = true;
options.releaseOnNframe = 100000;
[~,~,summary]=loadDCIMGchunks(fullfile(metadata.mainFolder, dcimgFileList(1).name), options);

%%
[bpFilter]=findBestFilterParameters(metadata.h5Path);

save(fullfile(metadata.savePath, 'metadata.mat'), 'metadata', 'options', 'bpFilter');


