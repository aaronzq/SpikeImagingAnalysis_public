clear all;

metadata.mainFolder='C:\Users\Z\Documents\SLab\20260129\obj16X08W_ASAP6c_M1\roi8';
metadata.outputFolder = 'results';
metadata.fps = 1000; % hz
metadata.exposure = 0.996;  % ms
metadata.bias = 100;
metadata.fullwellcap = 15000;
metadata.bitdepth = 16;
metadata.quantumeff = 0.8;

% Create folders
metadata.h5Path = fullfile(metadata.mainFolder, metadata.outputFolder, "dataset.h5");
metadata.savePath = fullfile(metadata.mainFolder, metadata.outputFolder);
if ~exist(metadata.savePath, 'dir')
    mkdir(metadata.savePath);
end

fileList = dir(metadata.mainFolder);
isFile = ~[fileList.isdir];
fileList = fileList(isFile);

[~,~,ext] = cellfun(@fileparts,{fileList.name},'UniformOutput',false);
tifMask = strcmpi(ext,'.tif');
imgFiles = fileList(tifMask);
[~,order] = sort({imgFiles.name});
imgNames = {imgFiles(order).name};

%% Load & Convert tiff data

options.cropROI = [585, 0, 300, 200];
options.h5Path = metadata.h5Path;
options.frameRange = [1,10000];
options.binning = 2;
options.scale_factor = 1/options.binning;
% read tiff files into workspace and save in mat
% temp = imread(fullfile(mainFolder, mocoFolder, imgNames{1}));
temp = zeros([options.cropROI(4), options.cropROI(3)]);
temp = imresize(temp, options.scale_factor, "bilinear");
dataset = single(zeros([size(temp,1), size(temp,2), diff(options.frameRange)+1]));
for t=options.frameRange(1):options.frameRange(2)
    temp = double(imread(fullfile(metadata.mainFolder, imgNames{t})));
    temp = (temp-metadata.bias) .* metadata.fullwellcap ./ 2^(metadata.bitdepth) ./ metadata.quantumeff .* (options.binning)^2;
    temp = temp(options.cropROI(2)+1:options.cropROI(2)+options.cropROI(4), options.cropROI(1)+1:options.cropROI(1)+options.cropROI(3));  
    dataset(:,:,t) = imresize(temp, options.scale_factor, "bilinear");
    t
end

h5append(options.h5Path, single(dataset)); 

clear dataset % release memory
%%
[bpFilter]=findBestFilterParameters(metadata.h5Path);

save(fullfile(metadata.savePath, 'metadata.mat'), 'metadata', 'options', 'bpFilter');


