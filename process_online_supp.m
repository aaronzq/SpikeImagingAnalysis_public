clear all;
% installSIA();

roiName = '20260219\obj16X08W_ASAP6c_M1\roi3';
h5Path = fullfile('C:\Users\Z\Documents\SLab', roiName, 'results/dataset.h5');

load(strrep(h5Path,'dataset.h5','metadata.mat'), 'bpFilter', 'metadata', 'options');
fps = metadata.fps;
binning = options.binning;
range = [1,200000];
savePath = fullfile('C:\Users\Z\Documents\SLab', roiName, 'results');
%% Manual clean cell candidates
if isempty(range)
    path=char(strrep(h5Path,'.h5', '_bp_moco.h5'));
else
    path=char(strrep(h5Path,'.h5', ['_bp_moco' num2str(range(1)) '_' num2str(range(2)) '.h5']));
end
movieDims = h5info(path).Datasets.Dataspace.Size;
dataset_sample = h5read(path, '/mov', [1 1 1], [movieDims(1) movieDims(2) 1000]);
imwrite(uint16(mean(dataset_sample,3)), fullfile(savePath, "avg.tif"));

datasetAVG = imread(fullfile(savePath, "avg.tif")); % avg.tif has to be created manually from dataset_bp_moco.h5
figure; show_cell_overlay(output, 1:size(output.spatial_weights,3), 'base', datasetAVG, 'contour_thresh', 0.3, 'label', 1);
cleanExtractFiles(savePath, datasetAVG);

