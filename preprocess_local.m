clear all;
installSIA();

roiName = '20260219\obj16X08W_ASAP6c_M1\roi1';
h5Path = fullfile('C:\Users\Z\Documents\SLab', roiName, 'results/dataset.h5');

load(strrep(h5Path,'dataset.h5','metadata.mat'), 'bpFilter', 'metadata', 'options');
fps = metadata.fps;
binning = options.binning;
%%
if isfile('log.txt')
    delete('log.txt');
end
diary log.txt

try
    %% Bandpass and motion correction
    bandPassMovieChunk(h5Path, bpFilter);
    path=char(strrep(metadata.h5Path,'.h5', '_bp.h5'));
    motionCorr1Movie(path,'nonRigid', false,'isRawInput',false,'dcRemoval',false);

    %% Detrend
    path=char(strrep(h5Path,'.h5', '_bp_moco.h5'));
    detrending(path, 'samplingRate', fps,'spatialChunk', true);

    %% Extract (demix)
    path=char(strrep(h5Path,'.h5','_bp_moco_dtr.h5'));
    tic;output=runEXTRACT(path,'polarityGEVI','pos','cellRadius',10,'removeBackground',true,'method','robust');toc;  

    disp('Preprocess finished.');
    out = 1;

catch exception
    fprintf('Error occurred: %s\n', exception.message);
    out = 0; % Indicate failure
end

%% Manual clean cell candidates
path=char(strrep(h5Path,'.h5', '_bp_moco.h5'));
movieDims = h5info(path).Datasets.Dataspace.Size;
dataset_sample = h5read(path, '/mov', [1 1 1], [movieDims(1) movieDims(2) 1000]);
imwrite(uint16(mean(dataset_sample,3)), fullfile(metadata.savePath, "avg.tif"));

datasetAVG = imread(fullfile(metadata.savePath, "avg.tif")); % avg.tif has to be created manually from dataset_bp_moco.h5
figure; show_cell_overlay(output, 1:size(output.spatial_weights,3), 'base', datasetAVG, 'contour_thresh', 0.3, 'label', 1);
cleanExtractFiles(metadata.savePath, datasetAVG);

