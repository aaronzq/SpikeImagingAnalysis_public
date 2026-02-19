%% Setup file path
clear all

metadata.mainFolder = 'C:\Users\Z\Documents\SLab\20260129\obj16X08W_ASAP6c_M1\roi3';
metadata.mocoFolder = 'demotion';
metadata.saveFolder = 'output';

metadata.fps = 1000; % hz
metadata.exposure = 0.996;  % ms
metadata.bias = 100;
metadata.fullwellcap = 15000;
metadata.bitdepth = 16;
metadata.quantumeff = 0.8;
metadata.roi = [455 0 840 195]; % in the original coordinate system [x0, y0, w, h], x0, y0 starts at 0
metadata.scaling = 1; % image scaling after roi cropping

metadata.savePath = fullfile(metadata.mainFolder, metadata.saveFolder);
if ~exist(metadata.savePath, 'dir')
    mkdir(metadata.savePath);
end

h5file = fullfile(metadata.savePath, 'dataset.h5');

%% Read tiff movie from disk and save in h5

fileList = dir(fullfile(metadata.mainFolder, metadata.mocoFolder));
isFile = ~[fileList.isdir];
fileList = fileList(isFile);

[~,~,ext] = cellfun(@fileparts,{fileList.name},'UniformOutput',false);
tifMask = strcmpi(ext,'.tif');
imgFiles = fileList(tifMask);
[~,order] = sort({imgFiles.name});
imgNames = {imgFiles(order).name};

% read tiff files into workspace and save in mat
% temp = imread(fullfile(mainFolder, mocoFolder, imgNames{1}));
dataset = single(zeros([metadata.roi(4), metadata.roi(3), length(imgNames)]));
for t=1:length(imgNames)
    temp = double(imread(fullfile(metadata.mainFolder, metadata.mocoFolder, imgNames{t})));
    temp = (temp-metadata.bias) .* metadata.fullwellcap ./ 2^(metadata.bitdepth) ./ metadata.quantumeff; % subtract the CMOS bias and convert to photons 
    dataset(:,:,t) = temp(metadata.roi(2)+1:metadata.roi(2)+metadata.roi(4), metadata.roi(1)+1:metadata.roi(1)+metadata.roi(3));
    if t==1
        figure; imagesc(dataset(:,:,t)); colormap(gca, gray); colorbar;
        while true
            a = input("Confirm crop (keyboard y/n): ", 's');
            switch a
                case 'y'
                    break;
                case 'n'
                    return;
                otherwise
                    disp('Type y or n to proceed. Try again.')
            end
        end
        close gcf
    end
    t
end
% Optionally save the dataset to a .mat file
% save(fullfile(metadata.savePath, 'dataset.mat'), 'dataset', 'metadata', '-v7.3');

datasetAVG = uint16(mean(dataset,3));
imwrite(datasetAVG, fullfile(metadata.savePath, "avg.tif"));

disp('Saving data as h5 file...')
if ~exist(fullfile(metadata.savePath, 'dataset.h5'))
    h5create(h5file,'/movie',size(dataset),'Datatype','single');
end
h5write(h5file, '/movie', single(dataset)); 
% write all metadata into h5
metaFields = fieldnames(metadata);
for i = 1:numel(metaFields)
    name = metaFields{i};
    value = metadata.(name);
    attrName = name;
    try
        if ischar(value) || isstring(value)
            h5writeatt(h5file,'/movie',attrName,char(value));
        elseif isnumeric(value) || islogical(value)
            % numeric or logical: store as double array attribute
            h5writeatt(h5file,'/movie',attrName,double(value));
        elseif iscell(value)
            % convert cell array of strings to single string
            try
                h5writeatt(h5file,'/movie',attrName,cellfun(@char,value,'UniformOutput',false));
            catch
                h5writeatt(h5file,'/movie',attrName,mat2str(value));
            end
        else
            % fallback: store as string representation
            h5writeatt(h5file,'/movie',attrName,mat2str(value));
        end
    catch ME
        warning('Could not write metadata field "%s": %s', name, ME.message);
    end
end
disp('Data succesfully saved.')

clear dataset % release memory

%% DETREND the movie

detrending(h5file, 'samplingRate', metadata.fps,'spatialChunk', false);


%% DEMIX THE MOVIE TO FIND SINGLE NEURONS

path=strrep(h5file,'.h5','_dtr.h5');
tic;output=runEXTRACT(fullfile(metadata.savePath, 'dataset_dtr.h5'),'polarityGEVI','pos','cellRadius',20,'removeBackground',true,'method','robust');toc;
datasetAVG = imread(fullfile(metadata.savePath, "avg.tif"));
figure; show_cell_overlay(output, 1:size(output.spatial_weights,3), 'base', datasetAVG, 'contour_thresh', 0.3, 'label', 1);


%% MANUAL CURING OF EXTRACT OUTPUT
% the function will run recursively through all EXTRACT output file found
% in all subfolders from the parent path 'mainFolder'
datasetAVG = imread(fullfile(metadata.savePath, "avg.tif"));
cleanExtractFiles(metadata.savePath, datasetAVG);


%% Only use the spatial filters from EXTRACT and get the raw signals. 
d = dir(fullfile(metadata.savePath, 'DemixingEXTRACT', 'dataset_dtr', '*_clean*.mat'));
if isempty(d)
    error('No _clean .mat file found in %s', fullfile(metadata.savePath, 'DemixingEXTRACT', 'dataset_dtr'));
end
cleanMatFile = fullfile(d(1).folder, d(1).name);
load(cleanMatFile);
datasetAVG = imread(fullfile(metadata.savePath, "avg.tif"));
figure; show_cell_overlay(output, output.cellID, 'base', datasetAVG, 'contour_thresh', 0.3, 'label', 1);

figure; plot(output.temporal_weights(:,output.cellID));

[traces] = getRawTraces(h5file, output.spatial_weights(:,:,output.cellID));
figure; plot(traces);

%% Calculate f0, dff, and detect spikes

id = 3;


t = 1/(metadata.fps)*(1:size(traces,1));
f0 = movmedian(traces(:,id),500);
% dff = traces(:,id)./f0-1;
dff = output.temporal_weights(:,output.cellID(id));
figure; 
plot(traces(:,id)); hold on; plot(f0,'r'); title(['Cell ID: ' num2str(output.cellID(id))])


hi = 0.02;
lo = 0.018;

dff_f = schmitt(dff, hi, lo);
spikes1 = diff([0;dff_f])>0.5;
spikes2 = diff([0;dff_f])<-0.5;
spikes_ind_up = find(spikes1);
spikes_ind_down= find(spikes2);
spikes_ind = zeros(size(spikes_ind_up));

%adjust spike position to the spike peak
for s = 1:length(spikes_ind_up) 
    temp = dff(spikes_ind_up(s):spikes_ind_down(s));
    spikes_ind(s) = spikes_ind_up(s) + find(temp==max(temp)) - 1;
end

%remove spikes with durations less than threshold
spike_duration_threshold = 1;
for s = 1:length(spikes_ind_up)
    if spikes_ind_down(s)-spikes_ind_up(s) < spike_duration_threshold
            spikes_ind(s)=0;
    end
end
spikes_ind = spikes_ind(spikes_ind>0);

dff_stat = dff(spikes_ind);
f0_stat = f0(spikes_ind);



figure;
plot(t, dff, 'Color',[0.0,0.0,0.0]); hold on;
plot(t(spikes_ind),dff(spikes_ind),'r.', 'MarkerSize',10); hold off;
title(['Cell ID: ' num2str(output.cellID(id))])
set(gcf, 'Color', 'w');
% figure;
% histogram(dff_stat);
% title(['Cell ID: ' num2str(output.cellID(id))])

% Box plot for correlation distribution
figure; boxplot(dff_stat, 'Colors', 'k', 'Orientation', 'vertical', 'Symbol',' ')
hold on;

hs = scatter(ones(length(dff_stat),1), dff_stat, ones(length(dff_stat),1)*8,...
    repmat([255, 40, 5]/255,[length(dff_stat),1]),...
    'filled', 'MarkerFaceAlpha', 1, 'jitter','on','JitterAmount',0.1);

h = findobj(gca,'Tag','Box');
hh=h.Parent;
% patch(get(h(1),'XData'),get(h(1),'YData'),[0.26,0.74,0.14],'FaceAlpha',.5);
for ci=1:length(hh.Children)
    set(hh.Children(ci), 'LineWidth', 1.5);
end
set(gcf, 'Color', 'w');
set(gcf, 'Position', [258.3333333333333,456.3333333333333,414.6666666666667,377.3333333333333])
title(['Cell ID: ' num2str(output.cellID(id))])
hold off


waveforms = {};
wid = 1;
wave_halfwidth = 5;
for ss = 1:length(spikes_ind)
    
    if (spikes_ind(ss)-wave_halfwidth > 0) && (spikes_ind(ss)+wave_halfwidth <= size(traces,1))
        
        temp = dff(spikes_ind(ss)-wave_halfwidth:spikes_ind(ss)+wave_halfwidth);
        waveforms{wid} = temp;
        plot(temp); hold on
        wid=wid+1;
    end
end



waveforms = cell2mat(waveforms);
waveforms = waveforms';
figure; plot(mean(waveforms,1))


options.handle     = figure;
options.color_area = [0 0 0]./255;    % Black theme
options.color_line = [ 0 0 0]./255;
% options.color_area = [128 193 219]./255;    % Blue theme
% options.color_line = [ 52 148 186]./255;
% options.color_area = [240 191 74]./255;    % Orange theme
% options.color_line = [240 191 74]./255;
options.alpha      = 0.3;
options.line_width = 2;
options.error      = 'std';
options.x_axis = 1000/1000*(1:size(waveforms,2));
plot_areaerrorbar(waveforms,options)
set(gcf, 'Color', 'w');
set(gcf, "Position", [289,480,506,363]);
ax = gca;
% ax.TickDir = 'none';
axis tight
ax.LineWidth = 1;
xlim([1 1000/1000*size(waveforms,2)])
% ylim([-0.05 0.20])

hold on;
mean_wave = mean(waveforms,1);
plot(options.x_axis(round(0.5*wave_halfwidth)+1:round(1.5*wave_halfwidth)+1), mean_wave(round(0.5*wave_halfwidth)+1:round(1.5*wave_halfwidth)+1), 'k.', 'MarkerSize', 20, 'MarkerEdgeColor',[0 0 0]./255)
