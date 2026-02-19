clear all;

mainFolder='C:\Users\Z\Documents\SLab\20260210\obj16X08W_ASAP6c_M1\roi2';
outputFolder = 'results';
h5Path = fullfile(mainFolder, outputFolder, "dataset.h5");
savePath = fullfile(mainFolder, outputFolder);
if ~exist(savePath, 'dir')
    mkdir(savePath);
end


dcimgFileList = dir(fullfile(mainFolder, '**\*.dcimg'));
a= {dcimgFileList.folder}';
[~,idx]=unique(a,'stable');
dcimgFileList=dcimgFileList(idx);

if isempty(dcimgFileList)
    error('no extract output file detected in any subfolder')
end
options.parallel = false;
[movie,~,summary]=loadDCIMG(fullfile(mainFolder, dcimgFileList(1).name),[1,10], options);
imwrite(uint16(movie(:,:,1)), fullfile(mainFolder, outputFolder, 'InitialFrame.tif'));
for i = 2:size(movie,3)
    imwrite(uint16(movie(:,:,i)), fullfile(mainFolder, outputFolder, 'InitialFrame.tif'), 'WriteMode', 'append');
end
%% Load & Convert .dcimg data


options.cropROI = [702, 0, 350, 200];

% options.resize = true;
% options.scale_factor = 0.5;
% [movie,~,summary]=loadDCIMG(fullfile(mainFolder, dcimgFileList(1).name),[1,100], options);

options.h5Path = fullfile(mainFolder, outputFolder, "dataset.h5");
options.frameRange = [1,200000];
options.binning = 2;
options.parallel = true;
options.imshow = true;
[movie,~,summary]=loadDCIMGchunks(fullfile(mainFolder, dcimgFileList(1).name), options);


%% Bandpass and motion correction

[bpFilter]=findBestFilterParameters(h5Path);
bandPassMovieChunk(h5Path, bpFilter);
path=char(strrep(h5Path,'.h5', '_bp.h5'));
motionCorr1Movie(path,'nonRigid', false,'isRawInput',false,'dcRemoval',false);

%% Detrend
fps = 1000;

path=char(strrep(h5Path,'.h5', '_bp_moco.h5'));
detrending(path, 'samplingRate', fps,'spatialChunk', true);

%% Extract (demix)
path=char(strrep(h5Path,'.h5','_bp_moco_dtr.h5'));
tic;output=runEXTRACT(path,'polarityGEVI','pos','cellRadius',15,'removeBackground',true,'method','robust');toc;
datasetAVG = imread(fullfile(savePath, "avg.tif"));
figure; show_cell_overlay(output, 1:size(output.spatial_weights,3), 'base', datasetAVG, 'contour_thresh', 0.3, 'label', 1);

%% Manual clean cell candidates
datasetAVG = imread(fullfile(savePath, "avg.tif"));
cleanExtractFiles(savePath, datasetAVG);

%% Only use the spatial filters from EXTRACT and get the raw signals. 
d = dir(fullfile(savePath, 'DemixingEXTRACT', 'dataset_bp_moco_dtr', '*_clean*.mat'));
if isempty(d)
    error('No _clean .mat file found in %s', fullfile(metadata.savePath, 'DemixingEXTRACT', 'dataset_dtr'));
end
cleanMatFile = fullfile(d(1).folder, d(1).name);
load(cleanMatFile);
datasetAVG = imread(fullfile(savePath, "avg.tif"));
figure; show_cell_overlay(output, output.cellID, 'base', datasetAVG, 'contour_thresh', 0.3, 'label', 1);
figure; plot(output.temporal_weights(:,output.cellID));
path=char(strrep(h5Path,'.h5', '_bp_moco.h5'));
[traces] = getRawTraces(path, output.spatial_weights(:,:,output.cellID));
figure; plot(traces);

%%
%% Calculate f0, dff, and detect spikes

id = 2;


t = 1/(fps)*(1:size(traces,1));
f0 = movmedian(traces(:,id),1000);
% dff = traces(:,id)./f0-1;
dff = output.temporal_weights(:,output.cellID(id));
figure; 
plot(traces(:,id)); hold on; plot(f0,'r'); title(['Cell ID: ' num2str(output.cellID(id))])


hi = 0.035;
lo = 0.03;

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
dp_stat = dff_stat .* sqrt(f0_stat.*fps.*(3e-3)/2);


figure;
plot(t, dff, 'Color',[0.0,0.0,0.0]); hold on;
plot(t(spikes_ind),dff(spikes_ind),'r.', 'MarkerSize',10); hold off;
title(['Cell ID: ' num2str(output.cellID(id))])
set(gcf, 'Color', 'w');
% figure;
% histogram(dff_stat);
% title(['Cell ID: ' num2str(output.cellID(id))])

%%

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
%%

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



