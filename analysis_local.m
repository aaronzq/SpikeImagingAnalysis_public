clear all;
installSIA();

roiName = '20260129\obj16X08W_ASAP6c_M1\roi10';
h5Path = fullfile('C:\Users\Z\Documents\SLab', roiName, 'results/dataset.h5');

load(strrep(h5Path,'dataset.h5','metadata.mat'), 'bpFilter', 'metadata', 'options');
fps = metadata.fps;
binning = options.binning;
% savePath = metadata.savePath;
savePath = fullfile('C:\Users\Z\Documents\SLab', roiName, 'results');
time_range = [];
% time_range = [1,200000];
% time_range = [200001,400000];
% time_range = [400001,600000];
% time_range = [600001,800000];
% time_range = [800001,1000000];
% time_range = [1000001,1200000];
%% Only use the spatial filters from EXTRACT and get the raw signals. 
d = dir(fullfile(savePath, 'DemixingEXTRACT', 'dataset_bp_moco_dtr', '*_clean*.mat'));
if isempty(d)
    error('No _clean .mat file found in %s', fullfile(savePath, 'DemixingEXTRACT', 'dataset_dtr'));
end
cleanMatFile = fullfile(d(1).folder, d(1).name);
load(cleanMatFile);
datasetAVG = imread(fullfile(savePath, "avg.tif"));
figure; show_cell_overlay(output, output.cellID, 'base', datasetAVG, 'contour_thresh', 0.3, 'label', 1);
figure; plot(output.temporal_weights(:,output.cellID));
if isempty(time_range)
    path=char(strrep(h5Path,'.h5', '_bp_moco.h5'));
else
    path=char(strrep(h5Path,'.h5', ['_bp_moco' num2str(time_range(1)) '_' num2str(time_range(2)) '.h5']));
end
[traces_photon] = getRawTraces(path, output.spatial_weights(:,:,output.cellID));
figure; plot(traces_photon);
%% Calculate f0, dff, and detect spikes


id = 3;


t = 1/(fps)*(1:size(traces_photon,1));
f0 = movmedian(traces_photon(:,id),1000);
dff = traces_photon(:,id)./f0-1;


% dff = output.temporal_weights(:,output.cellID(id));
figure; 
plot(t,traces_photon(:,id),'k'); hold on; plot(t,f0,'r','LineWidth',2); title(['Cell ID: ' num2str(output.cellID(id))]);
ylabel('Photons'); xlabel('Time(s)')


hi = 0.06;
lo = 0.035;

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
% title(['Cell ID: ' num2str(id)])
set(gcf, 'Color', 'w');
% figure;
% histogram(dff_stat);
% title(['Cell ID: ' num2str(output.cellID(id))])


% Box plot for correlation distribution
figure; boxplot(dp_stat, 'Colors', 'k', 'Orientation', 'vertical', 'Symbol',' ')
hold on;

hs = scatter(ones(length(dp_stat),1), dp_stat, ones(length(dp_stat),1)*8,...
    repmat([255, 40, 5]/255,[length(dp_stat),1]),...
    'filled', 'MarkerFaceAlpha', 1, 'jitter','on','JitterAmount',0.1);

h = findobj(gca,'Tag','Box');
hh=h.Parent;
% patch(get(h(1),'XData'),get(h(1),'YData'),[0.26,0.74,0.14],'FaceAlpha',.5);
for ci=1:length(hh.Children)
    set(hh.Children(ci), 'LineWidth', 1.5);
end
set(gcf, 'Color', 'w');
% set(gcf, 'Position', [258.3333333333333,456.3333333333333,414.6666666666667,377.3333333333333])
title(['Cell ID: ' num2str(output.cellID(id))])
hold off
%%
% 
% waveforms = {};
% wid = 1;
% wave_halfwidth = 5;
% for ss = 1:length(spikes_ind)
% 
%     if (spikes_ind(ss)-wave_halfwidth > 0) && (spikes_ind(ss)+wave_halfwidth <= size(traces,1))
% 
%         temp = dff(spikes_ind(ss)-wave_halfwidth:spikes_ind(ss)+wave_halfwidth);
%         waveforms{wid} = temp;
%         plot(temp); hold on
%         wid=wid+1;
%     end
% end
% 
% 
% 
% waveforms = cell2mat(waveforms);
% waveforms = waveforms';
% figure; plot(mean(waveforms,1))
% 
% 
% options.handle     = figure;
% options.color_area = [0 0 0]./255;    % Black theme
% options.color_line = [ 0 0 0]./255;
% % options.color_area = [128 193 219]./255;    % Blue theme
% % options.color_line = [ 52 148 186]./255;
% % options.color_area = [240 191 74]./255;    % Orange theme
% % options.color_line = [240 191 74]./255;
% options.alpha      = 0.3;
% options.line_width = 2;
% options.error      = 'std';
% options.x_axis = 1000/1000*(1:size(waveforms,2));
% plot_areaerrorbar(waveforms,options)
% set(gcf, 'Color', 'w');
% set(gcf, "Position", [289,480,506,363]);
% ax = gca;
% % ax.TickDir = 'none';
% axis tight
% ax.LineWidth = 1;
% xlim([1 1000/1000*size(waveforms,2)])
% % ylim([-0.05 0.20])
% 
% hold on;
% mean_wave = mean(waveforms,1);
% plot(options.x_axis(round(0.5*wave_halfwidth)+1:round(1.5*wave_halfwidth)+1), mean_wave(round(0.5*wave_halfwidth)+1:round(1.5*wave_halfwidth)+1), 'k.', 'MarkerSize', 20, 'MarkerEdgeColor',[0 0 0]./255)
% 
% 
% 
