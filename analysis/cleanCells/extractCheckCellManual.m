function [cellList,figH]=extractCheckCellManual(output,avg, varargin)

% can input your own cell list.
% % options.cellList=[];
options.waitForUser=true;
options.frameRate=1000;
options.positionFig=[0,0,1920,1080]; 

%% UPDATE OPTIONS
if nargin>1
    options=getOptions(options,varargin);
end

%%
if isempty(output.temporal_weights)
    cellList=[];
    figH=figure();
    return;
end

spatial=full(output.spatial_weights);
% temporal=wdenoise(double(output.temporal_weights),1);
temporal = double(output.temporal_weights);
% fs=output.fs;
fs=options.frameRate;
time=getTime(temporal,fs);
decade=1;
p=0;
cellList=[];
h=figure('defaultaxesfontsize',12,'color','w','Position',options.positionFig);
while p<size(temporal,2)
    % if p == 80
    %     aaa=1;
    % end
    nUnits=min(10,size(temporal,2)-p);
    clf(h);
    for i=1:nUnits
        subplot(10,10,[10*(i-1)+1 10*(i-1)+2])
        % imagesc(spatial(:,:,p+i))
        show_cell_overlay(output, p+i, 'base', avg);
        title(num2str(p+i))
        %         axis off
        subplot(10,10,[10*(i-1)+3 10*(i-1)+10])
        plot(time,100*temporal(:,p+i))
        xlim([0 time(end)])
        ylabel('dF/F (%)')
        % title('Press any key to enter your favorite cell number')
    end
    xlabel('Time (s)')
    figure(h);
    if options.waitForUser
        prompt = {'which cells to keep? [as a list, just space e.g. 1 5 8 9 15'};
        dlgtitle = 'Input';
        dims = [2 50];
        % pause 
        figure(h);
        answer = inputdlg(prompt,dlgtitle,dims);
        cellList=[cellList str2num(answer{1})];
    end
    decade=decade+1;
    p=10*(decade-1);
    figure(h);
end

close all

% double check with the user
figH=figure('defaultaxesfontsize',12,'color','w','Name','summary of selected cells','Position',options.positionFig);
nUnits=max(10,numel(cellList));
for i=1:numel(cellList)
    subplot(nUnits,10,[10*(i-1)+1 10*(i-1)+2])
    show_cell_overlay(output, cellList(i), 'base', avg);
    title(num2str(cellList(i)))
    ax = subplot(nUnits,10,[10*(i-1)+3 10*i]);
    yyaxis(ax,'right')
    ax.YColor = 'k';
    plot(time,100*temporal(:,cellList(i)), 'k')
    ylabel('dF/F (%)')
    yyaxis(ax,'left')
    ax.YColor = 'none';
    ax.YTick = [];
    xlim([0 time(end)])
    axis tight
    % title('Press any key to enter your favorite cell number')
end
xlabel('Time (s)')

answer = questdlg('Are you happy with you cell choice? (press any key to continue)', ...
    'Cell selection Summary', ...
    'Yes Happy','No do it again','Yes Happy');
% Handle response
switch answer
    case 'Yes Happy'
        disp('outputting results..')
%         close all
    case 'No do it again'
        close all
        disp('Lets run cell-check again')
        [cellList,figH]=extractCheckCellManual(output, avg);
end
%dont close figure for PDF saving
end