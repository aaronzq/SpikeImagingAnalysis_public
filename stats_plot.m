% clear all;
installSIA();

% read table from Excel
nCell = 9;
sheetName = {'LED50', 'LED20', 'LED5'};
powerName = {'26 mW/mm2', '16 mW/mm2', '8 mW/mm2'};
powerLevels = [26.02, 16.28, 7.96];
statName = {"dF/F0", "F0", "d'"};
Stats = {};
for idPower = 1:length(sheetName)
    StatsPower = {}
    for idCell = 1:nCell
        T = readtable(fullfile("C:\Users\Z\Documents\SLab\20260129\obj16X08W_ASAP6c_M1", sprintf("cell%d.xlsx", idCell)), 'Sheet', sheetName{idPower});
        StatsPower{idCell} = T;
    end
    Stats{idPower} = StatsPower;
end

% Stats{1}, Stats{2} and Stats{3} are tables for differnt power
% Stats{1}{cellID} are tables for different cells at certain power
% Stats{1}{cellID}{1} for df/f0, Stats{1}{cellID}{2} for f0, Stats{1}{cellID}{3} for d' 
 
%% per cell df/f0 and d' study

cellID = 9;

% dF/F0, from 8, 16 to 26 mW/mm2
d{1} = Stats{3}{cellID}{:,1};
d{2} = Stats{2}{cellID}{:,1};
d{3} = Stats{1}{cellID}{:,1};
data = [d{1}; d{2}; d{3}];
dataCategory = categorical(repelem([powerName(3);powerName(2);powerName(1)],[size(Stats{3}{cellID}{:,1},1),size(Stats{2}{cellID}{:,1},1),size(Stats{1}{cellID}{:,1},1)]));
dataCategory = reordercats(dataCategory, [powerName(3); powerName(2); powerName(1)]);

figure; 
v = violinplot(data, dataCategory);
for iv = 1:3
    v(iv).EdgeColor = [0.3, 0.3, 0.3]; % Set edge color for the violin plots
    v(iv).BoxColor = [0.3, 0.3, 0.3];
    v(iv).BoxWidth = 0.03;
end
set(gcf, 'Color', 'w');
% set(gcf, 'Position', [258.3333333333333,456.3333333333333,414.6666666666667,377.3333333333333])
title(['Cell ID: ' num2str(cellID)]);
ylabel('dF/F0', 'FontWeight', 'bold');
xlabel('Power', 'FontWeight', 'bold');
ylim([0.01, 0.1]);
hold off

% d', from 8, 16 to 26 mW/mm2
d{1} = Stats{3}{cellID}{:,3};
d{2} = Stats{2}{cellID}{:,3};
d{3} = Stats{1}{cellID}{:,3};
data = [d{1}; d{2}; d{3}];
dataCategory = categorical(repelem([powerName(3);powerName(2);powerName(1)],[size(Stats{3}{cellID}{:,1},1),size(Stats{2}{cellID}{:,1},1),size(Stats{1}{cellID}{:,1},1)]));
dataCategory = reordercats(dataCategory, [powerName(3); powerName(2); powerName(1)]);

figure; 
v = violinplot(data, dataCategory);
for iv = 1:3
    v(iv).EdgeColor = [0.3, 0.3, 0.3]; % Set edge color for the violin plots
    v(iv).BoxColor = [0.3, 0.3, 0.3];
    v(iv).BoxWidth = 0.03;
end
yline(6, '--k', 'LineWidth', 1.5);

set(gcf, 'Color', 'w');
% set(gcf, 'Position', [258.3333333333333,456.3333333333333,414.6666666666667,377.3333333333333])
title(['Cell ID: ' num2str(cellID)]);
ylabel("d'", 'FontWeight', 'bold');
xlabel('Power', 'FontWeight', 'bold');
ylim([0, 30]);
hold off

disp([num2str((mean(d{1})-median(d{1}))/std(d{1})), ' ' , num2str((mean(d{2})-median(d{2}))/std(d{2})), ' ' num2str((mean(d{3})-median(d{3}))/std(d{3}))]);

%% inter-cell study, df/f0

dff0Stats = zeros(nCell, length(powerName));
for powerID=1:length(powerName)
    for cellID=1:nCell
        
        dff0Stats(cellID, powerID) = mean(Stats{powerID}{cellID}{:,1});
    
    end
end
statsData = mean(dff0Stats,1);
statsCategory = categorical(repelem([powerName(1);powerName(2);powerName(3)],[1,1,1]));
statsCategory = reordercats(statsCategory, [powerName(3); powerName(2); powerName(1)]);

figure; hb = bar(statsCategory, statsData); hold on;
C = orderedcolors("gem");
hb.FaceColor = 'flat';
hb.CData = C(1:size(hb.CData,1),:);
hb.FaceAlpha = 0.7;

statsData = dff0Stats(:);
statsCategory = categorical(repelem([powerName(1);powerName(2);powerName(3)],[nCell,nCell,nCell]));
statsCategory = reordercats(statsCategory, [powerName(3); powerName(2); powerName(1)]);
statsCategoryOrder = double(statsCategory);
hs = scatter(statsCategoryOrder, statsData, 30*ones(size(statsData)),...
    repmat([0, 0, 0]/255, size(statsData)),...
     'MarkerFaceAlpha', 1, 'jitter','on','JitterAmount',0.1);

sigstar({[1,2],[2,3],[1,3]},[0.17,0.66,0.19]);
hold off;

ylabel('dF/F0', 'FontWeight', 'bold');
xlabel('Power', 'FontWeight', 'bold');
%% inter-cell study, f0

f0Stats = zeros(nCell, length(powerName));
for powerID=1:length(powerName)
    for cellID=1:nCell
        temp = Stats{powerID}{cellID}{:,2};
        f0Stats(cellID, powerID) = 0.5*(max(temp)+min(temp));
    end
end
figure;
hold on;

cmap = parula(nCell);
legendHandles = gobjects(nCell+1,1); % +1 for overall fit

% ----- Plot each cell individually -----
for i = 1:nCell
    
    x = powerLevels;
    y = f0Stats(i,:);
    
    % Raw line
    hLine = plot(x, y, '-', 'Color', cmap(i,:), 'LineWidth', 1.5);
    
    % Solid scatter points
    scatter(x, y, 60, cmap(i,:), 'filled');
    
    % Individual linear fit (extrapolated to 0)
    p = polyfit(x, y, 1);
    xFit = linspace(0, max(x), 100);
    yFit = polyval(p, xFit);
    
    plot(xFit, yFit, ':', 'Color', cmap(i,:), 'LineWidth', 1.5);
    
    legendHandles(i) = hLine;
end

% ----- Overall fit using ALL data points -----
xAll = repelem(powerLevels(:), nCell);
yAll = f0Stats(:);

pOverall = polyfit(xAll, yAll, 1);

xFitAll = linspace(0, max(powerLevels), 200);
yFitAll = polyval(pOverall, xFitAll);

hOverall = plot(xFitAll, yFitAll, 'k:', 'LineWidth', 3);
legendHandles(nCell+1) = hOverall;

% ----- Display equation (top-left corner) -----
slope = pOverall(1);
intercept = pOverall(2);

eqnStr = sprintf('$F_0 = %.3f x %+ .3f$', slope, intercept);

text(0.05, 0.95, eqnStr, ...
    'Units','normalized', ...
    'Interpreter','latex', ...
    'FontSize',12, ...
    'VerticalAlignment','top');

% ----- Labels & formatting -----
xlabel('Power density (mW/mm2)', 'FontWeight', 'bold');
ylabel('F0 (Photons)', 'FontWeight', 'bold');
title('F0 versus power density');

legendStrings = arrayfun(@(i) sprintf('Cell %d', i), 1:nCell, 'UniformOutput', false);
legendStrings{end+1} = 'Overall fit';

legend(legendHandles, legendStrings, 'Location', 'best');

box on;
hold off;

%%
% Box plot for correlation distribution
figure; boxplot(data, g, 'Colors', 'k', 'Orientation', 'vertical', 'Symbol',' ');
hold on;

for pid=1:3
    hs = scatter(pid*ones(size(d{pid},1), 1), d{pid}, ones(size(d{pid},1),1)*8,...
        repmat([255, 40, 5]/255,[size(d{pid},1),1]),...
        'filled', 'MarkerFaceAlpha', 1, 'jitter','on','JitterAmount',0.1);
end



%%
h = findobj(gca,'Tag','Box');
hh=h.Parent;
% patch(get(h(1),'XData'),get(h(1),'YData'),[0.26,0.74,0.14],'FaceAlpha',.5);
for ci=1:length(hh.Children)
    set(hh.Children(ci), 'LineWidth', 1.5);
end
set(gcf, 'Color', 'w');
% set(gcf, 'Position', [258.3333333333333,456.3333333333333,414.6666666666667,377.3333333333333])
title(['Cell ID: ' num2str(cellID)])
hold off





%%

% for each variable/column, read values until first NaN and store in a cell
nVars = width(T);
tabledata = cell(1, nVars);
for k = 1:nVars
    col = T{:, k};                    % extract column as array (numeric or other)
    if isnumeric(col) || islogical(col)
        nanIdx = find(isnan(col), 1);
    else
        % treat empty strings or missing as NaN-equivalent for non-numeric
        nanIdx = find(ismissing(col) | (string(col) == ""), 1);
    end
    if isempty(nanIdx)
        tabledata{k} = col(:);     % all values
    else
        tabledata{k} = col(1:nanIdx-1);
    end
end

g1 = repmat({'Power 26 mW/mm2'},size(tabledata{cellID},1),1);
g2 = repmat({'Power 16 mW/mm2'},size(tabledata{cellID+3},1),1);
g3 = repmat({'Power 8 mW/mm2'},size(tabledata{cellID+6},1),1);
g = [g1; g2; g3];

% Box plot for correlation distribution
figure; boxplot(data, g, 'Colors', 'k', 'Orientation', 'vertical', 'Symbol',' ')
hold on;

for pid=1:3
    hs = scatter(pid*ones(size(tabledata{cellID+pid*3-3},1), 1), tabledata{cellID+pid*3-3}, ones(size(tabledata{cellID+pid*3-3},1),1)*8,...
        repmat([255, 40, 5]/255,[size(tabledata{cellID+pid*3-3},1),1]),...
        'filled', 'MarkerFaceAlpha', 1, 'jitter','on','JitterAmount',0.1);
end


h = findobj(gca,'Tag','Box');
hh=h.Parent;
% patch(get(h(1),'XData'),get(h(1),'YData'),[0.26,0.74,0.14],'FaceAlpha',.5);
for ci=1:length(hh.Children)
    set(hh.Children(ci), 'LineWidth', 1.5);
end
set(gcf, 'Color', 'w');
% set(gcf, 'Position', [258.3333333333333,456.3333333333333,414.6666666666667,377.3333333333333])
title(['Cell ID: ' num2str(cellID)])
hold off
