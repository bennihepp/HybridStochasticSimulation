inputfilepath = 'data/';
inputfilename = 'heatShockResponse_distribution_set2';
inputfilename = '../../simulation';
outputfilepath = 'plots/';
outputfilename = inputfilename;
useTransparency = true;
drawStdDevOutlines = false;
linkAxes = false;
writeOutput = false;

S = load([inputfilepath, inputfilename,'.mat'], 'plots', 'rows', 'cols');

figure;
clf;
orient landscape;
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperType', 'A4');
paperSize = get(gcf, 'PaperSize');
% set(gcf, 'Renderer', 'opengl');

clear plotScales;
% plotScales = [2000, 1];

xlimMin = 0;
ylimMin = 0;
xlimMax = 0;
ylimMax = 0;
axArr = [];
rows = S.rows;
cols = S.cols;
plotIndices = 1:length(S.plots);
%rows = 2;
%cols = 1;
%plotIndices = [1,10];
l = 1;
for i=plotIndices
    plt = S.plots(i);
    ax2 = subplot(rows, cols, l);
    axArr(end + 1) = ax2;
    if l == 1
        ax1 = ax2;
    end
    if l == 4
        linkaxes([ax1, ax2], 'y');
    end
    xlabel('time in seconds');
    ylabel('copy numbers');
    title(plt.title);
    trajIndices = 1:length(plt.trajectories);
%     if i == 6
%         trajIndices = [2, 3];
%     end
    colorOrder = distinguishable_colors(length(trajIndices));
    labels = cell(1,length(trajIndices));
    plotHandles = zeros(1,length(trajIndices));
    m = 1;
    for j=trajIndices
        traj = plt.trajectories(j);
        tSeries = plt.tSeries;
        xSeries = traj.xSeries;
        if exist('plotScales')
            zSeries = plotScales(j) .* xSeries;
        else
            zSeries = traj.xSeries;
        end
        hold on;
        colorIndex = mod(m - 1, length(colorOrder)) + 1;
        color = colorOrder(colorIndex,:);
        fadeColor = ones(size(color)) - 0.1 * (ones(size(color)) - color);
        if isfield(traj, 'xSeriesStdDev')
            xSeriesStdDev = traj.xSeriesStdDev;
            if exist('plotScales')
                zSeriesStdDev = plotScales(j) .* xSeriesStdDev;
            else
                zSeriesStdDev = xSeriesStdDev;
            end
            upperBound = zSeries + zSeriesStdDev;
            lowerBound = zSeries - zSeriesStdDev;
            T = [tSeries; flipud(tSeries)];
            Z = [lowerBound; flipud(upperBound)];
            if useTransparency
                fill(T, Z, color, 'FaceAlpha', 0.1, 'EdgeColor', 'none');
            else
                fill(T, Z, fadeColor, 'EdgeColor', 'none');
            end
            if drawStdDevOutlines
                plot(tSeries, zSeries + zSeriesStdDev, ':', 'Color', color);
                plot(tSeries, zSeries - zSeriesStdDev, ':', 'Color', color);
            end
        end
        plotHandles(m) = plot(tSeries, zSeries, 'Color', color);
        labels{m} = traj.name;
        hold off;
        m = m + 1;
    end
    labels = char(labels);
    legend(plotHandles, labels);
    %legend(plotHandles, labels, 'Location', 'Best');
    xlimits = xlim(ax2);
    ylimits = ylim(ax2);
    if xlimits(2) > xlimMax
        xlimMax = xlimits(2);
    end
    if ylimits(2) > ylimMax
        ylimMax = ylimits(2);
    end
    if xlimits(1) < xlimMin
        xlimMin = xlimits(1);
    end
    if ylimits(1) > ylimMin
        ylimMin = ylimits(1);
    end
    l = l + 1;
end
%for i=1:1
%    linkaxes([axArr(i),axArr(1+i)],'y');
%end
if linkAxes
    for i=plotIndices
        subplot(rows, cols, i);
        xlim([xlimMin, xlimMax]);
        ylim([ylimMin, ylimMax]);
    end
end

if writeOutput
    %opts = struct('width', 7, 'height', 7, 'Resolution', 600, 'Color', 'CMYK');
    % 11.7 x 8.3
    %opts = struct('Resolution', 600, 'Color', 'CMYK');
    %exportfig(gcf, [outputfilepath, outputfilename, '.eps'], opts);
    %exportfig(gcf, [outputfilepath, outputfilename, '.pdf'], opts);

    saveas(gcf(), [outputfilepath, outputfilename, '.fig'], 'fig');
    print(gcf(), [outputfilepath, outputfilename, '.pdf'], '-dpdf', '-r600', '-cmyk', '-painters');
    print(gcf(), [outputfilepath, outputfilename, '.eps'], '-depsc2', '-r600', '-cmyk', '-painters');
    % matfig2pgf('filename', [outputfilepath, outputfilename, '.pgf'], 'figwidth', paperSize(1));
    matlab2tikz('filename', [outputfilepath, outputfilename, '.tikz'], ...
                'checkForUpdates', false, 'showInfo', false);
    %             'width', [int2str(round(paperSize(1))), 'cm']
    %             'standalone', true
end
