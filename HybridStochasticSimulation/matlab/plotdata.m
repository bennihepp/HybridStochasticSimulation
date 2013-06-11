inputfilepath = 'data/';
inputfilename = 'stochasticFocusingNetwork3';
outputfilepath = 'plots/';
outputfilename = inputfilename;

S = load([inputfilepath, inputfilename,'.mat'], 'plots', 'rows', 'cols');

figure;
clf;
orient landscape;
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperType', 'A4');

xlimMin = 0;
ylimMin = 0;
xlimMax = 0;
ylimMax = 0;
for i=1:length(S.plots)
    plt = S.plots(i);
    ax2 = subplot(S.rows, S.cols, i);
    if i == 1
        ax1 = ax2;
    end
    xlabel('time in seconds');
    ylabel('copy numbers');
    title(plt.title);
    colorOrder = distinguishable_colors(length(plt.trajectories));
    labels = cell(1,length(plt.trajectories));
    plotHandles = zeros(1,length(plt.trajectories));
    for j=1:length(plt.trajectories)
        traj = plt.trajectories(j);
        tSeries = plt.tSeries;
        xSeries = traj.xSeries;
        hold on;
        colorIndex = mod(j - 1, length(colorOrder)) + 1;
        plotHandles(j) = plot(tSeries, xSeries, 'Color', colorOrder(colorIndex,:));
        labels{j} = traj.name;
        if isfield(traj, 'xSeriesStdDev')
            xSeriesStdDev = traj.xSeriesStdDev;
            plot(tSeries, xSeries + xSeriesStdDev, ':', 'Color', colorOrder(colorIndex,:));
            plot(tSeries, xSeries - xSeriesStdDev, ':', 'Color', colorOrder(colorIndex,:));
        end
        hold off;
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
end
for i=1:length(S.plots)
    subplot(S.rows, S.cols, i);
    xlim([xlimMin, xlimMax]);
    ylim([ylimMin, ylimMax]);
end

%opts = struct('width', 7, 'height', 7, 'Resolution', 600, 'Color', 'CMYK');
% 11.7 x 8.3
%opts = struct('Resolution', 600, 'Color', 'CMYK');
%exportfig(gcf, [outputfilepath, outputfilename, '.eps'], opts);
%exportfig(gcf, [outputfilepath, outputfilename, '.pdf'], opts);

print(gcf(), [outputfilepath, outputfilename, '.pdf'], '-dpdf', '-r600', '-cmyk', '-painters');
print(gcf(), [outputfilepath, outputfilename, '.eps'], '-depsc2', '-r600', '-cmyk', '-painters');
