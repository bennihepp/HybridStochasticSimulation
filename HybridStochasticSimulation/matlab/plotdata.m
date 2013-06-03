inputfilepath = 'data/';
inputfilename = 'birthDeathTunnelModel1';
outputfilepath = 'plots/';
outputfilename = inputfilename;

S = load([inputfilepath, inputfilename,'.mat'], 'plots', 'rows', 'cols');

figure;
clf;
orient landscape;
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperType', 'A4');

for i=1:length(S.plots)
    plt = S.plots(i);
    subplot(S.rows, S.cols, i);
    xlabel('time in seconds');
    ylabel('copy numbers');
    title(plt.title);
    colorOrder = distinguishable_colors(length(plt.trajectories));
    labels = cell(length(plt.trajectories), 1);
    for j=1:length(plt.trajectories)
        traj = plt.trajectories(j);
        tSeries = plt.tSeries;
        xSeries = traj.xSeries;
        hold on;
        colorIndex = mod(j - 1, length(colorOrder)) + 1;
        plot(tSeries, xSeries, 'Color', colorOrder(colorIndex,:));
        labels{j} = traj.name;
        if isfield(traj, 'xSeriesStdDev')
            xSeriesStdDev = traj.xSeriesStdDev;
            plot(tSeries, xSeries + xSeriesStdDev, 'Color', colorOrder(colorIndex,:));
            plot(tSeries, xSeries - xSeriesStdDev, 'Color', colorOrder(colorIndex,:));
        end
        hold off;
    end
    labels = char(labels);
    legend(labels);
end

%opts = struct('width', 7, 'height', 7, 'Resolution', 600, 'Color', 'CMYK');
% 11.7 x 8.3
%opts = struct('Resolution', 600, 'Color', 'CMYK');
%exportfig(gcf, [outputfilepath, outputfilename, '.eps'], opts);
%exportfig(gcf, [outputfilepath, outputfilename, '.pdf'], opts);

print(gcf(), [outputfilepath, outputfilename, '.pdf'], '-dpdf', '-r600');
print(gcf(), [outputfilepath, outputfilename, '.eps'], '-dpsc2', '-r600');
