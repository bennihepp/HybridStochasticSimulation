function [hlist] = plotSimulationHistogramFunc(inputfilepath, inputfilename, ...
    outputfilepath, outputfileprefix, ...
    showPlots, numOfBins, useTransparency, drawStdDevOutlines, linkAxes, ...
    writeOutput, writeToOwnFolder, ...
    writeFig, writePdf, writeEps, writeTikz, writeTikzTex)

    hlist = [];

    S = load([inputfilepath, inputfilename], 'simulations');

    h = figure;
    if showPlots
        set(h, 'Visible', 'on');
    else
        set(h, 'Visible', 'off');
    end
    hlist(end+1) = h;
    clf;
    orient landscape;
    set(h, 'PaperUnits', 'centimeters');
    set(h, 'PaperType', 'A4');
    paperSize = get(h, 'PaperSize');
    % set(gcf, 'Renderer', 'opengl');

    clear plotScales;
    %plotScales = [1000, 1000, 1];

    xlimMin = 0;
    ylimMin = 0;
    xlimMax = 0;
    ylimMax = 0;
    axArr = [];
    rows = 1;
    cols = 1;

    distributionTimepointIndex = -1;

    l = 1;
    simulationIndex = 1;
    % simulationIndices = [1];
    if simulationIndex <= length(S.simulations);

        simulation = S.simulations{simulationIndex};

        numOfTrajectories = length(simulation);
        if numOfTrajectories > 0
            numOfSpecies = length(simulation(1).trajectories);
        else
            numOfSpecies = 0;
        end

        dataset = zeros(numOfSpecies, numOfTrajectories);
        trajectoryIndices = 1:numOfTrajectories;
        speciesIndices = 1:numOfSpecies;
        titles = [];
        for s=speciesIndices
            name = simulation(1).trajectories(s).name;
            titles = horzcat(titles, name);
        end
    %     plotIndices = [3,4];
        for i=trajectoryIndices
%             traj = trajectories(i);
    %             ax2 = subplot(rows, cols, l);
    %             axArr(end + 1) = ax2;
    %     %         if l == 1
    %     %             ax1 = ax2;
    %     %         end
    %     %         if l == 4
    %     %             linkaxes([ax1, ax2], 'y');
    %     %         end
    %             xlabel('time in seconds');
    %             ylabel('copy numbers');
    %             title(plt.title);
    %             trajIndices = 1:length(plt.trajectories);
    %         %     if i == 6
    %         %         trajIndices = [2, 3];
    %         %     end
    %             colorOrder = distinguishable_colors(length(trajIndices));
    %             labels = cell(1,length(trajIndices));
    %             plotHandles = zeros(1,length(trajIndices));
    %             m = 1;
    %             tSeries = traj.tSeries;
            for s=speciesIndices
                xSeries = simulation(i).trajectories(s).xSeries;
                if distributionTimepointIndex < 0
                    dti = length(xSeries);
                else
                    dti = distributionTimepointIndex;
                end
                if exist('plotScales', 'var')
                    zSeries = plotScales(s) .* xSeries;
                else
                    zSeries = xSeries;
                end
                dataset(s, i) = zSeries(dti);
    %                 hold on;
    %                 colorIndex = mod(m - 1, length(colorOrder)) + 1;
    %                 color = colorOrder(colorIndex,:);
    %                 if isfield(traj, 'xSeriesStdDev')
    %                     fadeColor = ones(size(color)) - 0.1 * (ones(size(color)) - color);
    %                     xSeriesStdDev = traj.xSeriesStdDev;
    %                     if exist('plotScales', 'var')
    %                         zSeriesStdDev = plotScales(j) .* xSeriesStdDev;
    %                     else
    %                         zSeriesStdDev = xSeriesStdDev;
    %                     end
    %                     upperBound = zSeries + zSeriesStdDev;
    %                     lowerBound = zSeries - zSeriesStdDev;
    %                     T = [tSeries; flipud(tSeries)];
    %                     Z = [lowerBound; flipud(upperBound)];
    %                     if useTransparency
    %                         fill(T, Z, color, 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    %                     else
    %                         fill(T, Z, fadeColor, 'EdgeColor', 'none');
    %                     end
    %                     if drawStdDevOutlines
    %                         plot(tSeries, zSeries + zSeriesStdDev, ':', 'Color', color);
    %                         plot(tSeries, zSeries - zSeriesStdDev, ':', 'Color', color);
    %                     end
    %                 end
    %                 plotHandles(m) = plot(tSeries, zSeries, 'Color', color);
    %                 labels{m} = traj.name;
    %                 hold off;
    %                 m = m + 1;
    % %             labels = char(labels);
    % %             legend(plotHandles, labels);
    %             %legend(plotHandles, labels, 'Location', 'Best');
    %             xlimits = xlim(ax2);
    %             ylimits = ylim(ax2);
    %             if xlimits(2) > xlimMax
    %                 xlimMax = xlimits(2);
    %             end
    %             if ylimits(2) > ylimMax
    %                 ylimMax = ylimits(2);
    %             end
    %             if xlimits(1) < xlimMin
    %                 xlimMin = xlimits(1);
    %             end
    %             if ylimits(1) > ylimMin
    %                 ylimMin = ylimits(1);
    %             end
    %             l = l + 1;
            end
        end
    %for i=1:1
    %    linkaxes([axArr(i),axArr(1+i)],'y');
    %end
    % if linkAxes
    %     for i=plotIndices
    %         subplot(rows, cols, i);
    %         xlim([xlimMin, xlimMax]);
    %         ylim([ylimMin, ylimMax]);
    %     end
    % end

        for k=1:size(dataset, 1)
                    ax = subplot(rows, cols, l);
                    xlabel('copy number');
                    ylabel('relative occurence');
                    title(titles(k));
                    hist(dataset(k, :), numOfBins);
                %     if i == 6
                %         trajIndices = [2, 3];
                %     end
%                     colorOrder = distinguishable_colors(length(trajIndices));
%                     labels = cell(1,length(trajIndices));
%                     plotHandles = zeros(1,length(trajIndices));
%                     m = 1;
%                     tSeries = traj.tSeries;
        end

        if writeOutput
            %opts = struct('width', 7, 'height', 7, 'Resolution', 600, 'Color', 'CMYK');
            % 11.7 x 8.3
            %opts = struct('Resolution', 600, 'Color', 'CMYK');
            %exportfig(gcf, [outputfilepath, outputfileprefix, '.eps'], opts);
            %exportfig(gcf, [outputfilepath, outputfileprefix, '.pdf'], opts);

            if writeToOwnFolder
                outputfilepath = [outputfilepath, outputfileprefix, '/'];
                tf = isdir(outputfilepath);
                if ~tf
                    mkdir(outputfilepath);
                end
                outputfileprefix = 'plot';
            end

            if writeFig
                saveas(gcf(), [outputfilepath, outputfileprefix, '.fig'], 'fig');
            end
            if writePdf
                print(gcf(), [outputfilepath, outputfileprefix, '.pdf'], '-dpdf', '-r600', '-cmyk', '-painters');
            end
            if writeEps
                print(gcf(), [outputfilepath, outputfileprefix, '.eps'], '-depsc2', '-r600', '-cmyk', '-painters');
            end
            if writeTikz
                % matfig2pgf('filename', [outputfilepath, outputfileprefix, '.pgf'], 'figwidth', paperSize(1));
                matlab2tikz('filename', [outputfilepath, outputfileprefix, '.tikz'], ...
                            'checkForUpdates', false, 'showInfo', false);
                %             'width', [int2str(round(paperSize(1))), 'cm']
                %             'standalone', true
            end
            if writeTikzTex
                % matfig2pgf('filename', [outputfilepath, outputfileprefix, '.pgf'], 'figwidth', paperSize(1));
                matlab2tikz('filename', [outputfilepath, outputfileprefix, '.tex'], ...
                            'checkForUpdates', false, 'showInfo', false, 'standalone', true);
                %             'width', [int2str(round(paperSize(1))), 'cm']
                %             'standalone', true
            end
        end

    else
        numOfTrajectories = 0;
    %     l = 1;
    end

end

