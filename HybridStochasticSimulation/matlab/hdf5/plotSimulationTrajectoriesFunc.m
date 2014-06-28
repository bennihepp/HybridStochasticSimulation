function [hlist] = plotSimulationTrajectoriesFunc( ...
    tSeries, xSeries, speciesList, trajectoryList, ...
    outputfilepath, outputfileprefix, ...
    showPlots, useTransparency, drawStdDevOutlines, linkAxes, ...
    writeOutput, writeToOwnFolder, ...
    writeFig, writePdf, writeEps, writeTikz, writeTikzTex)

    hlist = [];

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

    if length(trajectoryList) == 0
        trajectoryList = 1:size(data,1);
    end
    % simulationIndices = [1];
    %rows = S.rows;
    %cols = S.cols;
    cols = 1;
    numOfPlots = length(trajectoryList);
    rows = ceil(numOfPlots / cols);
    % rows = 2;
    % cols = 2;
    l = 1;
    for i=trajectoryList
        ax2 = subplot(rows, cols, l);
        axArr(end + 1) = ax2;
%         if l == 1
%             ax1 = ax2;
%         end
%         if l == 4
%             linkaxes([ax1, ax2], 'y');
%         end
        xlabel('time in seconds');
        ylabel('copy numbers');
        %title(plt.title);
        %trajIndices = 1:length(plt.trajectories);
    %     if i == 6
    %         trajIndices = [2, 3];
    %     end
        colorOrder = distinguishable_colors(length(speciesList));
        %labels = cell(1,length(trajIndices));
        %plotHandles = zeros(1,length(trajIndices));
        m = 1;
        for j=speciesList
            hold on;
            colorIndex = mod(m - 1, length(colorOrder)) + 1;
            color = colorOrder(colorIndex,:);
            plot(tSeries, xSeries(i,j,:), 'Color', color);
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
    % if linkAxes
    %     for i=plotIndices
    %         subplot(rows, cols, i);
    %         xlim([xlimMin, xlimMax]);
    %         ylim([ylimMin, ylimMax]);
    %     end
    % end

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

end
