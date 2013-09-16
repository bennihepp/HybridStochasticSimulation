function save_plot(outputfilepath)

    writeFig = true;
    writePdf = true;
    writeEps = true;
    writeTikz = true;
    writeTikzTex = true;

    tf = isdir(outputfilepath);
    if ~tf
        mkdir(outputfilepath);
    end
    outputfileprefix = 'plot';

    if writeFig
        saveas(gcf(), [outputfilepath, '/', outputfileprefix, '.fig'], 'fig');
    end
    if writePdf
        print(gcf(), [outputfilepath, '/',outputfileprefix, '.pdf'], '-dpdf', '-r600', '-cmyk', '-painters');
    end
    if writeEps
        print(gcf(), [outputfilepath, '/',outputfileprefix, '.eps'], '-depsc2', '-r600', '-cmyk', '-painters');
    end
    if writeTikz
        % matfig2pgf('filename', [outputfilepath, outputfileprefix, '.pgf'], 'figwidth', paperSize(1));
        matlab2tikz('filename', [outputfilepath, '/',outputfileprefix, '.tikz'], ...
                    'checkForUpdates', false, 'showInfo', false);
        %             'width', [int2str(round(paperSize(1))), 'cm']
        %             'standalone', true
    end
    if writeTikzTex
        % matfig2pgf('filename', [outputfilepath, outputfileprefix, '.pgf'], 'figwidth', paperSize(1));
        matlab2tikz('filename', [outputfilepath, '/',outputfileprefix, '.tex'], ...
                    'checkForUpdates', false, 'showInfo', false, 'standalone', true);
        %             'width', [int2str(round(paperSize(1))), 'cm']
        %             'standalone', true
    end

end
