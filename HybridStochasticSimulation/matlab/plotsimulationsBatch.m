inputfilepath = '../jobs/outputs/';
outputfilepath = '../jobs/plots/';

includepattern = '^(?<filename>.+)\.(?<extension>mat|m)$'
excludepattern = '^(?<filename>.+)_traj\.(mat|m)$';

listing = dir(inputfilepath);

for i = 1:length(listing)

    file = listing(i);
    if (file.isdir ~= 0)
        continue;
    end

    %excludepattern = '^(?<filename>.+)_traj\.(mat|m)$';
    [matchStart, matchEnd, tokenIndices, matchStrings, ...
        tokenStrings, tokenName, splitStrings] = regexpi(file.name, excludepattern, 'tokens');
    if length(matchStart) > 0
        continue;
    end

    %includepattern = '^(?<filename>.+)\.(mat|m)$';
    [matchStart, matchEnd, tokenIndices, matchStrings, ...
        tokenStrings, tokenName, splitStrings] = regexpi(file.name, includepattern, 'tokens');
    if length(matchStart) == 0
        continue;
    end

    inputfilename = [tokenName.filename, '.', tokenName.extension];
    ['Plotting ', inputfilename, ' ...']
    outputfileprefix = tokenName.filename;

    useTransparency = true;
    drawStdDevOutlines = false;
    linkAxes = false;
    writeOutput = true;
    writeToOwnFolder = true;
    writeFig = true;
    writePdf = false;
    writeEps = false;
    writeTikz = true;
    writeTikzTex = true;

    plotsimulationsFunc(inputfilepath, inputfilename, ...
        outputfilepath, outputfileprefix, ...
        useTransparency, drawStdDevOutlines, linkAxes, ...
        writeOutput, writeToOwnFolder, ...
        writeFig, writePdf, writeEps, writeTikz, writeTikzTex);

    close;

    useTransparency = false;
    drawStdDevOutlines = true;
    writeFig = false;
    writePdf = true;
    writeEps = true;
    writeTikz = false;
    writeTikzTex = false;

    plotsimulationsFunc(inputfilepath, inputfilename, ...
        outputfilepath, outputfileprefix, ...
        useTransparency, drawStdDevOutlines, linkAxes, ...
        writeOutput, writeToOwnFolder, ...
        writeFig, writePdf, writeEps, writeTikz, writeTikzTex);

    close;

end
