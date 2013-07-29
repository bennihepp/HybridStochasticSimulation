inputfilepath = '../jobs/outputs/';
inputfilepattern = '^(?<filename>.+)\.(mat|m)$'
inputfilefilter = '^(?<filename>.+)_traj\.(mat|m)$';
outputfilepath = '../jobs/plots/';

listing = dir(inputfilepath);

for i = 1:length(listing)

    file = listing(i);
    if (file.isdir == 0)
        %inputfilefilter = '^(?<filename>.+)_traj\.(mat|m)$';
        [matchStart, matchEnd, tokenIndices, matchStrings, ...
            tokenStrings, tokenName, splitStrings] = regexpi(file.name, inputfilefilter, 'tokens');
        if length(matchStart) > 0
            continue;
        end
        %inputfilepattern = '^(?<filename>.+)\.(mat|m)$';
        [matchStart, matchEnd, tokenIndices, matchStrings, ...
            tokenStrings, tokenName, splitStrings] = regexpi(file.name, inputfilepattern, 'tokens');
        if length(matchStart) > 0

            inputfilename = tokenName.filename;
            ['Plotting ', inputfilename, ' ...']
            outputfilename = inputfilename;

            useTransparency = true;
            drawStdDevOutlines = false;
            linkAxes = false;
            writeOutput = true;
            writeToExtraFolder = true;
            writeFig = true;
            writePdf = false;
            writeEps = false;
            writeTikz = true;
            writeTikzTex = true;

            plotsimulationsFunc(inputfilepath, inputfilename, ...
                outputfilepath, outputfilename, ...
                useTransparency, drawStdDevOutlines, linkAxes, ...
                writeOutput, writeToExtraFolder, ...
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
                outputfilepath, outputfilename, ...
                useTransparency, drawStdDevOutlines, linkAxes, ...
                writeOutput, writeToExtraFolder, ...
                writeFig, writePdf, writeEps, writeTikz, writeTikzTex);

            close;

        end
    end

end
