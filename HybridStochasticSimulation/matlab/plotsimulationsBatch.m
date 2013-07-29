inputfilepath = '../jobs/outputs/';
outputfilepath = '../jobs/plots2/';

listing = dir(inputfilepath);

for i = 1:length(listing)

    file = listing(i);
    if (file.isdir == 0)
        pattern = '^(?<filename>.+)_traj\.(mat|m)$';
        [matchStart, matchEnd, tokenIndices, matchStrings, ...
            tokenStrings, tokenName, splitStrings] = regexpi(file.name, pattern, 'tokens');
        if length(matchStart) > 0
            continue;
        end
        pattern = '^(?<filename>.+)\.(mat|m)$';
        [matchStart, matchEnd, tokenIndices, matchStrings, ...
            tokenStrings, tokenName, splitStrings] = regexpi(file.name, pattern, 'tokens');
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

            plotsimulationsFunc(inputfilepath, inputfilename, ...
                outputfilepath, outputfilename, ...
                useTransparency, drawStdDevOutlines, linkAxes, ...
                writeOutput, writeToExtraFolder, ...
                writeFig, writePdf, writeEps, writeTikz);

            close;

            useTransparency = false;
            drawStdDevOutlines = true;
            writeFig = false;
            writePdf = true;
            writeEps = true;
            writeTikz = false;

            plotsimulationsFunc(inputfilepath, inputfilename, ...
                outputfilepath, outputfilename, ...
                useTransparency, drawStdDevOutlines, linkAxes, ...
                writeOutput, writeToExtraFolder, ...
                writeFig, writePdf, writeEps, writeTikz);

            close;

        end
    end

end
