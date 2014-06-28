inputfilepath = '../jobs/';
inputfilename = 'test2';
outputfilepath = '../jobs/';
outputfilename = inputfilename;

showPlots = true;
numOfBins = 20;
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

plotSimulationHistogramFunc(inputfilepath, inputfilename, ...
    outputfilepath, outputfilename, ...
    showPlots, numOfBins, useTransparency, drawStdDevOutlines, linkAxes, ...
    writeOutput, writeToOwnFolder, ...
    writeFig, writePdf, writeEps, writeTikz, writeTikzTex);

useTransparency = false;
drawStdDevOutlines = true;
writeFig = false;
writePdf = true;
writeEps = true;
writeTikz = false;
writeTikzTex = false;

plotSimulationHistogramFunc(inputfilepath, inputfilename, ...
    outputfilepath, outputfilename, ...
    showPlots, numOfBins, useTransparency, drawStdDevOutlines, linkAxes, ...
    writeOutput, writeToOwnFolder, ...
    writeFig, writePdf, writeEps, writeTikz, writeTikzTex);
