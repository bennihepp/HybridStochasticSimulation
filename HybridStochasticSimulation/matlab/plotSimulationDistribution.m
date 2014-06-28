inputfilepath = '../jobs/';
inputfilename = 'test';
outputfilepath = '../jobs/';
outputfilename = inputfilename;

showPlots = true;
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

plotSimulationDistributionFunc(inputfilepath, inputfilename, ...
    outputfilepath, outputfilename, ...
    showPlots, useTransparency, drawStdDevOutlines, linkAxes, ...
    writeOutput, writeToOwnFolder, ...
    writeFig, writePdf, writeEps, writeTikz, writeTikzTex);

useTransparency = false;
drawStdDevOutlines = true;
writeFig = false;
writePdf = true;
writeEps = true;
writeTikz = false;
writeTikzTex = false;

plotSimulationDistributionFunc(inputfilepath, inputfilename, ...
    outputfilepath, outputfilename, ...
    showPlots, useTransparency, drawStdDevOutlines, linkAxes, ...
    writeOutput, writeToOwnFolder, ...
    writeFig, writePdf, writeEps, writeTikz, writeTikzTex);
