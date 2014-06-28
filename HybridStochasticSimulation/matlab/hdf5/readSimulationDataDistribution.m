function [t, xSeries] = readSimulationDataDistribution(filename, dataset, t)
    tSeries = readSimulationTimepoints(filename, dataset);
    k = find(tSeries >= t, 1);
    t = tSeries(k);

    xInfo = readSimulationInfo(filename, [dataset, '/xSeries']);
    xSize = xInfo.Dataspace.Size;
    xStart = ones(1, length(xSize));
    xStart(3) = k;
    xCount = xSize;
    xCount = xCount(length(xCount):-1:1);
    xCount(3) = 1;
    [tSeries, xSeries] = readSimulationData(filename, dataset, xStart, xCount);
    xSeries = squeeze(xSeries);
end
