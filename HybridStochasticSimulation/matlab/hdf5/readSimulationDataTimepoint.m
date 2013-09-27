function [t, xSeries] = readSimulationDataTimepoint(filename, dataset, t)
    tSeries = h5read(filename, ['/simulations/', dataset, '/tSeries']);
    k = find(tSeries >= t, 1);
    t = tSeries(k);

    xInfo = readSimulationInfo(filename, [dataset, '/xSeries']);
    xSize = xInfo.Dataspace.Size;
    xStart = ones(1, length(xSize));
    xStart(3) = k;
    xCount = xSize;
    xCount(3) = 1;
    [tSeries, xSeries] = readSimulationData(filename, dataset, xStart, xCount);
    xSeries = squeeze(xSeries);
end
