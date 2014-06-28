function [tSeries, xSeries] = readSimulationDataSingle(filename, dataset, index)
    [tSeries, xSeries] = readSimulationData(filename, dataset, index, 1);
end
