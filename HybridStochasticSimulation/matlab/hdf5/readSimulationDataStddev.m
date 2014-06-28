function [tSeries, xStddev] = readSimulationDataStddev(xMean, filename, dataset, batchSize)
    if nargin < 4
        batchSize = 5000;
    end
    xInfo = readSimulationInfo(filename, [dataset, '/xSeries']);
    xSize = xInfo.Dataspace.Size;
    totalSize = xSize(length(xSize));
    batchSize = min(batchSize, xSize(length(xSize)));
    %start = 1;
    xStart = ones(1, length(xSize));
    xCount = xSize(length(xSize):-1:1);
    xCount(1) = batchSize;

%     xCount = xInfo.ChunkSize;
%     batchSize = xCount(3);
%     xCount = xCount(length(xCount):-1:1);

    while xStart(1) <= totalSize
        [tSeries, xSeriesTmp] = readSimulationData(filename, dataset, xStart, xCount);
        if ~exist('xStddevSq')
            dim = size(xSeriesTmp);
            dim(1) = 1;
            xStddevSq = zeros(dim);
            xMean = repmat(xMean, [size(xSeriesTmp, 1), 1, 1]);
        end
        xStddevSq = xStddevSq + sum((xSeriesTmp - xMean) .^ 2, 1);
        xStart(1) = xStart(1) + batchSize;
    end
    xStddev = sqrt(1/(totalSize - 1) * xStddevSq);
    %xSeries = squeeze(xSeries);
end
