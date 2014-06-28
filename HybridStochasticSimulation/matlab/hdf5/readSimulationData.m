function [tSeries, xSeries] = readSimulationData(filename, dataset, xStart, xCount)
    if nargin < 1
        filename = selectFile;
    end
    if nargin < 2
        info = readSimulationInfo(filename);
        parts = regexp(info.Name, '/', 'split');
        dataset = char(parts(end));
    end
    %tInfo = h5info(filename, ['/simulations/', dataset, '/tSeries']);
    %timepoints = tInfo.Dataspace.Size;
    xInfo = h5info(filename, ['/simulations/', dataset, '/xSeries']);
    xSize = xInfo.Dataspace.Size;
    if nargin < 3
        %start = 1;
        xStart = ones(1, length(xSize));
    elseif length(xStart) == 1
        tmp = ones(1, length(xSize));
        tmp(length(tmp)) = xStart;
        xStart = tmp;
    else
        xStart = xStart(length(xStart):-1:1);
    end
    %xStart = ones(1, length(xSize));
    %xStart(length(xStart)) = start;
    if nargin < 4
        %count = xSize(length(xSize));
        xCount = xSize - xStart;
    elseif length(xCount) == 1
        tmp = xSize;
        tmp(length(tmp)) = xCount;
        xCount = tmp;
    else
        xCount = xCount(length(xCount):-1:1);
    end
    %xCount = xSize;
    %xCount(length(xSize)) = count;
    tStart = xStart(1);
    tCount = xCount(1);
    tSeries = h5read(filename, ['/simulations/', dataset, '/tSeries'], tStart, tCount);
    xSeries = h5read(filename, ['/simulations/', dataset, '/xSeries'], xStart, xCount);
    xSeries = permute(xSeries, fliplr(1:ndims(xSeries)));
    if ndims(xSeries) < 3
        dims = size(xSeries);
        dims = horzcat(ones(3 - length(dims)), dims);
        xSeries = reshape(xSeries, dims);
    end
end
