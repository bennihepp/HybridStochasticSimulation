function [tSeries, xSeries] = readSimulationData(filename, dataset)
    if nargin < 1
        filename = selectFile;
    end
    if nargin < 2
        info = readSimulationInfo(filename);
        parts = regexp(info.Name, '/', 'split');
        dataset = char(parts(end));
    end
    tSeries = h5read(filename, ['/simulations/', dataset, '/tSeries']);
    xSeries = h5read(filename, ['/simulations/', dataset, '/xSeries']);
    xSeries = permute(xSeries, fliplr(1:ndims(xSeries)));
    if ndims(xSeries) < 3
        dims = ndims(xSeries);
        dims = horzcat(ones(3 - length(dims)), dims);
        xSeries = reshape(xSeries, dims);
    end
end
