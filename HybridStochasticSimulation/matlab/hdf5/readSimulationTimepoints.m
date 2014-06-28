function tSeries = readSimulationTimepoints(filename, dataset, tStart, tCount)
    if nargin < 1
        filename = selectFile;
    end
    if nargin < 2
        info = readSimulationInfo(filename);
        parts = regexp(info.Name, '/', 'split');
        dataset = char(parts(end));
    end
    tInfo = h5info(filename, ['/simulations/', dataset, '/tSeries']);
    tSize = tInfo.Dataspace.Size;
    if nargin < 3
        tStart = 1;
    end
    if nargin < 4
        tCount = tSize - tStart;
    end
    tSeries = h5read(filename, ['/simulations/', dataset, '/tSeries'], tStart, tCount);
end
