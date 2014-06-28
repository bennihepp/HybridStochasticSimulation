function info = readSimulationInfo(filename, simulationName)
    if nargin < 1
        filename = selectFile;
    end
    if nargin < 2
            info = h5info(filename, '/simulations');
            groupNames = cell(length(info.Groups),1);
            for i=1:length(info.Groups)
                groupName = info.Groups(i).Name;
                parts = regexp(groupName, '/', 'split');
                groupNames{i,1} = char(parts(end));
            end
            [Selection,ok] = listdlg('ListString', groupNames, 'SelectionMode', 'single', ...
                'Name', 'Select Simulation', 'OKString', 'Select');
            if ok == 0
                clear info;
                return;
            end
            i = Selection(1);
            simulationName = info.Groups(i).Name;
            ['Selected simulation "', simulationName, '"']
    else
        simulationName = ['/simulations/', simulationName];
    end
    info = h5info(filename, simulationName);
end
