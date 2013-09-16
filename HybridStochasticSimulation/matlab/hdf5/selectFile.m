function filename = selectFile
        [filename, pathname] = uigetfile({'*.h5', 'Select HDF5 file'; '*.*', 'All files'});
        if isequal(filename, 0)
                clear filename;
            return;
        end
        filename = [pathname, filename];
        ['Selected file"', filename, '"']
end
