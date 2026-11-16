function d = data_dir()
%DATA_DIR Folder holding the input CSV files.
d = fullfile(fileparts(mfilename('fullpath')), '..', 'data');
end
