function setup_paths()
%SETUP_PATHS Add all package folders to the MATLAB path.
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'core'));
d = dir(fullfile(root, 'sec*'));
for i = 1:numel(d), addpath(fullfile(root, d(i).name)); end
set(0, 'DefaultFigureVisible', 'off');
end
