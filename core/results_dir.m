function d = results_dir()
%RESULTS_DIR Folder for numerical output (CSV/MAT/TXT).
here = fileparts(mfilename('fullpath')); d = fullfile(here, '..', 'results');
if ~exist(d, 'dir'), mkdir(d); end
end
