function save_fig(fig, name)
%SAVE_FIG Save a figure as PDF and colour EPS in ../figures.
here = fileparts(mfilename('fullpath')); out = fullfile(here, '..', 'figures');
if ~exist(out, 'dir'), mkdir(out); end
set(fig, 'PaperPositionMode', 'auto');
print(fig, fullfile(out, [name '.pdf']), '-dpdf', '-bestfit');   % PDF
print(fig, fullfile(out, [name '.eps']), '-depsc');              % colour EPS (vector)
end
