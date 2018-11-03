function save_fig(gcf,filename,orientation)
% SAVE_FIG  save the figures in .pdf format
%
% save_fig(gcf,filename,orientation)
% filename          the name fo the file to be saved
% oriantation      'landscape' or 'portrait'
%
% SAVE_FIG by Behnaz Pirzamanbein bepi@dtu.dk

h=gcf;
set(h,'PaperOrientation',orientation);
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', filename);
end
