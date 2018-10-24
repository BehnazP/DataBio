function save_fig(gcf,filename,orientation)

%oriantation = {'landscape','portrait'}
    
h=gcf;
set(h,'PaperOrientation',orientation);
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', filename);
end