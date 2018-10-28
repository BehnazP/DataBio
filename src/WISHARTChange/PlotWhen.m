function  A = PlotWhen(fig,when,times,ROI)

% make the table of p-value of test statistics 

% input:
% pP		-		p-values from WiCHProb function
% times		-		number of time, 4th dimension of the data
% 
% output:
% a table for different mean  P-value

if isempty(fig)
	figure('Name','Results','NumberTitle','off');
end

cmap = parula(times);
cmap(1,:) = [0.8,0.8,0.8];
nsz = times;

% first change
[maxch,changefirst] = max(when,[],3);
changefirst(maxch == 0) = NaN;
subplot(131);
imagesc(changefirst)
colormap(gca,cmap(1:nsz,:));
caxis([0,nsz-1]);

if nsz > 6 
    a1 = colorbar;
    v = repmat(2:nsz-1,[2 1]);v = v(:)';
    labels = strsplit(sprintf('[t_{%d}, t_{%d}]\n',[1,v,times]),'\n');labels = labels([end,1:end-1]);labels{1}='NoChange';
    ticks = 0.5*nsz/(nsz+1):nsz/(nsz+1):nsz;
    indx = round(linspace(1,nsz,7));
    a1.Ticks = ticks(indx);
    a1.TickLabels = labels(indx);
else
    a1 = colorbar('eastoutside','YTick',0.5*nsz/(nsz+1):nsz/(nsz+1):nsz,'YTickLabel',int2str((1:nsz+1)'), 'YLim', [0, nsz-1]);
    v = repmat(2:nsz-1,[2 1]);v = v(:)';
    labels = strsplit(sprintf('[t_{%d}, t_{%d}]\n',[1,v,times]),'\n');labels = labels([end,1:end-1]);labels{1}='NoChange';
    a1.TickLabels = labels;
end

title('First change')
axis image
axis off
pos = get(gcf,'position');
set(gcf,'position',[10 pos(2) pos(3)*4 pos(4)]);

% most recent/last change
[maxfreq,changelast] = max(cumsum(when,3),[],3);
changelast(maxfreq == 0) = NaN;
subplot(132);
imagesc(changelast)
colormap(gca,cmap(1:nsz,:));
caxis([0,nsz-1]);
if nsz > 6 
    a3 = colorbar;
    a3.Ticks = ticks(indx);
    a3.TickLabels = labels(indx);
else
    a3 = colorbar('eastoutside','YTick',0.5*nsz/(nsz+1):nsz/(nsz+1):nsz,'YTickLabel',int2str((1:nsz+1)'), 'YLim', [0, nsz-1]);
    a3.TickLabels = labels;
end
title('Last change')
axis image
%set(gca,'xtick',[],'ytick',[]);			
axis off

% frequency of change
sumC = sum(when,3);
freq = sumC;
nsz = max(freq(:));
if nsz == 0
    nsz = 1;
end
cmap1 = parula(nsz+1);
cmap1(1,:) = [0.8,0.8,0.8];
freq(sumC==0) = 0;
subplot(133);
imagesc(freq)
colormap(gca,cmap1(1:nsz+1,:));
caxis([0, nsz]);
a2 = colorbar('eastoutside','YTick',0.5*nsz/(nsz+1):nsz/(nsz+1):nsz,'YTickLabel',int2str((1:nsz+1)'), 'YLim', [0, nsz]);
labels = strseq('',0:nsz)';labels{1}='NoChange';
a2.TickLabels = labels;
title('Frequency of change')
axis image
%set(gca,'xtick',[],'ytick',[]);
axis off

if nargin > 3
    plotROI(ROI);
end

% for saving give back 3layered image [firstchange,lastchange, frequency of change]
A = cat(3,changefirst,changelast,freq);
 
end