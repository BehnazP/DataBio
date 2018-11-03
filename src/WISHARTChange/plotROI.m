function plotROI(fig,ROI)
% PLOTROI plot the region of interest defined by ROI in the current figure figure
% plotROI(fig,ROI)
% PLOTROI by Behnaz Pirzamanbein bepi@dtu.dk, last version 2018-11-03

figure(fig)
for i = 1:3
    subplot(1,3,i)
    hold on
    line(ROI.x, ROI.y, 'Color', 'r','LineWidth',1);
end

end
