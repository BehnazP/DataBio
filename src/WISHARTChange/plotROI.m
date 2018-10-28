function plotROI(fig,ROI)

figure(fig)
for i = 1:3
    subplot(1,3,i)
    hold on
    line(ROI.x, ROI.y, 'Color', 'r','LineWidth',1);
end

end