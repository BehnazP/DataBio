function ROI = readROI(shape_file,fig)
if ischar(shape_file)
    S = shaperead(shape_file);

    %creat ROI
    ROI.x = (S.X-R.XWorldLimits(1))/10-0.5;
    ROI.y = (R.YWorldLimits(2)-S.Y)/10-0.5;

    ROI.mask = roipoly(sz(1),sz(2),x(1:end-1),y(1:end-1));
elseif nargin > 1
    msgbox({'From the Result figure choose the ROI by:','1- Move the mouse to the right subplot, i.e "Frequency of change"','2- Creat your ROI','3- Rigth click inside the ROI','4- Choose "creat mask"','5- Click OK'},'ROI','help')
    figure(fig);
    [ROI.mask,ROI.x,ROI.y] = roipoly;
end
end