function saveWhen(image,save_name,files)
% SAVEWHEN  save the image of PlotWhen in .pdf format
%
% saveWhen(image,save_name,files)
% image        -  a matrix with 3 bands: [first change,last change,frequency of the change]
%                 this is the results of function PlotWhen
% save_name    -  the name to save the images
% files      -  a cell structure containing all the images path and name
%               Example: in windows
%               {'C:\Download\201712.VV_0'},{'C:\Download\201712.VV_1'}
%               {'C:\Download\201712.VH_0'},{'C:\Download\201712.VV_1'} and
%               so on, the name SHOULD contain .'pol_name'_'time_name'.format%
% SAVEWHEN by Behnaz Pirzamanbein bepi@dtu.dk, last version 2018-11-03

[~, ~, ext] = fileparts(files{1});

switch ext
    case '.tif'
        info = geotiffinfo(files{1});
        geotiffwrite(save_name,image,info.SpatialRef,'CoordRefSysCode', ...
        info.GeoTIFFTags.GeoKeyDirectoryTag.ProjectedCSTypeGeoKey)
    otherwise
        nrow = size(image,1);
        ncol = size(image,2);
        hdrWrite(save_name,nrow,ncol,3,class(image),'bsq',1,0,0)
        A = fopen(save_name,'w');
        fwrite(A,image,class(image));
        fclose(A);
end
end
