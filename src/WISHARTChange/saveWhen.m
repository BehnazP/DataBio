function saveWhen(image,save_name,files)

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