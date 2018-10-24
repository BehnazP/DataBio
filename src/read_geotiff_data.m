function [imaget1,imaget2,info1,info2] = read_geotiff_data(image_t1,image_t2,band_name_t1,band_name_t2)
if ischar(image_t1)
    sz = size(image_t1,1);
    if isempty(sz)
        msgbox('Exit without selecting any file','Error','error')
        error('Exit without selecting any file')
    end

    info1 = geotiffinfo(image_t1);
    imaget1 = imread(image_t1);
    info2 = geotiffinfo(image_t2);
    imaget2 = imread(image_t2);
else
    [imaget1, info1 ] = GeoTifRead(image_t1,band_name_t1);
    [imaget2, info2 ] = GeoTifRead(image_t2,band_name_t2);
end
end