function [imaget1,imaget2,info1,info2] = read_optic_data(image_t1,image_t2,band_name_t1,band_name_t2)
% READ_OPTIC_DATA this function can read GEOTIFF and ENVI images
%
% input :
% image_t1         - a cell structure containing path of all the bands of image at time 1
% image_t2         - a cell structure containing path of all the bands of image at time 2
% band_name_t1     - a cell structure containing name of all the bands of image at time 1
% band_name_t2     - a cell structure containing name of all the bands of image at time 2
% 
% output :
% imaget1          - an image matrix at time 1 (nrow,ncol,nband)
% imaget2          - an image matrix at time 2 (nrow,ncol,nband)
% info1            - a structure depending on the image type at time 1; .hdr orgeotiff information
% info2            - a structure depending on the image type at time 2; .hdr orgeotiff information
% 
% Read_optic_data.m bepi@dtu.dk

if ischar(image_t1)
    [~, ~, ext1] = fileparts(image_t1);
    [~, ~, ext2] = fileparts(image_t2);
else
    [~, ~, ext1] = fileparts(image_t1{1});
    [~, ~, ext2] = fileparts(image_t2{1});
end
if ~isequal(ext1,ext2)
    msgbox('The datasets have different format','Error','error')
    error('The datasets have different format');
end

switch ext1
    case '.tif'
        [imaget1,imaget2,info1,info2] = read_geotiff_data(image_t1,image_t2,band_name_t1,band_name_t2);
    otherwise
        [imaget1,imaget2,info1,info2] = read_envi_data(image_t1,image_t2);
        if strcmp(info1.interleave, 'bip') == 0 || strcmp(info2.interleave, 'bip') == 0
            disp('The input images have to be in ENVI BIP format!')
            return
        end
end

end