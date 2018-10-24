function [lines,info,sizes] = read_optic_data_Line(image,band_name,varargin)
% READ_OPTIC_DATA_Line this function can read GEOTIFF and ENVI images
%
% input :
% image        - image path
% band_name    - a cell structure containing name of all the bands of image
% 
% output :
% line          - a vector containing a row of an image (1,ncol,nband)
% info          - a structure depending on the image type; .hdr orgeotiff information
% sizes         - [nrow,ncol,nband];
%
% Read_optic_data.m bepi@dtu.dk

if size(varargin) == 1 
    row = varargin{1};
else
    row = 1;
end

if ischar(image)
    [pathstr, name, ext] = fileparts(image);
    tmp = image ;
    sz = size(image,1);
else
    [pathstr, name, ext] = fileparts(image{1});
    tmp = image{1};
    sz = size(image,2);
end

switch ext
    case '.tif'
        info = geotiffinfo(tmp);
        nrow = info.Height;
        ncol = info.Width;
        nband = sz;
    otherwise
        fname = [name,'.hdr'];
        hdrfile = fullfile(pathstr,fname);
        info = envihdrread(hdrfile);
        nrow = info.lines;
        ncol = info.samples;
        [precision, machineformat] = envInfo(info);
        % check the format of the input images
        if strcmp(info.interleave, 'bip') == 0 
            disp('The input images have to be in ENVI BIP format!')
            return
        end
        nband = size(band_name,2);
end

if nargout > 2
    sizes = [nrow,ncol,nband];
    lines = 0;
    info = [];
    return
end

if sz == 1
    switch ext
        case '.tif'
            lines = imread(image,'PixelRegion',{[row row],[1,ncol]});
            lines = reshape(lines,[ncol,1]);
        otherwise
            fileINDX = fopen(image,'r');
            lines = fread(fileINDX,ncol*nband,precision, 0, machineformat);
            lines = reshape(lines,[ncol,nband]); 
            fclose(fileINDX);
    end
else
    lines = nan(1,ncol,nband);
    for i = 1:sz
        [~, name, ~] = fileparts(image{i});
        for j = 1: nband
            if contains(name,band_name{j})
                switch ext
                    case '.tif'
                        lines(1,:,j) = imread(image{i},'PixelRegion',{[row row],[1,ncol]});
                        break
                    otherwise
                        fileINDX = fopen(image{i},'r');
                        lines(1,:,j)= fread(fileINDX,ncol,precision, 0, machineformat);
                        fclose(fileINDX);
                        break
                end
            end
        end
    end
    lines = reshape(lines,[ncol,nband]);    
end

end