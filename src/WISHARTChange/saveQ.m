function saveQ(wc,save_name,files,ROI)
% SAVEQ  save the image MxNxtime-1 of Q; takes struct, wc, from WishartChange 
%           function as input and outputs (1-P)-values, i.e., no-change probabilities 
%           for region of interest (ROI)
%
% saveWhen(image,save_name,files)
% Input
% wc            -	is a struct from WiCHParallel function based on WishartChange function
% save_name     -   the name to save the images
% files         -   a cell structure containing all the images path and name
%                       Example: in windows
%                       {'C:\Download\201712.VV_0'},{'C:\Download\201712.VV_1'}
%                       {'C:\Download\201712.VH_0'},{'C:\Download\201712.VV_1'} and
%                   so on, the name SHOULD contain .'pol_name'_'time_name'.format
% ROI           -   is a region of interest
% SAVEQ by Behnaz Pirzamanbein bepi@dtu.dk, first version 2019-11-26

if ~isstruct(wc)
    error('wc must be struct as output from WiCH')
end

[n,m,ntime] = size(wc(1).P);

if nargin < 4
    ROI = true(size(wc(end).P(:,:,1)));
end

image = nan(n,m,ntime-1);

% no-change probability
for count = 1:(ntime - 1)
    tmp = 1 - wc(count).P(:,:,1);
    image(:,:,count) = tmp.*ROI;
end

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
