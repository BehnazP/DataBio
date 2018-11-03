function [images, info ] = read_data(files,pol_name,time_name)
%
% READ_DATA read different images formats: GeoTIFF and ENVI binary image coupled
%           with a header file and it has a special module for big data
% Input:
% files      -  a cell structure containing all the images path and name
%               Example: in windows
%               {'C:\Download\201712.VV_0'},{'C:\Download\201712.VV_1'}
%               {'C:\Download\201712.VH_0'},{'C:\Download\201712.VV_1'} and
%               so on, the name SHOULD contain .'pol_name'_'time_name'.format
%               series of k 2x2 or 3x3 complex covariance matrices
%               (polarimetric SAR images in covariance matrix form,
%               hhhh (real), hhhv (complex), hhvv (complex),
%                            hvhv (real),    hvvv (complex),
%                                            vvvv (real)
% pol_name   -   a char vector containing the name of polarization appear
%               after the name of image following a dot
%               Example: .VV .VH
% times_name -  a vector containing either a sequence of number 0:10 or a
%               personalize numbers such as 062, 064, 067 at the end of the
%               name of the files
% Output:
% images     -  a 4D matrix containing images (NxM) for P different polarization
%                and T different time points; images (NxMxPxT)
% READ_DATA by Behnaz Pirzamanbein bepi@dtu.dk, last version 2018-11-03

sz = size(files,2);
if isempty(sz)
    msgbox('Exit without selecting any file','Error','error')
    error('Exit without selecting any file')
end

sz_time = size(time_name,2);
sz_pol = size(pol_name,2);

% if sz ~= sz_pol*sz_time
%     error('The dimention of data miss matched')
% end

[pathstr, name, ext] = fileparts(files{1});
switch ext
    case '.tif'
        info = geotiffinfo(files{1});
        nrow = info.Height;
        ncol = info.Width;
    otherwise
        fname = [name,'.hdr'];
        hdrfile = fullfile(pathstr,fname);
        info = envihdrread(hdrfile);
        nrow = info.lines;
        ncol = info.samples;
end

images = nan(nrow,ncol,sz_pol,sz_time);
h = waitbar(0,'Processing...','Name','Reading Data');
for i = 1:sz
    [pathstr, name, ext] = fileparts(files{i});
    fname = [name,'.hdr'];
    hdrfile = fullfile(pathstr,fname);
    for j = 1:sz_pol
        if contains(name,pol_name{j})
            for k = 1:sz_time
                if endsWith(name,['_',time_name{k}])
                    switch ext
                        case '.tif'
                            info = geotiffinfo(files{i});
                            images(:,:,j,k) = imread(files{i});
                            break
                        case '.emi'
                            [images(:,:,j,k), info] = emiread(files{i},hdrfile);
                            break
                        otherwise
                            [images(:,:,j,k), info] = enviread(files{i},hdrfile);
                            break
                    end
                end
            end
        end
    end
    waitbar(i/sz,h)
end
close(h)
end
