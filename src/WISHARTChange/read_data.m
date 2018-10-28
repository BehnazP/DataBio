function [images, info ] = read_data(files,pol_name,time_name)
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