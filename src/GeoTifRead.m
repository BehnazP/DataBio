function [images, info ] = GeoTifRead(files,band_name)
sz = size(files,2);

if isempty(sz)
    msgbox('Exit without selecting any file','Error','error')
    error('Exit without selecting any file')
end

info = geotiffinfo(files{1});
nrow = info.Height;
ncol = info.Width;
nband = size(band_name,2);

images = nan(nrow,ncol,nband);
h = waitbar(0,'Processing...','Name','Reading Data');
for i = 1:sz
    [~, name, ~] = fileparts(files{i});
    for j = 1:nband
        if contains(name,band_name{j})
        images(:,:,j) = imread(files{i});
        end
    end
    waitbar(i/sz,h)
end
close(h)
end