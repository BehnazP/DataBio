function [imaget1,imaget2,hdr1,hdr2] = read_envi_data(image_t1,image_t2)

%%%%%%%%%%%%%%%%%% read ENVI %%%%%%%%%%%%%%%%%%%%%%%%%%%
hdr1=envihdrread([image_t1,'.hdr']);
hdr2=envihdrread([image_t2,'.hdr']);

% check BIP format
if strcmp(hdr1.interleave, 'bip') == 0
    fprintf('data 1 in BIP format is required\n')
    image_old = enviread(image_t1);
    image_t1 = [image_t1,'_bip'];
    matlabToEnvi(image_old, image_t1,'bip')
    hdr1 = envihdrread([image_t1,'.hdr']);
end

if strcmp(hdr2.interleave, 'bip') == 0
    fprintf('data 2 in BIP format is required\n')
    image_old = enviread(image_t2);
    image_t2 = [image_t2,'_bip'];
    matlabToEnvi(image_old, image_t2,'bip')
    hdr2 = envihdrread([image_t2,'.hdr']);
end
% check the data dimension
if isequal(hdr1.bands,hdr2.bands) == 0
    disp('\nThe datasets have different dimensions');
    return
end

imaget1 = envidataread(image_t1);
imaget2 = envidataread(image_t2);
end