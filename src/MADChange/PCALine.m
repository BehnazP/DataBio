function [ PC1,PC2 ] = PCALine( image_t1, image_t2 ,band_name_1,band_name_2)
% Perform the PCA in line (reading the image line by line) 
%
% Nicola Falco 
% nicolafalco@ieee.org
% 
% Signal Processing Lab, University of Iceland
% 11/10/2014
% --------------------------------- 
% Modified by :
% Behnaz Pirzamanbein
% bepi@dtu.dk
% Image Analysis and Computer Graphics section
% Department of Applied Mathematics and Computer Science
% Technical University of Denmark
% First version 09.01.2018
% it is now possible to use it in appdesginer
% ---------------------------------  

[~,~,sizes] = read_optic_data_Line(image_t1,band_name_1);
nrow = sizes(1);
ncol = sizes(2);
nband = sizes(3);
size_n = ncol*nband;


[cov_mat,mean_vec] = covmatEval(image_t1,band_name_1, image_t2,band_name_2);
cmat1 = cov_mat(1:nband,1:nband);
cmat2 = cov_mat(nband+1:end,nband+1:end);

[coeff1, eigenvalues1] = pcacov(cmat1);
[coeff2, eigenvalues2] = pcacov(cmat2);

avar1 = (mean_vec(1:nband));
avar2 = (mean_vec(nband+1:end));

%new_coeff = coeff(:,1:3);
w_eigen1 = eigenvalues1/sum(eigenvalues1);
w_eigen2 = eigenvalues2/sum(eigenvalues2);

for n1=1: size(coeff1,2)
    %new_coeff(:,i) = coeff(:,i);
    if sum(w_eigen1(1:n1)) > 0.99
        break;
    end
end

for n2=1: size(coeff2,2)
    %new_coeff(:,i) = coeff(:,i);
    if sum(w_eigen2(1:n2)) > 0.99
        break;
    end
end
n = max(n1,n2);

if (n < 5)
    if nband < 5
        n = nband;
    else
        n = 5;
    end
elseif (n > 10) 
    n = 10;
end

new_coeff1 = coeff1(:,1:n);
new_coeff2 = coeff2(:,1:n);

PC1='PC1';
PC2='PC2';

filePC1 = fopen(PC1,'w');
filePC2 = fopen(PC2,'w');

for r = 1 : nrow
    %%%%%%%%%%%%%%% the first line of the b-th band
    [line1,~] = read_optic_data_Line(image_t1,band_name_1,r);
    line1 = permute(line1,[2,1]);
    [line2,~] = read_optic_data_Line(image_t2,band_name_2,r);
    line2 = permute(line2,[2,1]);
    
    %%%%%%%%%%%%%%% mean subtraction
    line1 = bsxfun(@minus, line1,avar1');
    line2 = bsxfun(@minus, line2,avar2');
    
    line = line - repmat(mean_vec',1,ncol);
    line = bsxfun(@minus, line,mean_vec');

    
    %%%%%%%%%%%%%%% principal components
    Pc1 = new_coeff1' * line1;
    Pc2 = new_coeff2' * line2;
    
    % bip
    Pc1 = reshape(Pc1,1,ncol*n);
    fwrite(filePC1,Pc1,class(Pc1));
    
    Pc2 = reshape(Pc2,1,ncol*n);
    fwrite(filePC2,Pc2,class(Pc2));
end
fclose(filePC1);
fclose(filePC2);

hdrWrite(PC1,nrow,ncol,n,class(Pc1));
hdrWrite(PC2,nrow,ncol,n,class(Pc2));

[~, ~, ext] = fileparts(image);
if strcmp(ext,'.tif')
    pc1 = enviread(PC1);
    pc2 = enviread(PC2);
    info = geotiffinfo(image);
    geotiffwrite(PC1,pc1,info.SpatialRef,'CoordRefSysCode',info.GeoTIFFTags.GeoKeyDirectoryTag.ProjectedCSTypeGeoKey)
    geotiffwrite(PC2,pc2,info.SpatialRef,'CoordRefSysCode',info.GeoTIFFTags.GeoKeyDirectoryTag.ProjectedCSTypeGeoKey)
    delete(PC1);delete([PC1,'.hdr']); 
    delete(PC2);delete([PC2,'.hdr']);
    clear pc1
    clear pc2
end
end

