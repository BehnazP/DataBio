function ICMask = ICM(image_t1,image_t2,band_name_t1,band_name_t2,flag_mask,low_val,save_name,path_name)
% ICM: function that identifies strong changes and creates a mask of them
% ---------------------------------
%   ICM(image_t1,image_t2,band_name_t1,band_name_t2,flag_mask,low_val,save_name,path_name)
% ---------------------------------
% Inputs:
%
%   - image_t1              - string of the whole path of the image at data t1
%   - image_t2              - string of the whole path of the image at data t1
%   - band_name_t1          - name of the band exist in data t1
%   - band_name_t2          - name of the band exist in data t2
%   - flag_mask             - 1 to choose a strict mask thresholding based on EM algorithm
%                             2 to choose a relaxed mask thresholding based on EM algorithm
%                             3 to use the principal component instead of the original data set
%                             default value = 1
%   - low_val               - value in % to mask the low values of the histogram
%                             0 no mask
%                             defoult value = 0
%   - save_name             - string of the name used as suffix when the files are saved 
%                             insert '' or anything for no one or just skip it
%   - path_name             - string of the path where to save the files
%                             insert '' or anything to save in the current directory  
% ---------------------------------
% Otputs:
%
%   - ICMask                - initial change mask
%                             1 -> unchange
%                             0 -> change
% ---------------------------------
% Reference:
%
%   Marpu, P. R. ; Gamba, P. and Canty, M. J. ;
%   "Improving Change Detection Results of IR-MAD by Eliminating Strong Changes", 
%   IEEE Geosci. Remote Sens. Lett., vol. 8, no. 4, pp. 799?803, July 2011.
% --------------------------------- 
% original code by:
% Nicola Falco 
% nicolafalco@ieee.org
% 
% Prashanth Reddy Marpu
% prashanthmarpu@ieee.org
% 11/10/2015 last revision
% ---------------------------------
% Modified by :
% Behnaz Pirzamanbein
% bepi@dtu.dk
% Image Analysis and Computer Graphics section
% Department of Applied Mathematics and Computer Science
% Technical University of Denmark
% First version 09.01.2018
% it is now possible to use it in appdesginer
%----------------------------------

disp('------------------------------------------------');
disp('------------------------------------------------');
disp('ICM: function in progress . . .');
disp('------------------------------------------------');

if isequal(nargin,8)
    %%%%  Data Reading  %%%%%
    [imaget1,imaget2,~,~] = read_optic_data(image_t1,image_t2,band_name_t1,band_name_t2);
    rows = size(imaget1,1);
    cols = size(imaget1,2);
    nbands = size(imaget1,3);
else 
    disp('Incorrect number of inputs');
    return;
end

if flag_mask < 3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% Equalization Image %%%%%%%%%%%%%%%%
    
    disp('ICM: processing of the original dataset');

    fprintf('ICM: stretching of the first image');
    str_img1 = LStretch(imaget1);
 
    fprintf('\nICM: stretching of the second image');
    str_img2 = LStretch(imaget2);
    
    str_img1 = reshape(str_img1, rows*cols,nbands);
    str_img2 = reshape(str_img2, rows*cols,nbands);
    
    imaget1 = reshape(imaget1, rows*cols,nbands);
    imaget2 = reshape(imaget2, rows*cols,nbands);
    fprintf('\n');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% Random Samples Selection from max_diff %%%%%%
    diff_img = double(str_img1 - str_img2);
    
    if size(diff_img,2) > 1
        max_diff = (max(abs(diff_img')))';  % maximum changes
    else
        max_diff = abs(diff_img);
    end
    
    dim = 50000;
    data = zeros(dim,1);
    idx = randi([1 size(max_diff,1)],1,dim);
    for i = 1 : size(idx,2)
        data(i,1) = max_diff(idx(i));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%% EM_GM Algorithm %%%%%%%%%%%%%%%%%
    disp('ICM: gaussian distributions estimation');
    [W,M,V] = EM_GM(data(:),3,[],2,[],[]);
    
    [M, indi]=sort(M);
    V=V(indi);
    W=W(indi);
    m1 = M(1);
    m2 = M(2);
    m3 = M(3);
    
    std1 = sqrt(V(1));
    std2 = sqrt(V(2));
    std3 = sqrt(V(3));
    w1   = W(1);
    w2   = W(2);
    w3   = W(3);
    
    A1 = log((std1 ./ std2)*(w2/w1));
    sqrt1 = sqrt((m1-m2).^2 + ((2.*A1) .* (std1.^2 - std2.^2)));
    s1 = (((m2 .* std1.^2) - (m1 .* std2.^2) +  ((std1 .* std2)) .* sqrt1)) ./ (std1.^2 - std2.^2);
    s2 = (((m2 .* std1.^2) - (m1 .* std2.^2) -  ((std1 .* std2)) .* sqrt1)) ./ (std1.^2 - std2.^2);
    St1 = ((m1>m2).*(m1>s1).*(m2<s1).*s1) + ((m1>m2).*(m1>s2).*(m2<s2).*s2) + ...
        ((m2>m1).*(m2>s1).*(s1>m1).*s1) + ((m2>m1).*(m2>s2).*(s2>m1).*s2);
    
    A2 = log((std2 ./ std3)*(w3/w2));
    sqrt2 = sqrt((m2-m3).^2 + ((2.*A2) .* (std2.^2 - std3.^2)));
    s2 = (((m3 .* std2.^2) - (m2 .* std3.^2) +  ((std2 .* std3)) .* sqrt2)) ./ (std2.^2 - std3.^2);
    s3 = (((m3 .* std2.^2) - (m2 .* std3.^2) -  ((std2 .* std3)) .* sqrt2)) ./ (std2.^2 - std3.^2);
    St2 = ((m2>m3).*(m2>s2).*(m3<s2).*s2) + ((m2>m3).*(m2>s3).*(m3<s3).*s3) + ...
        ((m3>m2).*(m3>s2).*(s2>m2).*s2) + ((m3>m2).*(m3>s3).*(s3>m2).*s3);
    
elseif flag_mask == 3
    
    imaget1 = reshape(imaget1, rows*cols,nbands);
    imaget2 = reshape(imaget2, rows*cols,nbands);

    disp('ICM: processing of the principal components');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% Random Samples Selection from max_diff %%%%%%
    %diff_img = zeros(size(imaget1));
    diff_img = double(imaget1 - imaget2);
    
    nPCs = 1;
    pcs = PCeval(diff_img, nPCs);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%% EM_GM Algorithm %%%%%%%%%%%%%%%%%
    disp('ICM: gaussian distributions estimation');
    
    idx = randi([1 rows*cols],1,50000);
    % Random Samples Selection from max_diff %%%%%
    for b = 1 : nPCs
        fprintf('\ncomponent: %d ',b);
        a = pcs(:,b);
        a = a(:);
        a = a(idx);
        
        % EM_GM Algorithm %
        fprintf(' EM_GM Algorithm in progress \n');
        %figure(b);
        [~,M,~] = EM_GM(a,3,[],2,[],[]);
        
        [M, ~]=sort(M);
        m1(b) = M(1);
        m3(b) = M(3); 
    end
    clear a diff_img 
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Mask of the Stronger Changes %%%%%%%%%
fprintf('\n');
disp('ICM: thresholding');

if flag_mask < 3
    ICMask = zeros(size(max_diff,1),size(max_diff,2));
    
    if flag_mask == 1
        ICMask(max_diff >= St1) = 1;
    elseif flag_mask == 2
        ICMask(max_diff >= St2) = 1;
    end
    
elseif flag_mask == 3
    Sta = max(m1);
    Stb = min(m3);
    min_val = min(Sta,Stb);
    max_val = max(Sta,Stb);
    ICMask = zeros(size(pcs,1),size(pcs,2));
    ICMask(pcs < min_val | pcs > max_val)=1;
    ICMask = sum(ICMask,2);
    ICMask(ICMask < nPCs) = 0;
    ICMask(ICMask ~= 0) = 1;
end

%%%% add the zero value pixels to the mask %%%%
for b = 1 : nbands
    ICMask(imaget1(:,b) == 0) = 1;
    ICMask(imaget2(:,b) == 0) = 1;
end

%%% mask of the lowest values Changes %%%%
if low_val > 0
    maskLow1 = maskLow_fun (imaget1, low_val);
    maskLow2 = maskLow_fun (imaget2, low_val);
    vecLine = (maskLow1 & maskLow2);
    ICMask = ( vecLine | ICMask );
end

ICMask = cast(ICMask,'double');
ICMask(ICMask==0)=2;
ICMask(ICMask==1)=0;
ICMask(ICMask==2)=1;
   
ICMask = reshape(ICMask,rows,cols);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Data Saving %%%%%%%%%%%%%%%%%%%
disp('ICM: saving file');
if strcmp(save_name,'') == 1
    file_mask = [path_name,'ICMask'];
else
    file_mask = [path_name,'ICMask_',save_name];
end
matlabToEnvi(ICMask,file_mask,'bip');   
save([file_mask,'.mat'],'ICMask');

[~, ~, ext] = fileparts(image_t1{1});
if strcmp(ext,'.tif')
    info = geotiffinfo(image_t1{1});
    geotiffwrite(file_mask,ICMask,info.SpatialRef,'CoordRefSysCode',info.GeoTIFFTags.GeoKeyDirectoryTag.ProjectedCSTypeGeoKey)
    delete(file_mask);delete([file_mask,'.hdr']);
end

disp('ICM: process over');
disp('------------------------------------------------');
disp('------------------------------------------------');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% LStretch Function %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str_image=LStretch (input)

[nrow, ncol ,nband] = size(input);

%%%%%%%%%%%%%%%% histogram
max_val = max(input(:));
min_val = min(input(:));

rng = (max_val-min_val)/256;
x = min_val+rng/2:rng:max_val-rng/2;
elem = zeros(256,nband);

for b = 1 : nband
    inp=input(:,:,b);
    [elem(:,b), binc(:,b)] = hist(inp(:),x);
end

%%%%%%%%%%%%%%%% cumulative histogram
cumHist = (cumsum(elem))/(nrow*ncol);
for b = 1 : nband
    if isempty(x(find(cumHist(:,b) < 0.02,1,'last')))
        Gmin(b) = x(find(min(cumHist(:,b))));
    else
        Gmin(b) = x(find(cumHist(:,b) < 0.02,1,'last')) ;
    end
    Gmax(b) = x(find(cumHist(:,b) > 0.98,1,'first'));
end
gmin = 0;
gmax = 255;
str_image = zeros(nrow,ncol,nband);
%%%%%%%%%%%%%%%% linear stretching
for r = 1 : nrow
    line = input(r,:,:);
    line = reshape(line, ncol, nband );
    for b = 1 : nband
        im = line(:,b);
        
        for i = 1 : ncol
            Pin = im(i);
            if (Pin <= Gmin(b))
                str_image(r,i,b) = gmin; % 0
            elseif ((Pin > Gmin(b)) && (Pin < Gmax(b)))
                str_image(r,i,b) = ((Pin - Gmin(b)) * (gmax - gmin) / (Gmax(b) - Gmin(b)));
            elseif Pin >= Gmax(b)
                str_image(r,i,b) = gmax; % 255
            end
        end
    end
    
end
clear line


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% PCeval Function %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PCs = PCeval( data,nPCs )

% reading of the data 
features     = size(data,2);

% subtraction between each band and its own mean
for i=1: features
    data(:,i) = data(:,i) - mean(data(:,i));
end

% covariance matrix
cov_matrix = cov(data);

% PC analysis considering the PCs 
[coeff, ~] = pcacov(cov_matrix);

new_coeff = coeff(:,1:nPCs);
% w_eigen = eigenvalues/sum(eigenvalues);
% for i=1: size(coeff,2)
%     new_coeff(:,i) = coeff(:,i);
%     if sum(w_eigen(1:i)) > 0.99
%         break;
%     end
% end

PCs=(data*new_coeff);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MaskLow Function %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function maskLow = maskLow_fun (input, low_val)

maskLow = ones(size(input(:,1)));
for b = 1 : size(input,2)

    im = input(:,b);
    [elem, binc] = hist(im,256);
    cum = (cumsum(elem));
    cum=cum/size(im,1);

    Gmin = binc(find(cum < (low_val/100),1,'last'));
    
    maskLow(im > Gmin) = 0;

end