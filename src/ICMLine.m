function path_mask = ICMLine(image_t1,image_t2,band_name_t1,band_name_t2,flag_mask,low_val,save_name,path_name)
% ICMLine: function that identifies strong changes and creates a mask of them
% ---------------------------------
%   ICMLine(image_t1,image_t2,band_name_t1,band_name_t2,flag_mask,low_val,save_name,path_name)
% ---------------------------------
% Inputs:
%
%   - image_t1              - string of the whole path of the image at data t1
%   - image_t2              - string of the whole path of the image at data t2
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
%   - path_mask              - string of the whole path of the initial change mask
%                              1 -> unchange
%                              0 -> change
% ---------------------------------
% Reference:
%
%   Marpu, P. R. ; Gamba, P. and Canty, M. J. ;
%   "Improving Change Detection Results of IR-MAD by Eliminating Strong Changes", 
%   IEEE Geosci. Remote Sens. Lett., vol. 8, no. 4, July 2011.
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
%
% ---------------------------------


disp('------------------------------------------------');
disp('------------------------------------------------');
disp('ICMLine: function in progress . . .');
disp('------------------------------------------------');

if isequal(nargin,8)
    [~,~,sizes] = read_optic_data_Line(image_t1,band_name_t1);
    nrow = sizes(1);
    ncol = sizes(2);
    nband = sizes(3);
else
    disp('Incorrect number of inputs');
    return;
end

if flag_mask < 3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% Equalization Image %%%%%%%%%%%%%%%%
    
    disp('ICMLine: processing of the original dataset');

    disp('ICMLine: stretching of the first image');
    path_tmp1 = [path_name,'str_im1'];
    LStretch(image_t1,band_name_t1,path_tmp1);
    disp('ICMLine: stretching of the second image');
    path_tmp2 = [path_name,'str_im2'];
    LStretch(image_t2,band_name_t2,path_tmp2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% Random Samples Selection from max_diff %%%%%%
    %[max_diff, rand_data] = maxDiff(tmp1,tmp2, path_name, save_name);
    
    [hdr_in, precision_tmp, machineformat_tmp] = envihdrread([path_tmp1,'.hdr']);

    dim = 50000;
    v(:,1) = randi([1 nrow],1,dim);
    v(:,2) = randi([1 ncol],1,dim);
    [idx,IX] = sort(v(:,1));
    new = [v(IX,1),v(IX,2)];
    
    rand_data = [];
    fileTMP1 = fopen(path_tmp1,'r');
    fileTMP2 = fopen(path_tmp2,'r');
    path_diff = [path_name,'max_diff'];
    fileDIFF = fopen(path_diff,'w');
    for r = 1 : nrow
        line1 = fread(fileTMP1, ncol * nband, precision_tmp, 0, machineformat_tmp);
        line2 = fread(fileTMP2, ncol * nband, precision_tmp, 0, machineformat_tmp);
        
        % differential of the recalibrated images
        diff = reshape((line1 - line2),nband, ncol);
        
        % find the max difference
        if nband > 1
            max_diff = max(abs(diff));
        elseif nband == 1
            max_diff = abs(diff);
        end
        fwrite(fileDIFF,max_diff,class(max_diff));
        
        % select random data
        rand_data = [rand_data ,max_diff(new(idx==r,2))];
    end
    hdrWrite(path_diff,nrow,ncol,1,class(max_diff));
    fclose(fileTMP1);
    fclose(fileTMP2);
    fclose(fileDIFF);

    clear line1 line2 max_diff str_im1 str_im2
    delete(path_tmp1); delete([path_tmp1,'.hdr']); delete(path_tmp2); delete([path_tmp2,'.hdr']);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%% EM_GM Algorithm %%%%%%%%%%%%%%%%%
    
    disp('ICMLine: gaussian distributions estimation');
    [W,M,V] = EM_GM(rand_data(:),3,[],2,[],[]);
    
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
    
    A1 = log((std1 / std2)*(w2/w1));
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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% Principal Components %%%%%%%%%%%%%%
    
    disp('ICMLine: processing of the principal components');

    path_diff = [path_name,'diff'];
    path_PC = [path_name,'PC'];
    
    
    fileDIFF = fopen(path_diff,'w');
    for r = 1 : nrow
        line1 = read_optic_data_Line(image_t1,band_name_t1);
        line2 = read_optic_data_Line(image_t2,band_name_t2);
        diff = reshape((line1 - line2),1,nband*ncol);
        fwrite(fileDIFF,diff,class(diff));
    end
    fclose(fileDIFF);
    hdrWrite(path_diff,nrow,ncol,nband,class(diff));
    
    [cov_mat,mean_vec] = covmatEval(path_diff,band_name_t1);
    [coeff, ~] = pcacov(cov_mat);
    
    nPCs = 1;
    new_coeff = coeff(:,1:nPCs);
    
    hdrDIFF = envihdrread([path_diff,'.hdr']);
    [precisionDIFF, machineformatDIFF] = envInfo(hdrDIFF);
    fileDIFF = fopen(path_diff, 'r');
    filePC = fopen(path_PC,'w');
    for r = 1 : nrow
        %%%%%%%%%%%%%%% the first line of the b-th band
        line = fread(fileDIFF, ncol * nband, precisionDIFF, 0, machineformatDIFF);
        line = reshape(line, nband, ncol);
        
        %%%%%%%%%%%%%%% mean subtraction
        line = bsxfun(@minus, line,mean_vec');
        
        %%%%%%%%%%%%%%% principal components
        PC = new_coeff'*line;
        
        PCs = reshape(PC,1,ncol*nPCs);
        fwrite(filePC,PCs,class(PCs));
    end
    fclose(fileDIFF);
    fclose(filePC);
    hdrWrite(path_PC,nrow,ncol,nPCs,class(PCs));
    
    dim = 50000;
    v(:,1) = randi([1 nrow],1,dim);
    v(:,2) = randi([1 ncol],1,dim);
    [idx,IX] = sort(v(:,1));
    new = [v(IX,1),v(IX,2)];
    
    for b = 1 : nPCs
        filePC = fopen(path_PC,'r');
        [~,hdrPC] = enviread(path_PC,[path_PC,'.hdr']);
        [precisionPC, machineformatPC] = envInfo(hdrPC);
        rand_data = [];
        for r = 1 : nrow
            line = fread(filePC, nPCs*ncol, precisionPC, 0, machineformatPC);
            line = reshape(line, nPCs, ncol);
            rand_data = [rand_data ,line(b,new(idx==r,2))];
        end
        fclose(filePC);
        clear line
        
        % EM_GM Algorithm %
        fprintf('ICMLine: gaussian distributions estimation of the component %d\n',b);
        [~,M,~] = EM_GM(rand_data(:),3,[],2,[],[]);
        
        [M, ~]=sort(M);
        m1(b) = M(1);
        m3(b) = M(3);
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Low Values %%%%%%%%%%%%%%%%%%%%

if low_val > 0
    path_maskLow1 = [path_name,'Low1'];
    maskLow(image_t1,low_val,path_maskLow1,band_name_t1);
    
    path_maskLow2 = [path_name,'Low2'];
    maskLow(image_t2,low_val,path_maskLow2,band_name_t2);
    
    fileMASKL1 = fopen(path_maskLow1,'r');
    fileMASKL2 = fopen(path_maskLow2,'r');
    
    [hdrml, precision_ml, machineformat_ml] = envihdrread([path_maskLow1,'.hdr']);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Thresholding  %%%%%%%%%%%%%%%%%%
disp('ICMLine: thresholding');

if strcmp(save_name,'') == 1
    path_mask = [path_name,'ICMask'];
else
    path_mask = [path_name,'ICMask_',save_name];
end

if flag_mask < 3
    
    [hdrDIFF, precision_DIFF, machineformat_DIFF] = envihdrread([path_diff,'.hdr']);
    
    fileDIFF = fopen(path_diff,'r');
    fileMASK = fopen(path_mask,'w');
    for r = 1 : nrow
        % the first line of the b-th band
        lineDIFF = fread(fileDIFF, ncol, precision_DIFF, 0, machineformat_DIFF);
        
        [line1,~] = read_optic_data_Line(image_t1,band_name_t1);
        line1 = permute(line1,[2,1]);
        [line2,~] = read_optic_data_Line(image_t2,band_name_t2);
        line2 = permute(line2,[2,1]);
        
        ICMask = zeros(1,ncol);
        
        if flag_mask == 1
            ICMask(lineDIFF >= St1) = 1;
        elseif flag_mask == 2
            ICMask(lineDIFF >= St2) = 1;
        end
        
        for b = 1 : nband
            ICMask(line1(b,:) == 0) = 1;
            ICMask(line2(b,:) == 0) = 1;
        end
        
        if low_val > 0
            lineM1 = fread(fileMASKL1, ncol, precision_ml, 0, machineformat_ml)';
            lineM2 = fread(fileMASKL2, ncol, precision_ml, 0, machineformat_ml)';
            vecLine = (lineM1 & lineM2);
            ICMask = ( vecLine | ICMask );
        end
        
        ICMask = cast(ICMask,'double');
        ICMask(ICMask==0)=2;
        ICMask(ICMask==1)=0;
        ICMask(ICMask==2)=1;
        fwrite(fileMASK,ICMask,class(ICMask));
    end
    hdrWrite(path_mask,nrow,ncol,1,class(ICMask));
    fclose(fileDIFF);
    fclose(fileMASK);
    
    clear lineDIFF line1 line2;
    delete(path_diff); delete([path_diff,'.hdr']);
    
elseif flag_mask == 3
    
    Sta = max(m1);
    Stb = min(m3);
    min_val = min(Sta,Stb);
    max_val = max(Sta,Stb);
    
    [hdrPC, precisionPC, machineformatPC] = envihdrread([path_PC,'.hdr']);
    
    filePC = fopen(path_PC, 'r');
    fileMASK = fopen(path_mask,'w');
    
    for r = 1 : nrow
        %%%%%%%%%%%%%%% the first line of the b-th band
        linePC = fread(filePC, ncol * nPCs, precisionPC, 0, machineformatPC);
        linePC = reshape(linePC, nPCs, ncol);
        
        ICMask = zeros(nPCs,ncol);
        ICMask((linePC < min_val | linePC > max_val))=1;
        
        ICMask = sum(ICMask,1);
        ICMask(ICMask < nPCs) = 0;
        ICMask(ICMask ~= 0) = 1;
        
        [line1,~] = read_optic_data_Line(image_t1,band_name_t1);
        line1 = permute(line1,[2,1]);
        [line2,~] = read_optic_data_Line(image_t2,band_name_t2);
        line2 = permute(line2,[2,1]);
        
        for b = 1 : nband
            ICMask(line1(b,:) == 0) = 1;
            ICMask(line2(b,:) == 0) = 1;
        end
        
        if low_val > 0
            lineM1 = fread(fileMASKL1, ncol, precision_ml, 0, machineformat_ml)';
            lineM2 = fread(fileMASKL2, ncol, precision_ml, 0, machineformat_ml)';
            vecLine = (lineM1 & lineM2);
            ICMask = ( vecLine | ICMask );
        end
        %%%%%%%%%%%%%%%%%%%%%% Bezo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %hdrwrite just write in single or double so using uint8 doesnt
        %change anything
        ICMask = cast(ICMask,'double');
        ICMask(ICMask==0)=2;
        ICMask(ICMask==1)=0;
        ICMask(ICMask==2)=1;
        fwrite(fileMASK,ICMask,class(ICMask));
    end
    hdrWrite(path_mask,nrow,ncol,1,class(ICMask));
    fclose(fileMASK);
    fclose(filePC);    
    delete(path_PC); delete([path_PC,'.hdr']);
    
end
if low_val > 0
    fclose(fileMASKL1);
    fclose(fileMASKL2);
    delete(path_maskLow1); delete([path_maskLow1,'.hdr']); 
    delete(path_maskLow2); delete([path_maskLow2,'.hdr']); 
end

mask = enviread(path_mask);
save(path_mask, 'mask');
[~, ~, ext] = fileparts(image_t1{1});
if strcmp(ext,'.tif')
    info = geotiffinfo(image_t1{1});
    geotiffwrite(path_mask,mask,info.SpatialRef,'CoordRefSysCode',info.GeoTIFFTags.GeoKeyDirectoryTag.ProjectedCSTypeGeoKey)
end
clear mask

disp('ICMLine: process over');
disp('------------------------------------------------');
disp('------------------------------------------------');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% LStretch Function %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function LStretch (path_in,band_name,path_name)

[~,~,sizes] = read_optic_data_Line(path_in,band_name);
nrow = sizes(1);
ncol = sizes(2);
nband = sizes(3);

%%%%%%%%%%%%%%%% histogram
max_val = 0;
min_val = 0;

for r = 1 : nrow
    [line,~] = read_optic_data_Line(path_in,band_name);

    max_val_new = max(line(:)); %max value over all the bands
    max_val = max(max_val,max_val_new);
    
    min_val_new = min(line(:)); %min value over all the bands
    min_val = min(min_val,min_val_new);
end
clear line

rng = (max_val-min_val)/256;
x = min_val+rng/2:rng:max_val-rng/2;
elem = zeros(256,nband);

for r = 1 : nrow
    [line,~] = read_optic_data_Line(path_in,band_name);

    for b = 1 : nband
        elem_new(:,b) = hist(line(:,b),x);
    end
    elem = elem + elem_new;
end

clear line

%%%%%%%%%%%%%%%% cumulative histogram
cumHist = (cumsum(elem))/(nrow*ncol);
for b = 1 : nband
    if isempty(x(find(cumHist(:,b) < 0.02,1,'last')))
        Gmin(b) = x(find(min(cumHist(:,b))));
    else
        Gmin(b) = x(find(cumHist(:,b) < 0.02,1,'last')) ;
    end
    %Gmax(b) = x(find(cumHist(:,b) > 0.98,1,'first'));
    if isempty(x(find(cumHist(:,b) > 0.98,1,'first')))
        Gmax(b) = x(find(max(cumHist(:,b))));
    else
        Gmax(b) = x(find(cumHist(:,b) > 0.98,1,'first'));
    end
end
gmin = 0;
gmax = 255;

%%%%%%%%%%%%%%%% linear stretching

fileSTR = fopen(path_name,'w');
for r = 1 : nrow
    [line,~] = read_optic_data_Line(path_in,band_name,r);

    for b = 1 : nband
        im = line(:,b);
        
        for i = 1 : size(im,1)
            Pin = im(i);
            if (Pin <= Gmin(b))
                str_image(b,i) = gmin; % 0
            elseif ((Pin > Gmin(b)) && (Pin < Gmax(b)))
                str_image(b,i) = ((Pin - Gmin(b)) * (gmax - gmin) / (Gmax(b) - Gmin(b)));
            elseif Pin >= Gmax(b)
                str_image(b,i) = gmax; % 255
            end
        end
    end
    str = reshape(str_image,1,nband*ncol);
    fwrite(fileSTR,str,class(str));
end
fclose(fileSTR);
hdrWrite(path_name,nrow,ncol,nband,class(str));
clear line


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  maskLow Function %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function maskLow (path_in,low_val,path_name,band_name)

[~,~,sizes] = read_optic_data_Line(path_in,band_name);

nrow = sizes(1);
ncol = sizes(2);
nband = sizes(3);

%%%%%%%%%%%%%%%% histogram

max_val=0;
min_val=0;
for r = 1 : nrow
    [line,~] = read_optic_data_Line(path_in,band_name);
    
    max_val_new = max(line(:)); %max value over all the bands
    max_val = max(max_val,max_val_new);
 
    min_val_new = min(line(:)); %min value over all the bands
    min_val = min(min_val,min_val_new);
end
clear line

rng = (max_val-min_val)/256;
x = min_val+rng/2:rng:max_val-rng/2;
elem = zeros(256,nband);

for r = 1 : nrow
    [line,~] = read_optic_data_Line(path_in,band_name,r);
    
    for b = 1 : nband
        elem_new(:,b) = hist(line(:,b),x);
    end
    elem = elem + elem_new;
end

clear line

%%%%%%%%%%%%%%%% cumulative histogram
cumHist = (cumsum(elem))/(nrow*ncol);
for b = 1 : nband
    Gmin(b) = x(find(cumHist(:,b) < (low_val/100),1,'last')) ;
end

%%%%%%%%%%%%%%%% linear stretching

fileMASKL = fopen(path_name,'w');

for r = 1 : nrow
    [line,~] = read_optic_data_Line(path_in,band_name,r);
    line = permute(line,[2,1]);
    MaskLow = ones(1,ncol);
    
    for b = 1 : nband
        MaskLow(line(b,:) > Gmin(b)) = 0;
    end
    fwrite(fileMASKL,MaskLow,class(MaskLow));
end
fclose(fileMASKL);
hdrWrite(path_name,nrow,ncol,1,class(MaskLow));

clear line