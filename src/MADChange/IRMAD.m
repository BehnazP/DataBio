function [ file_chi2,file_mads,file_cv1,file_cv2,rho_out] = IRMAD(image_t1,image_t2,band_name_t1,band_name_t2,opts,save_name,path_name)
% IRMAD: 
% Iteratively Reweighted Multivariate Alteration Detection
% ---------------------------------
%   IRMAD(image_t1,image_t2,band_name_t1,band_name_t2,opts,save_name,path_name);
%
%   opts = [epsln,flag_mask,low_val,w,flag_save] vector of integer
% ---------------------------------
% Inputs:
%
%   - image_t1              - string of the whole path of the image at data t1
%   - image_t2              - string of the whole path of the image at data t2
%   - band_name_t1          - name of the band exist in data t1
%   - band_name_t2          - name of the band exist in data t2
%
%   - opts:
%     
%     1  - epsln            - epsilon for iterations. Insert 100 to apply 0 iterations
%                             default value 0.01. 
%     2  - flag_mask        - 1 to choose the mask threshold between the 1st and 2nd distribution 
%                             2 to choose the mask threshold between the 2nd and 3rd distribution
%                             3 to use the principal component instead of the original data set
%                             0 no mask
%                             default value = 0
%     3  - low_val          - value in % to mask the low values of the histogram
%                             0 no mask
%                             default value = 0
%     4  - w                - weights for stats calculation (available only for command line)
%                             0 no weights
%                             default value = 0
%     5  - flag_save        - 0 to save the only chi square and MAD variables files
%                           - 1 to save the chi square and the intermediate files
%                             default value = 0
%
%   - save_name             - string of the name used as suffix when the files are saved (optional)                        
%   - path_name             - string of the path where to save the files (optional)
% ---------------------------------
% Otputs: 
%
%   - file_chi2             - string of the whole path of the chi2
%   - file_mads             - string of the whole path of the mads
%   - file_cv1              - string of the whole path of canonical variates cv1
%   - file_cv2              - string of the whole path of canonical variates cv2
%
%   (stored in ENVI bip format and if the original data is GeoTiff it stored in Geotiff)
%
%   - chi2                  - Chi Square
%   - mads                  - the MAD variates
%   - cv1                   - canonical variates1 unit variance
%   - cv2                   - canonical variates2 unit variance
% ---------------------------------
% 
% Reference to the method
%
% A. A. Nielsen.;
% "The regularized iteratively reweighed MAD method for change detection in multi- and hyperspectral data"
% IEEE Transcation on Image Processing, 16(2):463-478,2007.
% ---------------------------------
% original code by:
% Nicola Falco 
% nicolafalco@ieee.org
% 
% Prashanth Reddy Marpu
% prashanthmarpu@ieee.org
% 11/10/2015 last revision
% ---------------------------------
% Inspired by the code of:
% Allan Aasbjerg Nielsen, Ph.D., M.Sc.
% aa@space.dtu.dk, www.imm.dtu.dk/~aa
% 19 Sep 2010
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
disp('IRMAD: function in progress . . .');
disp('------------------------------------------------');

if isequal(nargin,7)
    [imaget1,imaget2,~,~] = read_optic_data(image_t1,image_t2,band_name_t1,band_name_t2);
else
    disp('Incorrect number of inputs');
    return;
end

[nrow,ncol,nband] = size(imaget1);

% read the input options vector

epsln = opts(1);
if isequal(epsln,0)
    disp('epsln value no valid. Insert a value greater then 0');
    return;
end

flag_mask = opts(2);
if flag_mask > 3
    disp('flag_mask value no valid. Insert a value between 0 and 3');
    return;
end

low_val = opts(3);

w = opts(4);
if w == 0
    w = ones(nrow*ncol,1);    
else
    [nroww,ncolw] = size(w);
    if ~(ncol==ncolw && nrow==nroww)
        disp('input1, input2 and w do not match');
        return
    end
end

flag_save = opts(5);

% application of the initial change mask
if flag_mask > 0
    disp('IRMAD: built of the Initial Change Mask');
    mask = ICM(image_t1,image_t2,band_name_t1,band_name_t2,flag_mask,low_val,save_name,path_name); 
end

F = reshape(imaget1,nrow*ncol,nband);
G = reshape(imaget2,nrow*ncol,nband);
% mask application
if ~isequal(flag_mask,0)
    mask = reshape(mask,nrow*ncol,1);
    ids = find(mask);
    Fmask = F(ids,:);
    Gmask = G(ids,:);
end

if strcmp(save_name,'') == 1
    file_mads = [path_name,'mads'];
    file_cv1  = [path_name,'cv1'];
    file_cv2  = [path_name,'cv2'];
    file_chi2 = [path_name,'chi2_irmad'];
else
    file_mads = [path_name,'mads_',save_name];
    file_cv1  = [path_name,'cv1_',save_name];
    file_cv2  = [path_name,'cv2_',save_name];
    file_chi2 = [path_name,'chi2_irmad_',save_name];
end

% rho initialization
prev_rho = 100*ones(1,nband);
rho = [];

MAXITER = 30;

%-----------------------------------------------------------------------------------------
%%% to see the progress added by Bezo
H = waitbar(0,'Process in progress ...');
rho_out = [];
tic;
for it = 1 : MAXITER
    
    if ~isequal(flag_mask,0)
        w = w(ids,:);
        [cov_mat,mean_vec] = covw([Fmask Gmask],w);
    else
        [cov_mat,mean_vec] = covw([F G],w);
    end
    
    %   generalized eigen value problem
    %
    %   [Mfg * inv(Mgg) * Mgf] * Va = Lamda * Mff * Va
    %   [Mgf * inv(Mff) * Mfg] * Vb = Lamda * Mgg * Vb
    %   M1 = [Mfg * inv(Mgg) * Mgf]
    
    Mff = cov_mat(1:nband,1:nband);
    Mgg = cov_mat(nband+1:end,nband+1:end);
    Mfg = cov_mat(1:nband,nband+1:end);
    Mgf = Mfg';
    
    meanF = mean_vec(1:nband);
    meanG = mean_vec(nband+1:end);
    
    % eigenvectors Va
    M1 = Mfg / Mgg * Mgf;
    [Va,da] = eigs(M1,Mff);
   
    varU = Va' * Mff * Va;
    
    % eigenvectors normalization to get unit variance CVs
    Va = bsxfun(@rdivide,Va,sqrt(diag(varU))');
    
    % sign correction of the correlations Va
    if ~isequal(flag_mask,0)
        invstd = diag(1./std(Fmask));
    else
        invstd = diag(1./std(F));
    end
    
    sign_diag = diag(sign(sum(invstd*Mff*Va)));
    Va = Va*sign_diag;
    
    % eigenvectors Vb from Va
    Vb = Mgg \ Mgf * Va;
    varV = Vb' * Mgg * Vb;
    
    % eigenvectors normalization to get unit variance CVs
    Vb = bsxfun(@rdivide,Vb,sqrt(diag(varV))');
    
    % canonical correlations rho
    rho = (diag(sqrt(da))');
    
    rho_out = [rho_out rho'];
    
    if it==1, disp('IRMAD: Canonical correlations'); end
    disp(num2str(rho,' %0.6g'));
    
    % variance of mads
    varMads = 2*(1-rho);        
    
    % irmad computation
    cv1 = bsxfun(@minus,F,meanF) * Va;
    cv2 = bsxfun(@minus,G,meanG) * Vb;
    Mads = cv1 - cv2; 
    Chi2 = sum(bsxfun(@rdivide,Mads.*Mads,varMads),2);
   
    % thresholding comparison
    if max(abs(rho - prev_rho)) < epsln
        waitbar(MAXITER,H)
        break;
    end
    
    prev_rho = rho;
    %w = 1-gammainc(0.5*Chi2(:),0.5*nband);
    % added by Bezo for numerically more accurate results, check help function
    w = gammainc(0.5*Chi2(:),0.5*nband,'upper');
    
    %%% added by Bezo
    %wait bar for showing progress
    waitbar(it/MAXITER,H)

end
close(H)

% image w as probability of No Change
if ~isequal(flag_mask,0)
    wM = zeros(nrow*ncol,1);
    wM(ids,:) = w;
    w = reshape(wM,[nrow,ncol]);
else
    w = reshape(w,[nrow,ncol]);
end
figure(1)
imagesc(w)
caxis([0,1])
axis off
colorbar
colormap gray
title('Probability of no change')
save_fig(gcf,[path_name,'ProbNoChange'],'landscape')

cv1 = reshape(cv1,nrow,ncol,nband);
cv2 = reshape(cv2,nrow,ncol,nband);

Chi2 = reshape(Chi2,nrow,ncol);
Mads = reshape(Mads,nrow,ncol,nband);

%%%% Data Saving %%%%
disp('IRMAD: saving file');

matlabToEnvi(Chi2,file_chi2,'bip');
%save(file_chi2, 'Chi2');
matlabToEnvi(Mads,file_mads,'bip');
%save(file_mads, 'Mads');

[~, ~, ext] = fileparts(image_t1{1});
if strcmp(ext,'.tif')
    info = geotiffinfo(image_t1{1});
    geotiffwrite(file_chi2,Chi2,info.SpatialRef,'CoordRefSysCode',info.GeoTIFFTags.GeoKeyDirectoryTag.ProjectedCSTypeGeoKey)
    geotiffwrite(file_mads,Mads,info.SpatialRef,'CoordRefSysCode',info.GeoTIFFTags.GeoKeyDirectoryTag.ProjectedCSTypeGeoKey)
    delete(file_mads);delete([file_mads,'.hdr']);
end

if isequal(flag_save,1)
    matlabToEnvi(cv1,file_cv1,'bip');
    %save(file_cv1, 'cv1');
    matlabToEnvi(cv2,file_cv2,'bip');
    %save(file_cv2, 'cv2');
    if strcmp(ext,'.tif')
        info = geotiffinfo(image_t1{1});
        geotiffwrite(file_cv1,cv1,info.SpatialRef,'CoordRefSysCode',info.GeoTIFFTags.GeoKeyDirectoryTag.ProjectedCSTypeGeoKey)
        geotiffwrite(file_cv2,cv2,info.SpatialRef,'CoordRefSysCode',info.GeoTIFFTags.GeoKeyDirectoryTag.ProjectedCSTypeGeoKey)
        delete(file_cv1); delete([file_cv1,'.hdr']);
        delete(file_cv2); delete([file_cv2,'.hdr']);
    end
end

disp('IRMAD: process over');
disp('------------------------------------------------');
disp('------------------------------------------------');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% PCeval Function %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [PComp1,PComp2] = PCeval( data_in1, data_in2 )

% reading of the data 
[nrow,nband] = size(data_in1);

data1 = NaN(nrow,nband);
data2 = NaN(nrow,nband);

% subtraction between each band and its own mean
for i=1 : nband
    data1(:,i) = data_in1(:,i) - mean(data_in1(:,i));
    data2(:,i) = data_in2(:,i) - mean(data_in2(:,i));
end

% covariance matrix
cov_matrix1 = cov(data1);
cov_matrix2 = cov(data2);

% PC analysis considering the PCs 
[coeff1, eigenvalues1] = pcacov(cov_matrix1);
[coeff2, eigenvalues2] = pcacov(cov_matrix2);

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
    n = 5;
elseif (n > 10) 
    n = 10;
end

new_coeff1 = coeff1(:,1:n);
new_coeff2 = coeff2(:,1:n);

PComp1=(data1*new_coeff1);
PComp2=(data2*new_coeff2);