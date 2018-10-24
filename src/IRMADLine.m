function [path_chi2,path_mads,path_cv1,path_cv2,rho_out] = IRMADLine(image_t1,image_t2,band_name_t1,band_name_t2,opts,save_name,path_name)
% IRMADLine:
% Iteratively Re-weighted Multivariate Alteration Detection Line by line
% ---------------------------------
%   IRMADLine(image_t1,image_t2,band_name_t1,band_name_t2,opts,save_name,path_name);
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
%                             default value 0.05. 
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
% Outputs: 
%
%   - path_chi2             - string of the whole path of the chi2
%   - path_mads             - string of the whole path of the mads
%   - path_cv1              - string of the whole path of canonical variates cv1
%   - path_cv2              - string of the whole path of canonical variates cv2
%
%   (stored in ENVI bip format and if the orginal data is GeoTiff it also stored in Geotiff)
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
% "The regularized iteratively reweighed MAD method for change detection in multi- and hyper-spectral data"
% IEEE Transaction on Image Processing, 16(2):463-478,2007.
%----------------------------------
% original code by:
% Nicola Falco 
% nicolafalco@ieee.org
% 
% Prashanth Reddy Marpu
% prashanthmarpu@ieee.org
% 11/10/2015 last revision
%
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
disp('IRMADLine: function in progress . . .');
disp('------------------------------------------------');

if isequal(nargin,7)
    % read the inputs
    [~,~,sizes] = read_optic_data_Line(image_t1,band_name_t1);
    nrow = sizes(1);
    ncol = sizes(2);
    nband = sizes(3);
else
    disp('Incorrect number of inputs');
    return;
end

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
    w = ones(nrow,ncol);    
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
    disp('IRMADLine: built of the Initial Change Mask');
    flag_mask = ICMLine(image_t1,image_t2,band_name_t1,band_name_t2,flag_mask,low_val,save_name,path_name); 
end

% rho initialization
prev_rho = 100*ones(nband,1);
rho = [];
F = image_t1;
G = image_t2;

MAXITER = 30;

if strcmp(save_name,'') == 1
    path_mads = [path_name,'mads'];
    path_cv1 = [path_name,'cv1'];
    path_cv2 = [path_name,'cv2'];
    path_chi2 = [path_name,'chi2_irmad'];
else
    path_mads = [path_name,'mads_',save_name];
    path_cv1 = [path_name,'cv1_',save_name];
    path_cv2 = [path_name,'cv2_',save_name];
    path_chi2 = [path_name,'chi2_irmad_',save_name];
end

%-----------------------------------------------------------------------------------------
%added by Bezo 
H = waitbar(0,'Process in progress ...');
rho_out = [];
for it = 1 : MAXITER
    [cov_mat,mean_vec] = covmatEval(F,band_name_t1,G,band_name_t2,w,flag_mask);

    %   generalized eigen value problem
    %
    %   [Mfg * inv(Mgg) * Mgf] * Va = Lamda * Mff * Va
    %   [Mgf * inv(Mff) * Mfg] * Vb = Lamda * Mgg * Vb
    %   M1 = [Mfg * inv(Mgg) * Mgf]
    
    Mff = cov_mat(1:nband,1:nband);
    Mgg = cov_mat(nband+1:end,nband+1:end);
    Mfg = cov_mat(1:nband,nband+1:end);
    Mgf = Mfg';
    
    meanF = (mean_vec(1:nband))';
    meanG = (mean_vec(nband+1:end))';
    
    % eigenvectors Va
    M1 = Mfg / Mgg * Mgf;
    [Va,da] = eigs(M1,Mff);
        
    varU = Va' * Mff * Va;
    
    % eigenvectors normalization to get unit variance CVs
    Va = bsxfun(@rdivide,Va,sqrt(diag(varU))');
    
    % sign correction of the correlations Va
    sign_diag = diag(sign(sum(bsxfun(@rdivide,Mff*Va,stdEval(F,band_name_t1,flag_mask)))));
    Va = Va*sign_diag;
    
    % eigenvectors Vb from Va
    Vb = Mgg \ Mgf * Va;
    varV = Vb' * Mgg * Vb;
    
    % eigenvectors normalization to get unit variance CVs
    Vb = bsxfun(@rdivide,Vb,sqrt(diag(varV))');

    % canonical correlations rho
    rho = (diag(sqrt(da)));
    
    rho_out = [rho_out rho];
    
    if it==1, disp('IRMADLine: Canonical correlations'); end
    disp(num2str(rho',' %0.6g'));
    
    % variance of mads
    varMads = 2*(1-rho);
    
    % thresholding comparison
    if max(abs(rho - prev_rho)) < epsln
        clear Chi2;
        
        if ~isequal(flag_save,0)
            CV1 = fopen(path_cv1,'w');
            CV2 = fopen(path_cv2,'w');
        end
        MAD = fopen(path_mads,'w');
        CHI = fopen(path_chi2,'w');        
                
        disp('IRMADLine: saving file');
        for r = 1 : nrow
            [lineF,~] = read_optic_data_Line(F,band_name_t1,r);
            lineF = permute(lineF,[2,1]);
            [lineG,~] = read_optic_data_Line(G,band_name_t2,r);
            lineG = permute(lineG,[2,1]);
            
            cv1 = Va' * bsxfun(@minus,lineF,meanF);
            cv2 = Vb' * bsxfun(@minus,lineG,meanG);
            Mads = cv1 - cv2;
            Chi2 = sum(bsxfun(@rdivide,Mads.*Mads,varMads),1);
            
            if ~isequal(flag_save,0)
                fwrite(CV1,reshape(cv1,1,nband*ncol),class(cv1));
                fwrite(CV2,reshape(cv2,1,nband*ncol),class(cv2));
            end
            fwrite(CHI,Chi2,class(Chi2));
            fwrite(MAD,reshape(Mads, 1,nband*ncol),class(Mads));

        end
        fclose(CHI);
        fclose(MAD);
        if ~isequal(flag_save,0)
            fclose(CV1);
            fclose(CV2);
        end
        waitbar(MAXITER,H)
        break;
    end

    % irmad computation
    parfor r= 1 : nrow
        [lineF,~] = read_optic_data_Line(F,band_name_t1,r);
        lineF = permute(lineF,[2,1]);
        [lineG,~] = read_optic_data_Line(G,band_name_t2,r);
        lineG = permute(lineG,[2,1]);
        
        Mi = Va' * bsxfun(@minus,lineF,meanF) - Vb' * bsxfun(@minus,lineG,meanG);
        Chi2(r,:) = sum(bsxfun(@rdivide,Mi.*Mi,varMads),1);
    end
    
    % set the rho and weights for the forward iteration 
    prev_rho = rho;
    %w = 1-gammainc(0.5*Chi2(:),0.5*nband);
    % added by Bezo for numerically more accurate results, check help function
    w = gammainc(0.5*Chi2(:),0.5*nband,'upper');
    w = reshape(w,nrow,ncol);
    waitbar(it/MAXITER,H)
end
close(H)

% image w as probability of No Change
figure(1)
imagesc(w)
caxis([0,1])
axis off
colorbar
colormap gray
title('Probability of no change')
save_fig(gcf,[path_name,'ProbNoChange'],'landscape')

%%%% Data Saving %%%%
hdrWrite(path_chi2,nrow,ncol,1,class(Chi2));
Chi2 = enviread(path_chi2);
%save(path_chi2, 'Chi2');

hdrWrite(path_mads,nrow,ncol,1,class(Mads));
Mads = enviread(path_mads);
%save(path_mads, 'Mads');

[~, ~, ext] = fileparts(image_t1{1});
if strcmp(ext,'.tif')
    info = geotiffinfo(image_t1{1});
    geotiffwrite(path_chi2,Chi2,info.SpatialRef,'CoordRefSysCode',info.GeoTIFFTags.GeoKeyDirectoryTag.ProjectedCSTypeGeoKey)
    geotiffwrite(path_mads,Mads,info.SpatialRef,'CoordRefSysCode',info.GeoTIFFTags.GeoKeyDirectoryTag.ProjectedCSTypeGeoKey)
    delete(path_mads);delete([path_mads,'.hdr']);
end

if isequal(flag_save,1)
    hdrWrite(path_cv1,nrow,ncol,1,class(cv1));
    CV1 = enviread(path_cv1);
    %save(path_cv1, 'CV1');

    hdrWrite(path_cv2,nrow,ncol,1,class(cv2));
    CV2 = enviread(path_cv2);
    %save(path_mads, 'CV2');

    if strcmp(ext,'.tif')
        info = geotiffinfo(image_t1{1});
        geotiffwrite(path_cv1,CV1,info.SpatialRef,'CoordRefSysCode',info.GeoTIFFTags.GeoKeyDirectoryTag.ProjectedCSTypeGeoKey)
        geotiffwrite(path_cv2,CV2,info.SpatialRef,'CoordRefSysCode',info.GeoTIFFTags.GeoKeyDirectoryTag.ProjectedCSTypeGeoKey)
        delete(path_cv1); delete([path_cv1,'.hdr']);
        delete(path_cv2); delete([path_cv2,'.hdr']);
    end
end

disp('IRMADLine: process over');
disp('------------------------------------------------');
disp('------------------------------------------------');