function MainIRMAD(image_t1,image_t2,band_name_t1,band_name_t2,opts,save_name,path_name)
% ChangeDetection :
% The function allows to perform the change detection analysis IRMAD:
%
% IRMAD is a technique that performs a Iteratively Reweighted Multivariate Alteration Detection. The chi
% square is calculated applying the EM algorithm, which estimate the multidimensional gaussian mixture of
% both changed and unchanged classes.
% 
% ---------------------------------
% Nicola Falco
% nicolafalco@ieee.org
% &
% Prashanth Reddy Marpu
% prashanthmarpu@ieee.org
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
%   MainIRMAD(image_t1,image_t2,band_name_t1,band_name_t2,opts,save_name,path_name);
%
%   opts: [epsln, size_im, flag_mask, low_val, flag_save]
%
% ---------------------------------
% Inputs:
%
%   - image_t1          - string of the whole path of the image at data t1
%   - image_t2          - string of the whole path of the image at data t2
%   - band_name_t1      - string of the band name of the image at t1
%   - band_name_t2      - string of the band name of the image at t2
%
%           ----------------------------------------------------------------
%   - opts for IRMAD: (vector of integer)
%
%       1   - epsln     - epsilon for iterations.
%                         default value = 0.01
%
%       2   - size_im   - in based on the dimension of the data set, is possible to select a
%                         different processing strategy:
%                         0 -> the whole data set is read and save in the workspace
%                         1 -> the data set is read line by line from files(less memory
%                         used but it takes more time)
%                         default value = 0
%
%       3   - flag_mask - 1 to choose a strict mask thresholding based on EM algorithm
%                         2 to choose a relaxed mask thresholding based on EM algorithm
%                         3 to use the principal component instead of the original data set
%                         0 no mask
%                         default value = 0
%
%       4   - low_val   - value in % to mask the low values of the histogram
%                         0 no mask
%                         default value = 0
%
%
%       5   - flag_save - 0 to save the chi square files
%                       - 1 to save the chi square and the intermediate files
%
%
%   - save_name         - string of the name used as suffix when the files are saved.
%                         insert '' for no one
%   - path_name         - string of the path where to save the files.
%                         insert '' to save in the current directory.
%---------------------------------
% Output:
%
%   - ICMmask           - initial change mask (pre-processing)
%
%   - cv1               - canonical variates1 unit variance
%   - cv2               - canonical variates2 unit variance
%   - chi2_irmad        - Chi Square obtained by using cvs
%
%
% ---------------------------------
% Reference:
%
%   N. Falco, P. R. Marpu, and J. A. Benediktsson,
%   "A Toolbox for Unsupervised Change Detection Analysis,"
%   International Journal of Remote Sensing 37 (7): 1505?26, 2016. Doi:10.1080/01431161.2016.1154226.
%
%   N. Falco, P. R. Marpu, and J. A. Benediktsson,
%   "Comparison of ITPCA and IRMAD for automatic change detection using  initial change mask,"
%   IEEE International Geoscience and Remote Sensing Symposium (IGARSS '12), 2012.
%
%   P. R. Marpu, P. Gamba, and M. J. Canty;
%   "Improving Change Detection Results of IR-MAD by Eliminating Strong Changes",
%   IEEE Geosci. Remote Sens. Lett., vol. 8, no. 4, pp. 799?803, July 2011.
%
%   M. J. Canty;
%   "Classification, and Change Detection in Remote Sensing, With Algorithms for ENVI/IDL"
%   Taylor and Francis, Second revised edition, 2010.
%
%   M. J. Canty and A. A. Nielsen.;
%   "Automatic radiometric normalization of multitemporal of satellite imagery
%    with the iteratively re-weighted MAD transformation"
%   Remote Sensing of Environment, 112(3):1025-1036, 2008.
%
%   A. A. Nielsen.;
%   "The regularized iteratively reweighed MAD method for change detection in
%    multi- and hyperspectral data"
%   IEEE Transcation on Image Processing, 16(2):463-478,2007.
%
%   R. Wiemker, A. Speck, D. Kulbach, H. Spitzer, and B. Johann;
%   "Unsupervised robust change detection on multispectral imagery using spectral and spatial features"
%   Proceedings of the Third International Airborne Remote Sensing Conference and Exhibition, vol. 1, no. July 1997, pp. 7?10, 1997.
%
% ---------------------------------

p = gcp('nocreate');
if isempty(p)
    parpool;
end 

disp('/////////////////////////////////////////////////');
disp('/////////////////////////////////////////////////');
disp('. . . changeDetection: FUNCTION IN PROGRESS . . .');
disp('/////////////////////////////////////////////////');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------- initializing the inputs ----------------------------%  
if ~iscell(image_t1)
    tmp = textscan(image_t1,'%s','Delimiter',',')';
    image_t1 = tmp{:}';
end

if ~iscell(image_t2)
    tmp = textscan(image_t2,'%s','Delimiter',',')';
    image_t2 = tmp{:}';
end

if ~iscell(band_name_t1)
    tmp = textscan(band_name_t1,'%s','Delimiter',',')';
    band_name_t1 = tmp{:}';
end

if ~iscell(band_name_t2)
    tmp = textscan(band_name_t2,'%s','Delimiter',',')';
    band_name_t2 = tmp{:}';
end

if nargin < 5
    opts = [0.01,0,0,0,0];
    path_name = fullfile(pwd,'/');
    save_name  = [];
else
    if ~isnumeric(opts)
        tmp = textscan(opts,'%f','Delimiter',',')';
        opts = tmp{1}';
    end

    if ~ischar(save_name)
        if isempty(save_name)
            save_name  = [];
        else
            tmp = textscan(save_name,'%s')';
            save_name = tmp{:}';
        end
    end

    if ~ischar(path_name)
        if isempty(path_name)
            path_name = fullfile(pwd,'/');
        else
            tmp = textscan(path_name,'%s')';
            path_name = tmp{:}';
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------- IRMAD Method ------------------------------------%    

epsln = opts(1);
if isequal(epsln,0)
    disp('epsln value no valid. Insert a value greater then 0');
    return;
end

size_im = opts(2);
if size_im > 1
    disp('size_im value no valid. Insert a value between 0 and 1');
    return;
end

flag_mask = opts(3);
if flag_mask > 3
    disp('flag_mask value no valid. Insert a value between 0 and 3');
    return;
end

low_val = opts(4);

flag_save = opts(5);
if (isequal(flag_save, 0))==0
    fl_save = 1;
else
    fl_save = 0;
end

opts_vec = [epsln,flag_mask,low_val,0,fl_save];

save_name   = num2str(save_name);
path_name   = num2str(path_name);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------- IRMAD Computation ----------------%
[~, ~, ext] = fileparts(image_t1{1});

if isequal(size_im ,0)          % ---------SAVE FILE IN WORKSPACE---------

    [ file_chi2, ~, ~,~,rho ] = IRMAD(image_t1,image_t2,band_name_t1,band_name_t2,opts_vec,save_name,path_name);
    
    %%%%%%%%%%%%%% added by Bezo for user %%%%%%%%%%%%%%
    figure(2)
    [nband, niter] = size(rho);
    for i = 1:nband
        plot(rho(i,:), '.-');
        ylabel('Canonical correlations');
        xlabel('iteration')
        xlim([0,niter+0.1])
        ylim([-0.1,1.1])
        hold on
    end
    title('Canonical correlations')

    if strcmp(ext,'.tif')
        delete(file_chi2); delete([file_chi2,'.hdr']);
    end

elseif isequal(size_im,1)       % ------------- LINE by LINE-------------

    [ file_chi2, ~, ~,~,rho ] = IRMADLine(image_t1,image_t2,band_name_t1,band_name_t2,opts_vec,save_name,path_name);
    figure(2)
    [nband, niter] = size(rho);
    for i = 1:nband
        plot(rho(i,:), '.-');
        ylabel('Canonical correlations');
        xlim([0,niter+0.1])
        ylim([-0.1,1.1])
        hold on
    end
    
    if strcmp(ext,'.tif')
        delete(file_chi2); delete([file_chi2,'.hdr']);
    end
end

disp('/////////////////////////////////////////////////');
disp('. . . changeDetection - IRMAD: PROCESS OVER . . .');
disp('/////////////////////////////////////////////////');

