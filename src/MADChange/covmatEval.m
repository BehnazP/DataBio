function [cov_matrix,mean_vec] = covmatEval(varargin)

% Function that performs the covariance matrix of an updating data set.
% It is based on the provisional mean algorithm. The algorithm reads 
% each images line by line from files.
% ---------------------------------
% Syntax:
%   
%   covmatEval(image_t1,band_name_t1)
%   covmatEval(image_t1,band_name_t1,image_t2,band_name_t2)
%   covmatEval(image_t1,band_name_t1,image_t2,band_name_t2,weight,mask)
% ---------------------------------
% Input:
%   
%   -image_t1       string of the whole path of the image at data t1 (ENVI bip format)
%   -band_name_t1
%   -image_t2       string of the whole path of the image at data t2 (ENVI bip format)
%   -band_name_t2
%   -weight         initial weight matrix ( '' for none)
%   -mask           string of the whole path of the mask ( '' for none) (ENVI bip format)
% ---------------------------------
% Output:
%
%   -cov_matrix     covariance matrix
%   -mean_vec       mean vector 
% ---------------------------------  
% Nicola Falco 
% nicolafalco@ieee.org
% 
% 11/09/2011 first version
% 15/10/2015 last version
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

if size(varargin) < 2
    disp('COVMat must have at least 2 inputs: image_t1!')
    return
    
elseif size(varargin,2) == 2
    fl_im = 1;
    image_t1     = varargin{1};
    band_name_t1 = varargin{2};
    [~,~,sizes]  = read_optic_data_Line(image_t1,band_name_t1);
    weight       = ones(sizes(2));
    mask         = 0;
    
elseif size(varargin,2) == 4
    if iscellstr(varargin(5)) == 1
        fl_im = 2;
        image_t1     = varargin{1};
        band_name_t1 = varargin{2};
        image_t2     = varargin{3};
        band_name_t2 = varargin{4};
        [~,~,sizes]  = read_optic_data_Line(image_t1,band_name_t1);
        weight       = ones(sizes(2));
        mask         = 0;
    else 
        fl_im = 1;
        image_t1     = varargin{1};
        band_name_t1 = varargin{2};
        weight       = varargin{3};
        mask         = 0;
    end
    
elseif size(varargin,2) == 5
    if iscellstr(varargin(5)) == 1
        fl_im = 2;
        image_t1     = varargin{1};
        band_name_t1 = varargin{2};
        image_t2     = varargin{3};
        band_name_t2 = varargin{4};
        weight       = varargin{5};
        mask         = 0;
    else 
        fl_im = 1;
        image_t1     = varargin{1};
        band_name_t1 = varargin{2};
        weight       = varargin{3};
        mask         = varargin{4};
    end
    
elseif size(varargin,2) == 6
    fl_im = 2;
    image_t1     = varargin{1};
    band_name_t1 = varargin{2};
    image_t2     = varargin{3};
    band_name_t2 = varargin{4};
    weight       = varargin{5};
    mask         = varargin{6};
end

if ~isequal(mask,0)
    hdrmask = envihdrread([mask,'.hdr']);
    [precisionMask, machineformatMask] = envInfo(hdrmask);
end
if strcmp(weight,'') == 1
    [~,~,sizes] = read_optic_data_Line(input,band_name);
    weight = ones(sizes(2));
end


if fl_im == 1
    [~,~,sizes] = read_optic_data_Line(image_t1,band_name_t1);
    NN = sizes(3);
    sw = 0;
    n = 0;
    mn = zeros(1,NN);
    cov_matrix = zeros(1,NN^2);
elseif fl_im == 2
    [~,~,sizes] = read_optic_data_Line(image_t1,band_name_t1);
    NN = 2*sizes(3);
    sw = 0;
    n = 0;
    mn = zeros(1,NN);
    cov_matrix = zeros(1,NN^2);
end


if isequal(mask,0) % without mask
    for r = 1 : sizes(1) % rows
        
        % data to add
        if fl_im == 1
            [line1,~] = read_optic_data_Line(image_t1,band_name_t1,r);
            data = permute(line1,[2,1]);
            
        elseif fl_im == 2
            [line1,~] = read_optic_data_Line(image_t1,band_name_t1,r);
            [line2,~] = read_optic_data_Line(image_t2,band_name_t2,r);
            data = [permute(line1,[2,1]); permute(line2,[2,1])];
          
        end
        Ws = weight(r,:);
        
        % loop over all observation
        for i = 1 : sizes(2)%size(data,2)
            sw = sw + Ws(i);
            c = Ws(i)/sw;
            %if 1st Ws = 0 then c  = NaN!
            if isnan(c) == 1
                c = 0;
            end
            
            % mean
            for j = 1 : NN
%                 d(j) = (data((i-1)*NN+j)) - mn(j) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%% bezo to get the correct indexing %%%%%%%%%%%%
                d(j) = data(j,i) - mn(j) ;
                mn(j) = mn(j) + d(j)*c;
            end
            
            % weighted covariance
            for j = 1 : NN
                for k =1 : NN
                    cov_matrix((j-1)*NN+k) = cov_matrix((j-1)*NN+k) + d(j)*d(k)*(1-c)*Ws(i);
                end
            end
            
        end
%        n = n + sum(Ws);
        n = n + sum(Ws~=0);
    end
    
elseif ~isequal(mask,0) % with mask
    
    fileMASK = fopen(mask, 'r');
    for r = 1 : sizes(1) % rows
        
        % data to add
        mask_im = fread(fileMASK, hdrmask.samples, precisionMask, 0, machineformatMask)';
        if fl_im == 1
            [line1,~] = read_optic_data_Line(image_t1,band_name_t1,r);
            
            % masking
            data = permute(line1,[2,1]);
            ids = find(mask_im);
            data_mask = data(:,ids);
           
        elseif fl_im == 2
            
            [line1,~] = read_optic_data_Line(image_t1,band_name_t1,r);
            [line2,~] = read_optic_data_Line(image_t2,band_name_t2,r);
            
            % masking
            data = [permute(line1,[2,1]); permute(line2,[2,1])];
            ids = find(mask_im);
            data_mask = data(:,ids);
        
        end
        Ws = weight(r,ids);
        
        % loop over all observation
        for i = 1 : size(data_mask,2)
            sw = sw + Ws(i);
            c = Ws(i)/sw;
            
            % mean
            for j = 1 : NN
                d(j) = data_mask(j,i) - mn(j) ;
                mn(j) = mn(j) + d(j)*c;
            end
            
            % weighted covariance
            for j = 1 : NN
                for k =1 : NN
                    cov_matrix((j-1)*NN+k) = cov_matrix((j-1)*NN+k) + d(j)*d(k)*(1-c)*Ws(i);
                end
            end
            
        end
%        n = n + sum(Ws);
        n = n + sum(Ws~=0);
    end
    fclose(fileMASK);
end
%cov_matrix = reshape(cov_matrix/(n-1),NN,NN);
cov_matrix = reshape(cov_matrix * n/(sw*(n-1)),NN,NN);
mean_vec = mn;

clear data_mask data mask_im line1 line2  
end
