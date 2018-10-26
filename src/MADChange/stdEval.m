function [ std_vec ] = stdEval( image,band_name,varargin )
%STDEVAL: standard deviation
%
% input
%   image
%   band_name
%   mask
%
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
% Nicola Falco 
% nicolafalco@ieee.org
% 
% Signal Processing Lab, University of Iceland
% 11/10/2014
% ---------------------------------  

if nargin < 2
    disp('stdEval must have at least 2 inputs: image_t1!')
    return
    
elseif nargin == 2
    mask = 0; 
    
elseif nargin > 2
    mask = varargin{1}; 
    if ~isequal(mask,0)
        [hdrmask, precisionMask, machineformatMask]=envihdrread([mask,'.hdr']);
    end
end

% [hdr, precision, machineformat] = envihdrread([input,'.hdr']);
[~,~,sizes] = read_optic_data_Line(image,band_name);

mean_vec = meanEval(image,band_name,mask);
amount = 0;

if isequal(mask,0) % without mask
    for r = 1 : sizes(1) % rows
        
        % data to add
        [line,~] = read_optic_data_Line(image,band_name,r);
        line = permute(line,[2,1]);
        
        for b = 1 : sizes(3)
            data(b,:) = line(b,:) - mean_vec(b);
        end
        data = data.*data;
        amount = amount + sum(data,2);
    end
    std_vec = sqrt(amount/((sizes(2)*sizes(1))-1));
    
elseif ~isequal(mask,0) % with mask
    
    fileMASK = fopen(mask, 'r');
    for r = 1 : sizes(1) % rows
        
        % data to add
        mask_line = fread(fileMASK, hdrmask.samples, precisionMask, 0, machineformatMask)';
        
        [line,~] = read_optic_data_Line(image,band_name,r);
        line = permute(line,[2,1]);
        % masking
        data_mask = line(:,mask_line == 1);
        
        % mean subtraction
        for b = 1 : sizes(3)
            data_mask(b,:) = data_mask(b,:) - mean_vec(b);
        end
        data_mask = data_mask.*data_mask;
        amount = amount + sum(data_mask,2);
    end
    std_vec = sqrt(amount/((sizes(2)*sizes(1))-1));
    fclose(fileMASK);
end
end
