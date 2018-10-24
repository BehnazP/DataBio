function [mean_vec] = meanEval(image,band_name,varargin)
% Function that performs the mean of an updating data set.
% It is based on the provisional mean algorithm. The algorithm reads 
% each images line by line from a stored image.
% --------------------------------
% input
%   image
%   band_name
%   mask
%   weights
%
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
% 20/06/2014
%
% ---------------------------------  

if nargin < 2
    disp('COVMat must have at least 1 input: image!/n')
    return
    
elseif nargin == 2
    mask = 0;    
    %[hdr, precision, machineformat] = envihdrread([input,'.hdr']);
    [~,~,sizes] = read_optic_data_Line(image,band_name);
    weight = ones(sizes(1),sizes(2));
    
elseif nargin > 2
    mask = varargin{1};
    if ~isequal(mask,0)
        [hdrmask, precisionMask, machineformatMask]=envihdrread([mask,'.hdr']);
    end
    %[hdr, precision, machineformat] = envihdrread([input,'.hdr']);
    [~,~,sizes] = read_optic_data_Line(image,band_name);
    weight = ones(sizes(1),sizes(2));

elseif nargin > 3
    mask = varargin{1};
    if ~isequal(mask,0)
        [hdrmask, precisionMask, machineformatMask]=envihdrread([mask,'.hdr']);
    end
    %[hdr, precision, machineformat] = envihdrread([input,'.hdr']);
    [~,~,sizes] = read_optic_data_Line(image,band_name);
    weight = varargin{2};
    
end
NN = sizes(3);
sw = 0;
mn = zeros(NN,1);

if isequal(mask,0) % without mask
        
    for r = 1 : sizes(1) % rows
        
        % data to add
        [line,~] = read_optic_data_Line(image,band_name,r);
        line = permute(line,[2,1]);
        Ws = weight(r,:);
        
        % loop over all observation
        for i = 1 : size(line,2)
            sw = sw + Ws(i);
            c = Ws(i)/sw;
            % mean
            for j = 1 : NN
                d(j) = (line((i-1)*NN+j)) - mn(j) ;
                mn(j) = mn(j) + d(j)*c;
            end
        end
    end
    
elseif ~isequal(mask,0) % with mask
    
    fileMASK = fopen(mask, 'r');
    for r = 1 : sizes(1)% rows
        
        % data to add
        [line,~] = read_optic_data_Line(image,band_name,r);
        line = permute(line,[2,1]);
        mask_line = fread(fileMASK, hdrmask.samples, precisionMask, 0, machineformatMask)';
        
        % mask application
        Ws = weight(r,mask_line == 1);
        data_mask = line(:,mask_line == 1);
        
        % loop over all observation
        for i = 1 : size(data_mask,2)
            sw = sw + Ws(i);
            c = Ws(i)/sw;
            % mean
            for j = 1 : NN
                d(j) = (data_mask((i-1)*NN+j)) - mn(j) ;
                mn(j) = mn(j) + d(j)*c;
            end
        end
    end
    fclose(fileMASK);
end
mean_vec = mn;
end