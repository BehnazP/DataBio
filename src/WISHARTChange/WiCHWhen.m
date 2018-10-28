function when = WiCHWhen(wc,pval)

% WiCHWhen takes struct wc from wcRjl as input and outputs occurrence of change
%
% when = WiCHWhen(wc,0.9999);
%
% input:
% wc 	-	is struct from WiCH, 
% pval 	-   is probability threshold for change prob, 
%			(default 0.99)
% output:
% when 	-	is 3D matrix (nrow x ncol x ntimes) showing if lnR shows change after each time point for each pixel

% (when) Copyright 2016
% Allan Aasbjerg Nielsen, PhD
% alan@dtu.dk, http://people.compute.dtu.dk/alan
% 2 May 2016

% modified by Behnaz Pirzamanbein
% bepi@dtu.dk
% Image Analysis and Computer Graphics section
% Department of Applied Mathematics and Computer Science
% Technical University of Denmark
% First version 24 May 2018


if ~isstruct(wc)
    error('wc must be struct as output from wcRjl')
end

if nargin<2
    pval = 0.99;
end

% if size(wc(1).P,3) == 1 % two time points only: copy lnQ to lnR
%     wc(1).P(:,:,2) = wc(1).P;
%     wc(1).P(:,:,1) = wc(1).P;
% end

[nrow,ncol,ntime] = size(wc(1).P);

when = zeros(nrow,ncol,ntime-1);

if nrow ~= 1
    h = waitbar(0,'Processing...','Name','WiCHWhen');
end
for col = 1:ncol
    if nrow ~= 1
        waitbar(col/ncol,h)
    end
    for row = 1:nrow
        time = 0;
        count = 1;
        while and( count < ntime , time < ntime )
            for j = 2:( ntime - count + 1 ) % no check of lnQ
                time = time + 1;
                if wc(count).P(row,col,j) > pval
                    when(row,col,time) = true;
                    count = count + j - 1;
                    break;
                end
            end
        end
    end
end
if nrow ~=1
    close(h)
end

end
