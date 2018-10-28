function pP = WiCHProb(wc,ROI)

% print_Ptable takes struct from WiCH as input and outputs table of
% mean (1-P)-values, i.e., no-change probabilities for region of interest, ROI
%
% WiCHProb(wc,ROI,opt);
%
% input
% wc 	-	is struct from WiCH, 
% ROI 	-	is region of interest,
%
% output:
% pP 	-	the mean P-value of the area

% (c) Copyright 2016
% Allan Aasbjerg Nielsen, PhD
% alan@dtu.dk, http://people.compute.dtu.dk/alan
% 14-17 Jan 2016

% modified by Behnaz Pirzamanbein
% bepi@dtu.dk
% Image Analysis and Computer Graphics section
% Department of Applied Mathematics and Computer Science
% Technical University of Denmark
% First version 24 May 2018

if ~isstruct(wc)
    error('wc must be struct as output from WiCH')
end

[~,~,ntime] = size(wc(1).P);

if nargin < 2
    ROI = true(size(wc(end).P(:,:,1)));
end

pP = nan(ntime,ntime-1);

% no-change probability
for count = 1:(ntime - 1)
    Paux = 1 - wc(count).P(:,:,1);
    Paux = Paux(ROI);
    pP(1,count) = mean(Paux);
end
for count = 1:(ntime - 1)
    m = 2;
    for r = (1+count):ntime
        Paux = 1 - wc(count).P(:,:,m);
        Paux = Paux(ROI);
        pP(r,count) = mean(Paux);
        m = m + 1;
    end
end

end