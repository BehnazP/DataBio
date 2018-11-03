function wc = WiCHParallel(X,noL,pol)
%
% WICHPARALLEL is a Wishart based change detection in time series of polarimetric
%              SAR images that computes the function WishartChange in parallel scheme
%
% wc = WiCHParallel(X,noL,npols)
%
% Input
% images     -  a 4D matrix containing images (NxM) for P different polarization
%                and T different time points; images (NxMxPxT)
% noL        - number of looks for X
%	pol	       - 9 = 'full',
%              5 = 'azim',
%              4 = 'dual',
%              3 = 'fdiag',
%              2 = 'ddiag',
%              1 = 'single'
%
% Output
% wc         - structure with
%             f       - number of degrees of freedom for Q and Rj test statistics
%             rho     - constant terms for calculating the probabilities; rho(1)
%                       correspond to ln(Q) and the rest correspond to lnRj
%             omega2  - constant terms for calculating the probabilities; omega(1)
%                       correspond to ln(Q) and the rest correspond to lnRj
%             P       - probabilities of finding smaller value of -2*rho(i)*lnR(:,:,i), i=1,...,k
%
% Reference
% Knut Conradsen, Allan Aasbjerg Nielsen and Henning Skriver (2016):
% "Determining the Points of Change in Time Series of Polarimetric SAR Data".
% IEEE Transactions on Geoscience and Remote Sensing 54(5), 3007-3024.
% http://www.imm.dtu.dk/pubdb/p.php?6825
% https://bit.ly/2J7otCk
%
% (c) Copyright 2014-2017
% Allan Aasbjerg Nielsen, PhD
% alan@dtu.dk, http://people.compute.dtu.dk/alan
% 29 Oct 2014 - 11 Aug 2017
%
% modified by Behnaz Pirzamanbein
% bepi@dtu.dk
% Image Analysis and Computer Graphics section
% Department of Applied Mathematics and Computer Science
% Technical University of Denmark
% First version 24 May 2018
% last version 2018-11-03

if nargin<2 || nargin>4
    help WiCH
    error('WiCH: wrong input');
end
if nargout>1
    help WiCH
    error('WiCH: wrong output');
end

[~,~,~,ntime] = size(X);
iter = ntime-1;
h = parfor_progressbar(iter,'Processing...','Name','WiCH');

wc = struct('f',nan,'rho',nan,'omega2',nan,'P',nan);
parfor i = 1:iter
    h.iterate(i/(ntime-1))
    [wciP,wcif,wcirho,wciomega2] = WishartChange(X(:,:,:,i:ntime),noL,pol);
    wc(i).f = single(wcif);
    wc(i).rho = single(wcirho);
    wc(i).omega2 = single(wciomega2);
    wc(i).P = single(wciP);
    wciP = [];
end
close(h)
end
