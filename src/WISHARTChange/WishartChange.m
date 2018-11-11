function [P,ff,rho,omega2,lnR,sumX,logDet] = WishartChange(X,noL,pol)
%
% Calculate ln(Q) and ln(Rj) for change detection in polarimetric SAR data
%
% detect difference/change in k 3x3, 2x2 or 1x1 complex Wishart distributed
% matrices used for change detection in polarimetric SAR images
%
% [P,ff,rho,omega2,lnR,sumX,logDet] = WishartChange(X,noL,pol);
%
% Input
% images     - a 4D matrix containing images (NxM) for P different polarization
%                and T different time points; images (NxMxPxT)
% noL        - number of looks for images
% pol	     - 9 = 'full',
%              5 = 'azim',
%              4 = 'dual',
%              3 = 'fdiag',
%              2 = 'ddiag',
%              1 = 'single'
% Output
% lnR     - logarithm of test quantities Q and Rj,
%           For example: ln(Q)   is lnR(:,:,1),
%                        ln(R_2) is lnR(:,:,2) and etc.
%           Note: lnQ is the sum of the lnRs, i.e., lnQ = sum(lnR(:,:,2:end),3)
% ff      - number of degrees of freedom for Q and Rj test statistics
% rho     - constant terms for calculating the probabilities; rho(1)
%           correspond to ln(Q) and the rest correspond to lnRj
% omega2  - constant terms for calculating the probabilities; omega(1)
%           correspond to ln(Q) and the rest correspond to lnRj
% P       - probabilities of finding smaller value of
%             -2*rho(i)*lnR(:,:,i), i=1,...,k
%           or for npol=1,2 & 3 probabilities of finding smaller value of
%             -2*rho(i)*lnR(:,:,i), i=1,...,k
% sumX    - sum of the noL*X (for efficient updating) for the last time
% logDet  - logarithm of the determinant of sumX for the last time
%
% Reference
% Knut Conradsen, Allan Aasbjerg Nielsen and Henning Skriver (2016):
% "Determining the Points of Change in Time Series of Polarimetric SAR Data".
% IEEE Transactions on Geoscience and Remote Sensing 54(5), 3007-3024.
% https://doi.org/10.1109/TGRS.2015.2510160
%
% (c) Copyright 2014-2016
% Allan Aasbjerg Nielsen, PhD
% alan@dtu.dk, http://people.compute.dtu.dk/alan
% 6 Dec 2016
%
% modified by Behnaz Pirzamanbein
% bepi@dtu.dk
% Image Analysis and Computer Graphics section
% Department of Applied Mathematics and Computer Science
% Technical University of Denmark
% First version 11 April 2018
% last version 2018-11-03

if nargin>3
    error('WishartChange: wrong input');
end
if nargout>7
    error('WishartChange: wrong output');
end

[nrows,ncols,npol,ntime] = size(X);

if ntime < 2
    error('Data should be at least for 2 times')
end

X = noL.*X;

switch pol
    case 1
        p = 1;
        f = p^2;
    case 2
        p = 1+1;
        f = 1^2+1^2;
    case 3
        p = 1+1+1;
        f = 1^2+1^2+1^2;
    case 4
        p = 2;
        f = p^2;
    case 5
        p = 2 + 1;
        f = 2^2 + 1^2;
    case 9
        p = 3;
        f = p^2;
    otherwise
        error('Wrong number of polarization')
end

%%%%% Bezo, what about the rest of pol?
fQ = (ntime-1)*f; % we think this is true for 'dd', 'fd', and 'a' also

% block diagonal case: P-values for -2*ln(R) - and not for -2*rho*ln(R)
if pol == 2 || pol == 3 || pol == 1
    rho = ones(ntime,1);
    omega2 = zeros(ntime,1);
else
    rho = nan(ntime,1);
    omega2 = nan(ntime,1);

    rho(1) = 1 - (2*p^2-1)*(ntime/noL-1/noL/ntime)/(6*(ntime-1)*p);
    omega2(1) = -(ntime-1)*p^2/4*(1-1/rho(1))^2 + ...
                p^2*(p^2-1)*(ntime/noL^2-1/noL^2/ntime^2)/(24*rho(1)^2);

    for j=2:ntime
        rho(j) = 1 - (2*f-1)/(6*p*noL)*(1+1/(j*(j-1)));
        omega2(j) = -f/4*(1-1/rho(j))^2 ...
                    + f*(f-1)*(1+(2*j-1)/(j^2*(j-1)^2))/(24*noL^2*rho(j)^2);
    end
end

dfQ = fQ/2; % this is 0.5*f; for fewer multiplications in call of gammainc
dfQp4 = dfQ + 2; % this is 0.5*(f+4)
df = f/2; % this is 0.5*f; for fewer multiplications in call of gammainc
dfp4 = df + 2; % this is 0.5*(f+4)

ff = f*ones(ntime,1);
ff(1) = fQ;

sumX = nan(nrows,ncols,npol);
lnR = nan(nrows,ncols,ntime);

for j=2:ntime
    cst = p*(j*log(j)-(j-1)*log(j-1));
    sumX = sum(X(:,:,:,1:(j-1)),4);
    logDet = log(wishart_det(sumX,pol));
    lnR(:,:,j) = cst + (j-1)*logDet;
    logDet = log(wishart_det(X(:,:,:,j),pol));
    lnR(:,:,j) = lnR(:,:,j) + logDet;
    sumX = sumX + X(:,:,:,j);
    logDet = log(wishart_det(sumX,pol));
    lnR(:,:,j) = noL*(lnR(:,:,j) - j*logDet);
end

if nargout > 5
    sumX = sum(X,4);
    logDet = log(wishart_det(sumX,pol));
end
clear X

if ntime == 2
    lnR(:,:,1) = lnR(:,:,2); % index "2" is not a mistake!
    omega2(2) = omega2(1);
    rho(2) = rho(1);
    ff(2) = ff(1);
    P = (1 - omega2(1))*gammainc(-rho(1)*lnR(:,:,2),dfQ) + ...
             omega2(1) *gammainc(-rho(1)*lnR(:,:,2),dfQp4);
    %%%%%%%%%%%%%%%%%%BEZO added to read the WiCHProb%%%%%%%%%%%%%%%%%%%%%
    P(:,:,1) = P;
    P(:,:,2) = P;
else
    lnR(:,:,1) = sum(lnR(:,:,2:end),3); % this is ln(Q)
    P = nan(nrows,ncols,ntime);
    for j=[2:ntime 1]
        if j==1, df = dfQ; dfp4 = dfQp4; end
        P(:,:,j) = (1 - omega2(j))*gammainc(-rho(j)*lnR(:,:,j),df) + ...
                        omega2(j) *gammainc(-rho(j)*lnR(:,:,j),dfp4);
    end
end

if ~isreal(P) % this is because som lnR/lnQ values are >0 %%%% BEZO numerical instability?
    warning('*** some P-values are complex, only real part is output, max(abs(imag(P(:)))): ***')
    P = real(P);
end
end
