function [det negdet] = wishart_det(X,pol)
%
% WISHART_DET calculates determinants of complex covariance polarimetric SAR data
%
% [det negdet] = wishart_det(X,pol);
%
% Input
% X       - covariance matrix input at first time point;
%	pol	    - 9 = 'full',
%           5 = 'azim',
%           4 = 'dual',
%           3 = 'fdiag',
%           2 = 'ddiag',
%           1 = 'single'%
%
% Output
% det     - determinant of X
% negdet  - identifies row and col number of pixels with negative determinants
%
% Note:     hhhh (real), hhhv (complex), hhvv (complex),
%                        hvhv (real),    hvvv (complex),
%                                        vvvv (real)
% 'full' full
% ShhShh* ShhShv* ShhSvv*    k     a    rho   | 1   2,3   4,5
% ShvShh* ShvShv* ShvSvv* =  a*    ksi  b     |     6     7,8
% SvvShh* SvvShv* SvvSvv*    rho*  b*   zeta  |           9
% det = k x ksi x zeta + a x b x rho* + rho x a* x b* - rho^2 x ksi - b^2 x k - zeta x a^2
% a x b x rho* = a* x b* x rho = 2 x Re(a x b x rho*)
%
% 'azim' azimuthal symmetry
% ShhShh* ShhShv* ShhSvv*    k     0    rho  |  1      2,3
% ShvShh* ShvShv* ShvSvv* =  0     ksi  0    |     4
% SvvShh* SvvShv* SvvSvv*    rho*  0    zeta |         5
% det = ksi x (k x zeta - rho^2)
%
% 'fdiag' full diagonal
% ShhShh* ShhShv* ShhSvv*    k     0    0    |  1
% ShvShh* ShvShv* ShvSvv* =  0     ksi  0    |     2
% SvvShh* SvvShv* SvvSvv*    0     0    zeta |        3
% det = k x ksi x zeta
%
% 'dual' dual polarimetry
% ShhShh* ShhSvv*    k     rho   | 1   2,3
% SvvShh* SvvSvv*    rho*  ksi   |     4
% or
% ShhShh* ShhShv*    k     a     | 1   2,3
% ShvShh* ShvShv*    a*    zeta  |     4
% or
% ShvShv* ShvSvv*    ksi   b     | 1   2,3
% SvvShv* SvvSvv*    b*    zeta  |     4
% det  = k x ksi - rho^2;
%
% 'ddiag' dual diagonal
% ShhShh* ShhSvv*    k     0     |  1
% SvvShh* SvvSvv*    0     zeta  |     2
% det  = k x zeta
%
% (c) Copyright 2012
% Allan Aasbjerg Nielsen
% aa@space.dtu.dk, www.imm.dtu.dk/~aa
% 23 Jan 2012
%
% modified by Behnaz Pirzamanbein
% bepi@dtu.dk
% Image Analysis and Computer Graphics section
% Department of Applied Mathematics and Computer Science
% Technical University of Denmark
% First version 27 May 2018
% last version 2018-11-03

if  ndims(X)<2
    error('wishart_det: wrong dimensionality of input');
end

switch pol
	case 9 % det = k*ksi*zeta+a*b*rho*+rho*a**b*-rho^2*ksi-b^2*k-zeta*a^2
           % a*b*rho* = a**b**rho = 2*Re(a*b*rho*)
		det = X(:,:,1).*X(:,:,4).*X(:,:,6) ...
			+ 2*(real(X(:,:,3)).*(real(X(:,:,2)).*real(X(:,:,5))- imag(X(:,:,2)).*imag(X(:,:,5)))+ ...
              imag(X(:,:,3)).*(imag(X(:,:,2)).*real(X(:,:,5))+real(X(:,:,2)).*imag(X(:,:,5)))) ...
			- X(:,:,4).*(real(X(:,:,3)).^2+imag(X(:,:,3)).^2) ...
			- X(:,:,1).*(real(X(:,:,5)).^2+imag(X(:,:,5)).^2) ...
			- X(:,:,6).*(real(X(:,:,2)).^2+imag(X(:,:,2)).^2);
	case 5
		det = X(:,:,3).*(X(:,:,1).*X(:,:,4) - real(X(:,:,2)).^2 - imag(X(:,:,2)).^2);
	case 4
		det = X(:,:,1).*X(:,:,3) - real(X(:,:,2)).^2 - imag(X(:,:,2)).^2;
	case 3
		det = X(:,:,1).*X(:,:,2).*X(:,:,3);
	case 2
		det = X(:,:,1).*X(:,:,2);
	case 1
		det = X(:,:,1);
	otherwise
	    error('wishart_det: input (9, 5, 3), (4, 2) or 1 band(s), real (and complex)');
end

negdet.r = [];
negdet.c = [];

if min(det(:))<=0
    warning('*** det <= 0, bad ***')
   [negdet.r, negdet.c] = find(det<=0);
end
end
