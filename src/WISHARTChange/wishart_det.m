function [det negdet] = wishart_det(X,pol)
% [Xdet Xnegdet] = wishart_det(X);
%
% Calculate determinants of complex covariance polarimetric SAR data
%
% the order of the pol: HHHH, HHHV, HHVV, HVHV,HVVV, VVVV
%
% full
% ShhShh* ShhShv* ShhSvv*    k     a    rho   | 1   2,3     4,5
% ShvShh* ShvShv* ShvSvv* =  a*    ksi  b     |     6       7,8
% SvvShh* SvvShv* SvvSvv*    rho*  b*   zeta  |             9
% det = k*ksi*zeta+a*b*rho*+rho*a**b*-rho^2*ksi-b^2*k-zeta*a^2
% a*b*rho* = a**b**rho = 2*Re(a*b*rho*)

% azimuthal symmetry
% ShhShh* ShhShv* ShhSvv*    k     0    rho  |  1           2,3
% ShvShh* ShvShv* ShvSvv* =  0     ksi  0    |     4
% SvvShh* SvvShv* SvvSvv*    rho*  0    zeta |              5
% det = ksi*(k*zeta-rho^2)

% diagonal only
% ShhShh* ShhShv* ShhSvv*    k     0    0    |  1
% ShvShh* ShvShv* ShvSvv* =  0     ksi  0    |     2
% SvvShh* SvvShv* SvvSvv*    0     0    zeta |        3
% det = k*ksi*zeta

% dual polarimetry
% ShhShh* ShhSvv*    k     rho   | 1   2,3
% SvvShh* SvvSvv*    rho*  ksi   |     4
% or
% ShhShh* ShhShv*    k     a     | 1   2,3
% ShvShh* ShvShv*    a*    zeta  |     4
% or
% ShvShv* ShvSvv*    ksi   b     | 1   2,3
% SvvShv* SvvSvv*    b*    zeta  |     4
% det  = k*ksi-rho2;

% diagonal only
% ShhShh* ShhSvv*    k     0     |  1
% SvvShh* SvvSvv*    0     zeta  |     2
% det  = k*zeta

% Input
%   X       - covariance matrix input at first time point;
%             must be 3-D with either
%               9 (full polarimetry),
%               5 (full, azimuthal symmetry),
%               3 (full, diagonal only),
%               4 (dual polarimetry),
%               2 (dual, diagonal only), or
%               1 band(s);
%             reading order is alphabetical in indexes with
%             real parts before imaginary parts (as read by
%             freadenvisar), fx. for full polarimetry:
%               hhhh (real), hhhv (complex), hhvv (complex),
%               hvhv (real), hvvv (complex), vvvv (real)
%
% Output
%   Xdet    - determinant of X
%   Xnegdet - identifies row and col number of pixels with negative
%             determinants

% (c) Copyright 2012
% Allan Aasbjerg Nielsen
% aa@space.dtu.dk, www.imm.dtu.dk/~aa
% 23 Jan 2012

% modified by Behnaz Pirzamanbein
% bepi@dtu.dk
% Image Analysis and Computer Graphics section
% Department of Applied Mathematics and Computer Science
% Technical University of Denmark
% First version 27 May 2018

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
