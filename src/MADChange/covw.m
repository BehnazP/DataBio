function [dispersion,meanw] = covw(x,w,varargin)

% [covaw meanw] = covw(x,w,varargin)
%
%   COVW Covariance matrix with weights.
%
%   COVW(W,X), if X is a vector, returns the variance.  For matrices,
%   where each row is an observation, and each column a variable,
%   COVW(X,W) is the variance-covariance matrix.  DIAG(COVW(X,W)) is
%   a vector of variances for each column, and SQRT(DIAG(COVW(X,W)))
%   is a vector of standard deviations. 
%
%   Observations are weighted with vector of weights, W, which has all
%   nonnegative values.  COVW(X,W) gives the weighted estimate of the
%   variance-covariance matrix.
%   
%   COVW(X,W) normalizes by (N-1) where N is the number of
%   observations.
%
%   COVW(X,W,1) normalizes by N and produces the second
%   moment matrix of the observations about their mean.
%   COVW(X,W,0) is the same as COVW(X,W).
%
%   The weighted mean is removed from each column before calculating the
%   result.

%   J. Little 5-5-86
%   Revised 6-9-88 LS 3-10-94 BJ
%   Copyright 1984-2002 The MathWorks, Inc. 
%   $Revision: 5.16 $  $Date: 2002/06/05 17:06:38 $
%
%   Modified for weights
%   Allan Aasbjerg Nielsen
%   aa@imm.dtu.dk
%   10 Feb 2005
%
%   Modified
%   Behnaz Pirzamanein
%   bepi@dtu.dk
%   11 Jan 2018

if nargin<2
    error('Not enough input arguments.'); 
end
if nargin>3
    error('Too many input arguments.'); 
end

w = w(:); % w is now a column vector
[m,n] = size(x);
mw = size(w,1);
if length(x)==(m*n)
    x = x(:); % x is now a column vector
    m = size(x,1);
end
if m~=mw
    error('x and w must match.')
end 
if ~isempty(w(w<0)) 
    error('Weights must be nonnegative.'); 
end

flag = 0;
% Check for covw(x,w) or covw(x,w,flag)
if (nargin==3)
    flag = varargin{end};
end

if m == 1  % Handle special case
    dispersion = 0;
else
    sumw = sum(w);
    meanw = sum(bsxfun(@times,w,x)) / sumw;
    xc = bsxfun(@minus, x, meanw);  % Remove weighted mean
    xc = bsxfun(@times,xc,sqrt(w));
    dispersion = xc' * xc / sumw;
    if flag==0
        dispersion = m/(m-1) * dispersion;
    end
end