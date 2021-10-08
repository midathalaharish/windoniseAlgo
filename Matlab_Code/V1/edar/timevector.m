function T = timevector(X,fs,t0)
%TIMEVECTOR   Uniformly sampled time vector
%   T = TIMEVECTOR(X,FS) returns a uniformly sampled time vector matching
%   the samples in vector X, starting at T = 0. The sample rate is FS Hz.
%   T = TIMEVECTOR(N,FS) returns N uniformly samples.
%   T = TIMEVECTOR(...,T0) returns a time vector starting at T = T0 [s].

% Copyright 2016-2020 Sigma Connectivity AB
% Author: Claes Hovmalm
% $Id: timevector.m 19 2020-07-06 12:55:15Z hovm1011 $
% $HeadURL: file:///C:/Users/hovm1011/OneDrive%20-%20Sigma%20AB/Tools/MATLAB_repository/trunk/chutilities/timevector.m $

narginchk(2,3)

if ~isscalar(X)
   if ndims(X)~=2
       error('Matrix must be two-dimensional.')
   end
   [n,iDim] = max(size(X)) ;
else
    % Return error if X is not a positive integer
    if ~(isreal(X) && isfinite(X) && (0<X) && (round(X)==X))
        error('Scalar input arguments must be positive integers')
    end
    n = X ;
    iDim = 1 ;
end

if nargin<3
    t0 = 0 ;
end

T = (0:(n-1))'/fs + t0 ;

if iDim==2 % T shall be a row vector
    T = T' ;
end
