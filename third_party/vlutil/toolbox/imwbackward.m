function wI = imwbackward(varargin)
% IMWBACKWARD  Image backward warping
%   J = IMWBACKWARD(I, X, Y) returns the values of image I at
%   locations X,Y. X and Y are real matrices of arbitrary but
%   identical dimensions. I is bilinearly interpolated between
%   samples and extended with NaNs to the whole real plane.
%
%   By default, IMWBACKWARD assumes that the image I uses the standard
%   coordinate system (see WARP_OVERVIEW). IMWBACKWARD(XR,YR,I,X,Y)
%   assumes instead that I is defined on a rectangular grid
%   specified by the vectors XR and YR.
%
%   ALGORITHM. What IMWBACKWARD is less general than the MATLAB native
%   function INTERP2, but significantly faster.
%
%   See also WARP_OVERVIEW, IMWFORWARD, INTERP2.

if nargin < 5
  I = varargin{1} ;
  [M,N,K] = size(I) ;
  xr = 1:N ;
  yr = 1:M ;
  varargin = { varargin{2:end} } ;
else
  xr = varargin{1} ;
  yr = varargin{2} ;
  I  = varargin{3} ;
  [M,N,K] = size(I) ;
  varargin = { varargin{4:end} } ;
end
  
if K == 1
  wI = imwbackwardmx(xr, yr, I, varargin{:}) ;
else
  [M,N] = size(varargin{1}) ;
  wI = zeros(M,N,K) ;
  for k=1:K
    wI(:,:,k) = imwbackwardmx(xr, yr, squeeze(I(:,:,k)),  varargin{:}) ;
  end
end
