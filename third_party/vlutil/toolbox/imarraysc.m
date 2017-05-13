function [H,J] = imarraysc(A, varargin)
% IMARRAYSC  Scale and display an array of images
%   H = IMARRAYSC(A) behaves as IMARRAY(A), but scales the range
%   of each image to fill the interval [0,1].
%
%   In addition to the option-value paris accepted by IMARRAY(),
%   the function accepts also:
%
%   'CLim' [[]] 
%     If not empty, use the specified intensity range.
%
%   See also IMARRAY().

% AUTORIGHTS
% Copyright (C) 2006 Andrea Vedaldi
%       
% This file is part of VLUtil.
% 
% VLUtil is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2, or (at your option)
% any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software Foundation,
% Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

reverse=1 ;
spacing=0 ;
clim=[] ;

% Parse
if(length(varargin)>1 && isnumeric(varargin{1}))
  spacing=varargin{1} ;
  varargin = { varargin{2:end} } ;
end

lay=[] ;

for k=1:2:length(varargin)  
  switch lower(varargin{k})
    case 'layout'
      lay=varargin{k+1} ;
    case 'spacing'
      spacing=varargin{k+1} ;
    case 'clim'
      clim=varargin{k+1} ; 
    otherwise
      error(['Unknown option ''', varargin{k}, '''']) ;
  end
end

if ndims(A) <= 3
  d = 1 ;
  [Mi,Ni,K] = size(A) ;
else
  [Mi,Ni,d,K] = size(A) ;
  if(d ~= 3)
    error(['A should be either M*N*K or M*N*3*K']);
  end
end

if isempty(lay)
  N = ceil(sqrt(K)) ;
  M = ceil(K/N) ;
else
  M = lay(1) ;
  N = lay(2) ;
  K = min(K,M*N) ; 
end

cdata = zeros(Mi*M + spacing*(M-1), Ni*N + spacing*(N-1), d ) ;

for k=1:K
  p = k - 1 ;
  i = floor(p/N) ;
  if reverse
    i = M-1 - i ;
  end
  j = mod(p,N) ;
  irng = i*(Mi+spacing) + (0:Mi-1) + 1 ;
  jrng = j*(Ni+spacing) + (0:Ni-1) + 1 ;
  if d == 1
    J = A(:,:,k) ;
  else
    J = A(:,:,:,k) ;
  end
  if isempty( clim )
    J = J - min(J(:)) ;
    m = max(J(:)) ;
    if(m > eps)
      J = J / m ;
    end
  end
  cdata(irng,jrng,:) = J ;
end

if ~isempty(clim)
  H = imagesc(cdata,clim) ;
else
  H = imagesc(cdata) ;
end

J = cdata ;