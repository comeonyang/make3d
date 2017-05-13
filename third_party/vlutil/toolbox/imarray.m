function [H,J] = imarray(A,varargin)
% IMARRAY  Display array of images
%   IMARRAY(A) plots the array of images A. A can be either a M*N*K
%   array, storing one intensity image per layer, or a M*N*3*K array,
%   storing one RGB image per layer.
%
%   H=IMARRAY(...) returns the handle H of the axis where the image is
%   plotted.
%
%   [H,J]=IMARRAY(...) returns the plotted image J as well.
%
%   The function accepts the following option-value pairs:
%
%   'Spacing' [0]
%     Orlates the images with a black border of the specified width.
%   
%   IMARRAYSC(A, 'Spacing', spacing) is another way of specifying the
%   spacing.
%
%   See also IMARRAYSC().

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

reverse = 1 ;
spacing = 0 ;
lay     = [] ;

for k=1:2:length(varargin)  
  switch varargin{k}
    case 'Layout'
      lay=varargin{k+1} ;
    case 'Spacing'
      sp=varargin{k+1} ;
    otherwise
      error(['Unknown parameter ''', varargin{k}, '''']) ;
  end
end

if ndims(A) <= 3
  d=1 ;
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

cdata = zeros(Mi*M + spacing*(M-1), Ni*N + spacing*(N-1), d) ;
for k=1:K
  p = k - 1 ;
  i = floor(p/N) ;
  if reverse
    i = M-1 - i ;
  end
  j = mod(p,N) ;
  irng = i*(Mi+spacing) + (0:Mi-1) + 1 ;
  jrng = j*(Ni+spacing) + (0:Ni-1) + 1 ;
  if(d == 1)
    J = A(:,:,k) ;
  else
    J = A(:,:,:,k) ;
  end
  cdata(irng,jrng,:) = J ;
end

H = image(cdata) ;
J = cdata ;