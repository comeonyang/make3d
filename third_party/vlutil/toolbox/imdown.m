function J = imdown(I,method)
% IMDOWN  Downsample image 
%   J = IMDOWN(I,'sample') downsamples the image I by half by
%   discarding each other pixel.
%   
%   IMDOWN(I,'avg') downsmples by averaging groups of four pixels.
%
%   The image is always converted to double format.
%   
%   See also IMUP().

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

I = rgb2double(I) ;

if nargin < 2
  method = 'sample' ;
end

switch method
  case 'sample'
    J = I(1:2:end,1:2:end,:) ;
    
  case 'avg'
    J = ...
        I(1:2:end-1,1:2:end-1,:) + ...
        I(2:2:end,1:2:end-1,:) + ...
        I(1:2:end-1,2:2:end,:) + ...
        I(2:2:end,2:2:end,:) ;
    J = J / 4 ;
end
