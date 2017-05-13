function I=rgb2double(I)
% RGB2DOUBLE  Convert RGB image to double format
%   I=RGB2DOUBLE(I) converts the RGB image J from class 'uint8',
%   'uint16' or 'double' to class 'double.
%
%   See RGB2XYZ().

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

switch class(I)
  case 'uint8'
    I = double(I)/(2^8-1) ;
    
  case 'uint16'
    I = double(I)/(2^16-1) ;
    
  case 'double'
  otherwise
    error(['Data type ''', class(I), ''' not supported.']) ;
end
