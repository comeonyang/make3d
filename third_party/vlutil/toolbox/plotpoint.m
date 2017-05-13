function h=plotpoint(V,varargin)
% PLOTPOINT  Plot 2 or 3 dimensional points
%   PLOTPOINT(V) plots the 2 or 3 dimensional points V. V is a 2xK or
%   3xK array, with one point per column.
%
%   H=PLOTPOINT(...) returns the handle H of the plot.
%
%   PLOTPOINT() is a simple wrapper around the PLOT() and PLOT3()
%   functions. By default, PLOTPOINT(V) plots the points with line
%   style '.'.  PLOTPOINT(V,...) does not use the default line style;
%   rather it passess any extra argument to the underlying plot
%   function.
%
%   See also PLOT(), PLOT3().

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

if length(varargin) == 0 
  varargin = {'.'};
end

switch size(V,1)
  case 2
    h=plot(V(1,:),V(2,:),varargin{:}) ;
  case 3 
    h=plot3(V(1,:),V(2,:),V(3,:),varargin{:}) ;
  otherwise
    error(['V must be either 2xK or 3xK.']) ;          
end
