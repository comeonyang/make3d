function H = histmarg(H, dims)
% HISTMARG  Marginal of histogram
%  H = HISTMARG(H, DIMS) marginalizes the historgram HISTMARG w.r.t
%  the dimensions DIM. This is done by summing out all dimensions
%  not listed in DIM and deleting them.
%
%  REMARK. If DIMS lists only one dimension, the returned H is a
%  column vector. This way of deleting dimensions is not always
%  consistent with the SQUEEZE function.

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

sz = size(H) ;

for d=setdiff( 1:length(sz), dims(:) )
  H = sum(H, d) ;
end

% Squeeze out marginalized dimensions
sz = sz(dims(:)) ;
sz = [sz ones(1,2-length(dims(:)))] ;
H = reshape(H, sz) ;
