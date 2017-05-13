% DRAWEDGELIST - plots pixels in edgelists
%
% Usage:      drawedgelist(edgelist, rowscols, randcol, figno)
%
%    edgelist   - Cell array of edgelists in the form
%                     { [r1 c1   [r1 c1   etc }
%                        ...
%                        rN cN]   ....]
%    rowscols -   2 element vector [rows cols] specifying the size of the
%                 image from which edges were detected (used to set size
%                 of plotted image) 
%    randcol    - Optional flag (1/0) to turn on/off random color coding for
%                 each edgelist so that it is easier to see how the edges
%                 have been broken up into separate lists.  Default is 0.
%    figno      - Optional figure number in which to display image.
%
% See also: EDGELINK, LINESEG

% Copyright (c) 2003-2004 Peter Kovesi
% School of Computer Science & Software Engineering
% The University of Western Australia
% http://www.csse.uwa.edu.au/
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.

% February  2003 - Original version
% September 2004 - Revised and updated

function drawedgelist(edgelist, rowscols, randcol, figno)
    
    if nargin == 4
	figure(figno);
    end
    
    if nargin == 2
	randcol = 0;
    end
    
    Nedge = length(edgelist);

    if randcol
	colourmp = hsv(Nedge);    % HSV colour map with Nedge entries
	colourmp = colourmp(randperm(Nedge),:);  % Random permutation
	for I = 1:Nedge
	    line(edgelist{I}(:,2), edgelist{I}(:,1),'color',colourmp(I,:));
	end	
    else
	for I = 1:Nedge
	    line(edgelist{I}(:,2), edgelist{I}(:,1));
	end	
    end

    % Check whether we need to expand bounds
    
    minx = 1; miny = 1;
    maxx = rowscols(2); maxy = rowscols(1);
    
    for I = 1:Nedge
	minx = min(min(edgelist{I}(:,2)),minx);
	miny = min(min(edgelist{I}(:,1)),miny);
	maxx = max(max(edgelist{I}(:,2)),maxx);
	maxy = max(max(edgelist{I}(:,1)),maxy);	
    end	    
    
    %    axis([1 rowscols(2) 1 rowscols(1)]);
    axis([minx maxx miny maxy]);
    axis('equal'); axis('ij');
