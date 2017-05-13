function o = graph_idraw(ifile, means, stds_low, stds_hi, minval, maxval)
% This writes to an idraw file a set of lines that graph a set of function values 
% along with error bars.  
% Values are assumed to lie in between MINVAL and MAXVAL.  
% If STDS has a zero value, no error bar is drawn.  (not implemented)
% stds_low is mean-std, stds_hi is mean+std.  We don't just pass stds 
% in case a log scale is used.

% For this to work, ifile must already contain a template idraw file, generated 
% with a line, that is then removed.  The end of the file must be removed, and must be
% added later.  It looks like:

%End %I eop

%showpage

%%%Trailer

%end

o = 1;
irec_lx = 50;
irec_ly = 50;
irec_ux = 500;
irec_uy = 500;
% We write to an image rectangle with l giving lower corner, u upper.
stdline = 5;
dotrad = 6;

numvals = size(means,2);
x_len = irec_ux - irec_lx;
xspace = x_len/(numvals + 1);
xvals = round(irec_lx + xspace*(1:numvals));

y_len = irec_uy - irec_ly;
yvals = round(irec_ly + y_len*(means - minval)./(maxval-minval));

if ~isempty(stds_low)
  ystdmin = round(irec_ly + y_len*(stds_low - minval)./(maxval-minval));
  ystdmax = round(irec_ly + y_len*(stds_hi - minval)./(maxval-minval));
end

fid = fopen(ifile, 'a');
for i = 1:numvals
   x = xvals(i);
   if ~isempty(stds_low)
     ymin = ystdmin(i);
     ymax = ystdmax(i);
     line_idraw(fid, x, ymin, x, ymax);
     line_idraw(fid, x-stdline, ymin, x+stdline, ymin);
     line_idraw(fid, x-stdline, ymax, x+stdline, ymax);
  end
  dot_idraw(fid, x, yvals(i), dotrad);
end
for i = 1:(numvals - 1)
  x = xvals(i);
  nextx = xvals(i+1);
  y = yvals(i);
  nexty = yvals(i+1);
  line_idraw(fid, x, y, nextx, nexty);
end
fclose(fid);
