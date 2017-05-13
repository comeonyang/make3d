function o = label_graph(ifile, xlabs, xminval, xmaxval, ylabs, yminval, ymaxval)
% Draw Axis for graph, on left and bottom of bounding rect.
% Add hash marks for labels.
% Also, write stuff to close file.

o = 1;
irec_lx = 50;
irec_ly = 50;
irec_ux = 500;
irec_uy = 500;
lab_length = 20;
% We write to an image rectangle with l giving lower corner, u upper.
% This must be same values as in graph_idraw.

fid = fopen(ifile, 'a');
line_idraw(fid, irec_lx, irec_ly, irec_lx, irec_uy);
line_idraw(fid, irec_lx, irec_ly, irec_ux, irec_ly);

x_len = irec_ux - irec_lx;
xvals = irec_lx + round(x_len*(xlabs - xminval)./(xmaxval-xminval));

y_len = irec_uy - irec_ly;
yvals = irec_ly + round(y_len*(ylabs - yminval)./(ymaxval-yminval));

for x = xvals
  line_idraw(fid, x, irec_ly, x, irec_ly - lab_length);
end

for y = yvals
  line_idraw(fid, irec_lx - lab_length, y, irec_lx, y);
end

fprintf(fid, '\n\nEnd %%I eop\n\nshowpage\n\n%%%%Trailer\n\nend');
fclose(fid);

