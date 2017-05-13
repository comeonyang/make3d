function o = line_idraw(fid, x1, y1, x2, y2)
% Write to ifile the .ps code to draw a line.  Assume file is open.

fprintf(fid, '\nBegin %%I Line\n%%I b 65535\n0 0 0 [] 0 SetB\n%%I cfg Black\n');
fprintf(fid, '0 0 0 SetCFg\n%%I cbg White\n1 1 1 SetCBg\nnone SetP %%I p n');
fprintf(fid, '%%I t\n[ 1 -0 -0 1 173 165 ] concat\n%%I\n');
fprintf(fid, '%d %d %d %d Line\n', x1, y1, x2, y2);
fprintf(fid, '%%I 1\nEnd\n');
