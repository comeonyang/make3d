function o = save_impute(filename, e, a, s, call, fo)

fid = fopen(filename,'a');

fprintf(fid, '%% %s\n\n', call);

fprintf(fid, 'e%i = ...\n', fo);
fprintf(fid, '[');
fprintf(fid, '%4.5f,   %4.5f,   %4.5f; ...\n',e');
fprintf(fid, ']\n\n');

fprintf(fid, 'a%i = ...\n', fo);
fprintf(fid, '[');
fprintf(fid, '%4.5f,   %4.5f,   %4.5f; ...\n',a');
fprintf(fid, ']\n\n');

fprintf(fid, 's%i = ...\n', fo);
fprintf(fid, '[');
fprintf(fid, '%4.5f; ...\n',s);
fprintf(fid, ']\n\n');

fclose(fid);

o = 1;
