function [es, as, shum_err] = impute_compare(frot,sigma,n)
% This will be a simple comparison between an imputation method of 
% Tomasi and Kanade, and my method.  

% matrix will have form:
% x x x
% x x x
% .
% .
% .
% x x x
% x x x x x x 
% x x x x x x 
% x x x x x x
% x x x x x x
%       x x x
%       x x x
%         .
%         .
%         .
%       x x x 


left_start_rows = 1;
left_end_rows = 18;
full_start_rows = left_end_rows + 1;
full_end_rows = full_start_rows + 3;
right_start_rows = full_end_rows + 1;
right_end_rows = full_end_rows + left_end_rows;

left_start_col = 1;
left_end_col = 4;
right_start_col = left_end_col + 1;
right_end_col = left_end_col + left_end_col;

es = [];
as = [];
shum_err = [];

for i=1:n
[M, pts] = unoccluded_motion(right_end_rows/2, right_end_col, frot, 0);
INC = [[ones(left_end_rows, left_end_col), ...
        zeros(left_end_rows, right_end_col - left_end_col)]; ...
       ones(full_end_rows - left_end_rows, right_end_col); ...
       [zeros(right_end_rows - full_end_rows, left_end_col), ...
        ones(right_end_rows - full_end_rows, right_end_col - left_end_col)]];

ERR = randn(size(M)).*sigma;
Merr = M + ERR;
[e1, Mapprox1, s1] = rankr(Merr, INC, 3, 1000, 0);
a1 = affine_error(Mapprox1, pts, 3);

M_left_top = Merr(left_start_rows:left_end_rows, left_start_col:left_end_col);
M_right_top = Merr(left_start_rows:left_end_rows, right_start_col:right_end_col);
M_left_full = Merr(full_start_rows:full_end_rows, left_start_col:left_end_col);
M_right_full = Merr(full_start_rows:full_end_rows, right_start_col:right_end_col);
M_left_bottom = Merr(right_start_rows:right_end_rows, left_start_col:left_end_col);
M_right_bottom = Merr(right_start_rows:right_end_rows, right_start_col:right_end_col);

IM_right_top = M_left_top*(M_left_full\M_right_full);
IM_left_bottom = M_right_bottom*(M_right_full\M_left_full);

IM = [[M_left_top, IM_right_top];[M_left_full, M_right_full];...
        [IM_left_bottom, M_right_bottom]];
Mapprox2 = approx_full_matrix(IM,3);
e2 = sum(sum((M.*INC - Mapprox2.*INC).^2));
a2 = affine_error(Mapprox2, pts, 3);

[shum_err1, shum_res1] = shum(Merr,INC,3,Mapprox1,100);
[shum_err2, shum_res2] = shum(Merr,INC,3,Mapprox2,100);

es = [es; [e1,e2]];
as = [as; [a1,a2]];
shum_err = [shum_err; [shum_err1,shum_err2]];
end
