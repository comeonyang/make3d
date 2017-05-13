%s_err2 = dlmread('transexps/s_err_.2', ',');
%s_err3 = dlmread('transexps/s_err_.3', ',');
%s_err4 = dlmread('transexps/s_err_.4', ',');
%s_err5 = dlmread('transexps/s_err_.5', ',');
%s_err6 = dlmread('transexps/s_err_.6', ',');
%s_err7 = dlmread('transexps/s_err_.7', ',');

%s_stab2 = ones(500,11);
%s_stab3 = ones(500,11);
%s_stab4 = ones(500,11);
%s_stab5 = dlmread('transexps/s_stab_.5', ',');
%s_stab6 = dlmread('transexps/s_stab_.6', ',');
%s_stab7 = dlmread('transexps/s_stab_.7', ',');

%[R_s2, means_s2, stds_s2] = trans_best(s_err2, s_stab2, 0);
%[R_s3, means_s3, stds_s3] = trans_best(s_err3, s_stab3, 0);
%[R_s4, means_s4, stds_s4] = trans_best(s_err4, s_stab4, 0);
%[R_s5, means_s5, stds_s5] = trans_best(s_err5, s_stab5, 0);
%[R_s6, means_s6, stds_s6] = trans_best(s_err6, s_stab6, 0);
%[R_s7, means_s7, stds_s7] = trans_best(s_err7, s_stab7, 0);

%R_s = [R_s2(2,:)',R_s3(2,:)',R_s4(2,:)',R_s5(2,:)',R_s6(2,:)',R_s7(2,:)'];

%means = [means_s2', means_s3', means_s4', means_s5', means_s6', means_s7'];
%stds = [stds_s2', stds_s3', stds_s4', stds_s5', stds_s6', stds_s7'];

%s_aff2 = dlmread('transexps/s_aff_.2', ',');
%s_aff3 = dlmread('transexps/s_aff_.3', ',');
%s_aff4 = dlmread('transexps/s_aff_.4', ',');
%s_aff5 = dlmread('transexps/s_aff_.5', ',');
%s_aff6 = dlmread('transexps/s_aff_.6', ',');
%s_aff7 = dlmread('transexps/s_aff_.7', ',');

%[R_s2, affmeans_s2, affstds_s2] = trans_best(s_aff2, s_stab2, 0);
%[R_s3, affmeans_s3, affstds_s3] = trans_best(s_aff3, s_stab3, 0);
%[R_s4, affmeans_s4, affstds_s4] = trans_best(s_aff4, s_stab4, 0);
%[R_s5, affmeans_s5, affstds_s5] = trans_best(s_aff5, s_stab5, 0);
%[R_s6, affmeans_s6, affstds_s6] = trans_best(s_aff6, s_stab6, 0);
%[R_s7, affmeans_s7, affstds_s7] = trans_best(s_aff7, s_stab7, 0);

%affmeans = [affmeans_s2', affmeans_s3', affmeans_s4', affmeans_s5', affmeans_s6', affmeans_s7'];
%affstds = [affstds_s2', affstds_s3', affstds_s4', affstds_s5', affstds_s6', affstds_s7'];

r_err2 = dlmread('transexps/r_err_.2', ',');
r_err3 = dlmread('transexps/r_err_.3', ',');
r_err4 = dlmread('transexps/r_err_.4', ',');
r_err5 = dlmread('transexps/r_err_.5', ',');
r_err6 = dlmread('transexps/r_err_.6', ',');
r_err7 = dlmread('transexps/r_err_.7', ',');

r_stab2 = dlmread('transexps/r_stab_.2', ',');
r_stab3 = dlmread('transexps/r_stab_.3', ',');
r_stab4 = dlmread('transexps/r_stab_.4', ',');
r_stab5 = dlmread('transexps/r_stab_.5', ',');
r_stab6 = dlmread('transexps/r_stab_.6', ',');
r_stab7 = dlmread('transexps/r_stab_.7', ',');

[R_r2, means_r2, stds_r2] = trans_best(r_err2, r_stab2, 0);
[R_r3, means_r3, stds_r3] = trans_best(r_err3, r_stab3, 0);
[R_r4, means_r4, stds_r4] = trans_best(r_err4, r_stab4, 0);
[R_r5, means_r5, stds_r5] = trans_best(r_err5, r_stab5, 0);
[R_r6, means_r6, stds_r6] = trans_best(r_err6, r_stab6, 0);
[R_r7, means_r7, stds_r7] = trans_best(r_err7, r_stab7, 0);

R_r = [R_r2(2,:)',R_r3(2,:)',R_r4(2,:)',R_r5(2,:)',R_r6(2,:)',R_r7(2,:)'];

means_r = [means_r2', means_r3', means_r4', means_r5', means_r6', means_r7'];
stds_r = [stds_r2', stds_r3', stds_r4', stds_r5', stds_r6', stds_r7'];

r_aff2 = dlmread('transexps/r_aff_.2', ',');
r_aff3 = dlmread('transexps/r_aff_.3', ',');
r_aff4 = dlmread('transexps/r_aff_.4', ',');
r_aff5 = dlmread('transexps/r_aff_.5', ',');
r_aff6 = dlmread('transexps/r_aff_.6', ',');
r_aff7 = dlmread('transexps/r_aff_.7', ',');

[R_r2, affmeans_r2, affstds_r2] = trans_best(r_aff2, r_stab2, 0);
[R_r3, affmeans_r3, affstds_r3] = trans_best(r_aff3, r_stab3, 0);
[R_r4, affmeans_r4, affstds_r4] = trans_best(r_aff4, r_stab4, 0);
[R_r5, affmeans_r5, affstds_r5] = trans_best(r_aff5, r_stab5, 0);
[R_r6, affmeans_r6, affstds_r6] = trans_best(r_aff6, r_stab6, 0);
[R_r7, affmeans_r7, affstds_r7] = trans_best(r_aff7, r_stab7, 0);

affmeans_r = [affmeans_r2', affmeans_r3', affmeans_r4', affmeans_r5', affmeans_r6', affmeans_r7'];
affstds_r = [affstds_r2', affstds_r3', affstds_r4', affstds_r5', affstds_r6', affstds_r7'];


