
%[err2_10,aff2_10,stab2_10] = comp_transmot(.2,20,20,.25,.5,.005,10);
%[err3_10,aff3_10,stab3_10] = comp_transmot(.3,20,20,.25,.5,.005,10);
%[err4_10,aff4_10,stab4_10] = comp_transmot(.4,20,20,.25,.5,.005,10);
%[err5_10,aff5_10,stab5_10] = comp_transmot(.5,20,20,.25,.5,.005,10);
%[err6_10,aff6_10,stab6_10] = comp_transmot(.6,20,20,.25,.5,.005,10);
%[err7_10,aff7_10,stab7_10] = comp_transmot(.7,20,20,.25,.5,.005,10);

%[err2_100,aff2_100,stab2_100] = comp_transmot(.2,20,20,.25,.5,.005,100);
%[err3_100,aff3_100,stab3_100] = comp_transmot(.3,20,20,.25,.5,.005,100);
%[err4_100,aff4_100,stab4_100] = comp_transmot(.4,20,20,.25,.5,.005,100);
%[err5_100,aff5_100,stab5_100] = comp_transmot(.5,20,20,.25,.5,.005,100);
%[err6_100,aff6_100,stab6_100] = comp_transmot(.6,20,20,.25,.5,.005,100);
%[err7_100,aff7_100,stab7_100] = comp_transmot(.7,20,20,.25,.5,.005,100);

%[err2_500,aff2_500,stab2_500] = comp_transmot(.2,20,20,.25,.5,.005,500);
%[err3_400,aff3_400,stab3_400] = comp_transmot(.3,20,20,.25,.5,.005,400);
%[err4_400,aff4_400,stab4_400] = comp_transmot(.4,20,20,.25,.5,.005,400);
%[err5_500,aff5_500,stab5_500] = comp_transmot(.5,20,20,.25,.5,.005,500);
%[err6_500,aff6_500,stab6_500] = comp_transmot(.6,20,20,.25,.5,.005,500);
%[err7_500,aff7_500,stab7_500] = comp_transmot(.7,20,20,.25,.5,.005,500);

box_frame;
boxframe = boxframe./80;

%[r_err2_100,r_aff2_100,r_stab2_100] = comp_transmot(.2,8,40,0,0,0,100,boxframe);
%[r_err3_100,r_aff3_100,r_stab3_100] = comp_transmot(.3,8,40,0,0,0,100,boxframe);
%[r_err4_100,r_aff4_100,r_stab4_100] = comp_transmot(.4,8,40,0,0,0,100,boxframe);
%[r_err5_100,r_aff5_100,r_stab5_100] = comp_transmot(.5,8,40,0,0,0,100,boxframe);
%[r_err6_100,r_aff6_100,r_stab6_100] = comp_transmot(.6,8,40,0,0,0,100,boxframe);


[r_err2_100b,r_aff2_100b,r_stab2_100b] = comp_transmot(.2,8,40,0,0,0,100,boxframe);
r_err2_200 = [r_err2_100; r_err2_100b];
r_aff2_200 = [r_aff2_100; r_aff2_100b];
r_stab2_200 = [r_stab2_100; r_stab2_100b];
[r_err3_100b,r_aff3_100b,r_stab3_100b] = comp_transmot(.3,8,40,0,0,0,100,boxframe);
r_err3_200 = [r_err3_100; r_err3_100b];
r_aff3_200 = [r_aff3_100; r_aff3_100b];
r_stab3_200 = [r_stab3_100; r_stab3_100b];
[r_err4_100b,r_aff4_100b,r_stab4_100b] = comp_transmot(.4,8,40,0,0,0,100,boxframe);
r_err4_200 = [r_err4_100; r_err4_100b];
r_aff4_200 = [r_aff4_100; r_aff4_100b];
r_stab4_200 = [r_stab4_100; r_stab4_100b];
[r_err5_100b,r_aff5_100b,r_stab5_100b] = comp_transmot(.5,8,40,0,0,0,100,boxframe);
r_err5_200 = [r_err5_100; r_err5_100b];
r_aff5_200 = [r_aff5_100; r_aff5_100b];
r_stab5_200 = [r_stab5_100; r_stab5_100b];
[r_err6_100b,r_aff6_100b,r_stab6_100b] = comp_transmot(.6,8,40,0,0,0,100,boxframe);
r_err6_200 = [r_err6_100; r_err6_100b];
r_aff6_200 = [r_aff6_100; r_aff6_100b];
r_stab6_200 = [r_stab6_100; r_stab6_100b];
[r_err7_100b,r_aff7_100b,r_stab7_100b] = comp_transmot(.7,8,40,0,0,0,100,boxframe);
r_err7_200 = [r_err7_100; r_err7_100b];
r_aff7_200 = [r_aff7_100; r_aff7_100b];
r_stab7_200 = [r_stab7_100; r_stab7_100b];
