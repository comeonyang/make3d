for n = [1000]
n
[e2, a2, s2] = comp_transmot(0.2,20,20,.25,.5,0.005,400,n);
save_res('transres.m', e2, a2, s2, ...
   sprintf('comp_transmot(0.2,20,20,.25,.5,0.005,400,%d);', n), .2);
[e3, a3, s3] = comp_transmot(0.3,20,20,.25,.5,0.005,400,n);
save_res('transres.m', e3, a3, s3, ...
   sprintf('comp_transmot(0.3,20,20,.25,.5,0.005,400,%d);', n), .3);
[e4, a4, s4] = comp_transmot(0.4,20,20,.25,.5,0.005,400,n);
save_res('transres.m', e4, a4, s4, ...
   sprintf('comp_transmot(0.4,20,20,.25,.5,0.005,400,%d);', n), .4);
[e5, a5, s5] = comp_transmot(0.5,20,20,.25,.5,0.005,400,n);
save_res('transres.m', e5, a5, s5, ...
   sprintf('comp_transmot(0.5,20,20,.25,.5,0.005,400,%d);', n), .5);
[e6, a6, s6] = comp_transmot(0.6,20,20,.25,.5,0.005,400,n);
save_res('transres.m', e6, a6, s6, ...
   sprintf('comp_transmot(0.6,20,20,.25,.5,0.005,400,%d);', n), .6);
[e7, a7, s7] = comp_transmot(0.7,20,20,.25,.5,0.005,400,n);
save_res('transres.m', e7, a7, s7, ...
   sprintf('comp_transmot(0.7,20,20,.25,.5,0.005,400,%d);', n), .7);
end
