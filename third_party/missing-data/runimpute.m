frot = .5;
ftrans = .25;
sigma = .005;

for n = [100]
n
%[e9,a9,s9] = impute_compare2(frot,ftrans,sigma,9,n);
%save_impute('imputeres.m', e9, a9, s9, ...
%   sprintf('impute_compare2(%d,%d,%d,9,%d);', frot, ftrans, sigma, n), 9);
[e11,a11,s11] = impute_compare2(frot,ftrans,sigma,11,n);
save_impute('imputeres.m', e11, a11, s11, ...
   sprintf('impute_compare2(%d,%d,%d,11,%d);', frot, ftrans, sigma, n), 11);
%[e13,a13,s13] = impute_compare2(frot,ftrans,sigma,13,n);
%save_impute('imputeres.m', e13, a13, s13, ...
%   sprintf('impute_compare2(%d,%d,%d,13,%d);', frot, ftrans, sigma, n), 13);
[e15,a15,s15] = impute_compare2(frot,ftrans,sigma,15,n);
save_impute('imputeres.m', e15, a15, s15, ...
   sprintf('impute_compare2(%d,%d,%d,15,%d);', frot, ftrans, sigma, n), 15);
%[e17,a17,s17] = impute_compare2(frot,ftrans,sigma,17,n);
%save_impute('imputeres.m', e17, a17, s17, ...
%   sprintf('impute_compare2(%d,%d,%d,17,%d);', frot, ftrans, sigma, n), 17);
%[e19,a19,s19] = impute_compare2(frot,ftrans,sigma,19,n);
%save_impute('imputeres.m', e19, a19, s19, ...
%   sprintf('impute_compare2(%d,%d,%d,19,%d);', frot, ftrans, sigma, n), 19);
end
